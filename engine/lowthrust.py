"""Low-thrust transfers with continuous thrust and changing mass (pykep Sims-Flanagan, independently checked).

A low-thrust journey has three phases, kept apart in the results:
  1. planet-side departure from a parking orbit to escape. When the escape burn lasts at most a tenth of the
     parking-orbit period (FINITE_BURN_LIMIT), the ship makes one Oberth burn of sqrt(2 mu/r) - sqrt(mu/r). Longer
     burns become a slow tangential spiral, whose velocity change from a circular orbit to escape equals the circular
     speed (Edelbaum); the phase's time and fuel follow from thrust and mass flow. At small bodies the ship starts at rest relative to
     the body. Optionally a launcher or kick stage supplies a departure excess speed (vinf_dep), not paid by the ship.
  2. heliocentric leg, in two steps. Sims-Flanagan (pykep.trajopt.sf_pl2pl) finds the route quickly from several
     random starts; it replaces each segment's thrust with one impulse, which is coarse for long segments. Its throttle
     profile and dates then seed pykep's zero-order-hold model (pykep.trajopt.zoh_pl2pl): thrust constant in magnitude
     and direction within each segment and continuous in time, integrated with the heyoka Taylor integrator together
     with the falling mass (dm/dt = -|T|/ve). SLSQP maximises the arrival mass subject to position, velocity and mass
     continuity. The zero-order-hold solution is the model's result.
  3. planet-side capture: the reverse burn or spiral into the destination's parking orbit (none at small bodies).

Which ships fly this model: those that need more than a day to burn a full tank (Ship.low_thrust). The model holds
for any thrust level, so high-thrust, high-Isp engines use it as well.

Independent check: the solution's thrust profile is integrated as truly continuous thrust (two-body gravity plus
thrust, dm/dt = -|T|/ve) with scipy's DOP853 integrator, separate from pykep's propagation, and compared with the
arrival planet's position and velocity and with the optimiser's final mass.
"""
import math, sys
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work')
import numpy as np

from engine.conventions import end_state, NAMES, parking_period_s
from engine.ship import G0, FINITE_BURN_LIMIT

MU_SUN = 1.32712440018e20            # m^3/s^2
# The zero-order-hold problem needs a non-zero range for the excess speeds at both ends (a zero-width bound leaves their
# direction variables without effect and stalls the solvers). 1e-4 in units of Earth's orbital speed allows up to about
# 3 m/s of numerical slack; the result reports the excess speeds the solution uses.
VINF_SLACK = 1e-4
AU = 149597870700.0
DAY = 86400.0

def spiral(node, ship, mass_kg):
    """Edelbaum spiral between the node's parking orbit and escape: (dv km/s, days, fuel kg, mass after kg)"""
    e = end_state(node)
    if e['kind'] != 'parking':
        return 0.0, 0.0, 0.0, mass_kg
    dv = math.sqrt(e['mu'] / e['r'])                                # km/s, circular speed of the parking orbit
    m1 = mass_kg * math.exp(-dv / ship.ve_kms)
    fuel = mass_kg - m1
    return dv, fuel / ship.mdot_kgs / DAY, fuel, m1

def _problem(origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep, vinf_arr, nseg):
    import pykep as pk, bodies
    atlas = bodies.load_atlas()
    pls, plf = bodies.body(origin, atlas), bodies.body(target, atlas)
    udp = pk.trajopt.sf_pl2pl(pls=pls, plf=plf, ms=m0, mu=MU_SUN, max_thrust=ship.thrust_n, veff=ship.isp_s * G0,
                              t0_bounds=list(t0_bounds), tof_bounds=list(tof_bounds), mf_bounds=[ship.dry_kg + ship.payload_kg, m0],
                              vinfs=vinf_dep, vinff=vinf_arr, nseg=nseg, cut=0.6, mass_scaling=m0, with_gradient=True)
    return udp, pls, plf

def sims_flanagan(origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep=0.0, vinf_arr=0.0, nseg=24, starts=8, seed=1):
    """feasible Sims-Flanagan legs, best (highest arrival mass) first: the fast first step. Bounds in MJD2000 and days."""
    import pygmo as pg
    udp, pls, plf = _problem(origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep, vinf_arr, nseg)
    prob = pg.problem(udp)
    prob.c_tol = [1e-7] * prob.get_nc()                             # scaled units: 1e-7 AU = 15 km, 1e-7 x 30 km/s
    algo = pg.algorithm(pg.nlopt('slsqp'))
    algo.extract(pg.nlopt).maxeval = 1500
    found = []
    # random starts with throttles scaled so that the first guess spends about 10 km/s: full random throttles on a
    # high-thrust ship push it onto an escape orbit from the Sun, where the propagation returns NaN
    scale = min(1.0, 10000.0 / (ship.thrust_n / m0 * 0.5 * (tof_bounds[0] + tof_bounds[1]) * DAY))
    for k in range(starts):
        pop = pg.population(prob, 1, seed=seed + k)
        x0 = np.array(pop.get_x()[0]); x0[8:-1] *= scale; pop.set_x(0, x0)
        pop = algo.evolve(pop)
        x = pop.champion_x
        if prob.feasibility_x(x):
            found.append(list(x))
    return [_sf_result(x, origin, target, ship, m0, nseg) for x in sorted(found, key=lambda x: -x[1])]

def _sf_result(x, origin, target, ship, m0, nseg):
    thr = np.array(x[8:-1]).reshape(nseg, 3)
    return {'origin': origin, 'target': target, 't0_mjd2000': x[0], 'tof_days': x[-1], 'm0_kg': m0, 'mf_kg': x[1],
            'fuel_kg': m0 - x[1], 'vinf_dep_vec_ms': x[2:5], 'vinf_arr_vec_ms': x[5:8], 'throttles': thr.tolist(), 'nseg': nseg,
            'max_throttle': float(np.linalg.norm(thr, axis=1).max()), 'dv_equivalent_kms': ship.ve_kms * math.log(m0 / x[1]),
            'thrust_n': (thr * ship.thrust_n).tolist(), 'x': list(x)}

def zero_order_hold(sf, origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep=0.0, vinf_arr=0.0, tol=1e-12):
    """refine a Sims-Flanagan leg in pykep's continuous-thrust zero-order-hold model, seeded with its throttle profile"""
    import pykep as pk, pygmo as pg, bodies
    atlas = bodies.load_atlas()
    pls, plf = bodies.body(origin, atlas), bodies.body(target, atlas)
    L, V = AU, math.sqrt(MU_SUN / AU)                  # integrator units: AU, mu = 1, mass relative to m0
    TIME = L / V
    t_nd = ship.thrust_n * TIME ** 2 / (m0 * L)
    ve_nd = ship.isp_s * G0 / V
    ta, ta_var = pk.ta.get_zoh_kep(tol), pk.ta.get_zoh_kep_var(tol)
    ta.pars[4] = 1.0 / ve_nd; ta_var.pars[4] = 1.0 / ve_nd
    nseg = sf['nseg']
    udp = pk.trajopt.zoh_pl2pl(pls=pls, plf=plf, ms=1.0, max_thrust=t_nd, t0_bounds=list(t0_bounds),
                               tof_bounds=[tof_bounds[0] * DAY / TIME, tof_bounds[1] * DAY / TIME],
                               mf_bounds=[(ship.dry_kg + ship.payload_kg) / m0, 1.0],
                               vinf_dep_bounds=[0.0, max(vinf_dep * 1000 / V, VINF_SLACK)],
                               vinf_arr_bounds=[0.0, max(vinf_arr * 1000 / V, VINF_SLACK)],
                               nseg=nseg, cut=0.6, tas=(ta, ta_var), inequalities_for_tc=True, L=L, V=V)
    prob = pg.problem(udp)
    prob.c_tol = [1e-9] * prob.get_nc()
    # seed: the Sims-Flanagan solution in zero-order-hold variables
    x = sf['x']
    def unit(v):
        n = np.linalg.norm(v); return (list(np.asarray(v) / n) if n > 0 else [1.0, 0.0, 0.0]), n
    idep, vd = unit(x[2:5]); iarr, va = unit(x[5:8])
    ctrl = []
    for u in np.array(x[8:-1]).reshape(nseg, 3):
        d, n = unit(u); ctrl += [min(n, 1.0)] + d
    seed = [(x[0] - t0_bounds[0]) * DAY / TIME, x[1] / m0, vd / V] + idep + [va / V] + iarr + ctrl + [x[-1] * DAY / TIME]
    lb, ub = prob.get_bounds()
    seed = list(np.clip(seed, lb, ub))
    ip = pg.ipopt(); ip.set_numeric_option('tol', 1e-9); ip.set_numeric_option('constr_viol_tol', 1e-9)
    ip.set_integer_option('max_iter', 1000)
    algo = pg.algorithm(ip)                          # SLSQP stalls on this problem; IPOPT converges from the seed
    pop = pg.population(prob); pop.push_back(seed)
    pop = algo.evolve(pop)
    z = list(pop.champion_x)
    if not prob.feasibility_x(z):
        return None
    t0 = t0_bounds[0] + z[0] * TIME / DAY
    tof = z[10 + 4 * nseg] * TIME / DAY
    mf = z[1] * m0
    seg = np.array(z[10:10 + 4 * nseg]).reshape(nseg, 4)
    thrust = [list(r[0] * ship.thrust_n * np.asarray(r[1:]) / np.linalg.norm(r[1:])) for r in seg]
    vdep = list(z[2] * V * np.asarray(z[3:6])); varr = list(z[6] * V * np.asarray(z[7:10]))
    leg = {'origin': origin, 'target': target, 'model': 'zero-order hold (continuous thrust)', 't0_mjd2000': t0, 'tof_days': tof,
           'm0_kg': m0, 'mf_kg': mf, 'fuel_kg': m0 - mf, 'nseg': nseg, 'thrust_n': thrust,
           'max_throttle': float(seg[:, 0].max()), 'dv_equivalent_kms': ship.ve_kms * math.log(m0 / mf),
           'vinf_dep_vec_ms': vdep, 'vinf_arr_vec_ms': varr,
           'vinf_dep_ms': float(np.linalg.norm(vdep)), 'vinf_arr_ms': float(np.linalg.norm(varr))}
    leg['check'] = independent_check(t0, tof, vdep, varr, thrust, mf, ship, m0, pls, plf)
    return leg

def heliocentric(origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep=0.0, vinf_arr=0.0, nseg=24, starts=8, seed=1):
    """best continuous-thrust leg: Sims-Flanagan search, then zero-order-hold refinement. None if infeasible."""
    best = None
    sfs = sims_flanagan(origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep, vinf_arr, nseg, starts, seed)
    for sf in sfs:                                    # every feasible Sims-Flanagan solution seeds a refinement
        leg = zero_order_hold(sf, origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep, vinf_arr)
        if leg is not None and (best is None or leg['mf_kg'] > best['mf_kg']):
            leg['seed_sims_flanagan_fuel_kg'] = sf['fuel_kg']; best = leg
    if best is not None:
        best['best_sims_flanagan_fuel_kg'] = sfs[0]['fuel_kg']     # cross-check: the other transcription's optimum
    return best

def independent_check(t0, tof_days, vdep, varr, thrust, mf, ship, m0, pls, plf):
    """integrate the piecewise-constant thrust profile with scipy (DOP853), independent of pykep's propagators, and compare
    with the arrival planet's state and the optimiser's arrival mass"""
    from scipy.integrate import solve_ivp
    tof = tof_days * DAY
    rs, vs = pls.eph(t0)
    y = np.array(list(rs) + [a + b for a, b in zip(vs, vdep)] + [m0], float)
    ve = ship.isp_s * G0
    nseg = len(thrust)
    dt = tof / nseg
    for i in range(nseg):
        T = np.asarray(thrust[i], float)
        def f(t, s, T=T):
            r = s[:3]; rn = np.linalg.norm(r)
            return np.concatenate([s[3:6], -MU_SUN * r / rn ** 3 + T / s[6], [-np.linalg.norm(T) / ve]])
        y = solve_ivp(f, (0, dt), y, method='DOP853', rtol=1e-11, atol=1e-6).y[:, -1]
    rf, vf = plf.eph(t0 + tof_days)
    vf = [a + b for a, b in zip(vf, varr)]
    return {'position_error_km': float(np.linalg.norm(y[:3] - rf) / 1000), 'velocity_error_ms': float(np.linalg.norm(y[3:6] - vf)),
            'mass_error_kg': float(y[6] - mf), 'integrated_final_mass_kg': float(y[6])}

def planet_side(node, ship, mass_kg):
    """between the node's parking orbit and escape (zero excess speed): an impulsive Oberth burn when it lasts at most
    FINITE_BURN_LIMIT of the parking-orbit period, otherwise the Edelbaum spiral.
    Returns (mode, dv km/s, days, fuel kg, mass after kg); mode is 'burn', 'spiral' or 'none'."""
    e = end_state(node)
    if e['kind'] != 'parking':
        return 'none', 0.0, 0.0, 0.0, mass_kg
    dv = math.sqrt(2 * e['mu'] / e['r']) - math.sqrt(e['mu'] / e['r'])
    m1 = mass_kg * math.exp(-dv / ship.ve_kms)
    if (mass_kg - m1) / ship.mdot_kgs <= FINITE_BURN_LIMIT * parking_period_s(node):
        return 'burn', dv, (mass_kg - m1) / ship.mdot_kgs / DAY, mass_kg - m1, m1
    return ('spiral',) + spiral(node, ship, mass_kg)

def transfer(origin, target, ship, t0_bounds, tof_bounds, vinf_dep=0.0, nseg=24, starts=8, seed=1):
    """full low-thrust transfer: departure spiral, heliocentric leg, capture spiral; mass carried through all three"""
    m0 = ship.mass(ship.fuel_capacity_kg)
    mode1, dv1, days1, fuel1, m1 = ('kick', 0.0, 0.0, 0.0, m0) if vinf_dep > 0 else planet_side(origin, ship, m0)
    # the window bounds the departure from the parking orbit; the heliocentric leg starts after the spiral
    t0_bounds = (t0_bounds[0] + days1, t0_bounds[1] + days1)
    leg = heliocentric(origin, target, ship, m1, t0_bounds, tof_bounds, vinf_dep, 0.0, nseg, starts, seed)
    if leg is None:
        return {'feasible': False, 'reason': 'no feasible heliocentric leg in the window',
                'phases': [{'phase': f'departure {mode1} from {end_state(origin)["label"]}', 'kind': 'planet-side', 'dv_kms': dv1,
                            'days': days1, 'fuel_kg': fuel1}]}
    mode3, dv3, days3, fuel3, m3 = planet_side(target, ship, leg['mf_kg'])
    fuel_total = fuel1 + leg['fuel_kg'] + fuel3
    feasible = fuel_total <= ship.fuel_capacity_kg + 1e-6
    return {'feasible': feasible, 'reason': '' if feasible else 'fuel runs out in the capture spiral',
            'phases': [
                {'phase': f'departure {mode1} from {end_state(origin)["label"]}' if dv1 else ('kick stage supplies the departure' if vinf_dep else 'departure at rest relative to ' + NAMES[origin]),
                 'kind': 'planet-side', 'dv_kms': dv1, 'days': days1, 'fuel_kg': fuel1},
                {'phase': f'heliocentric leg {NAMES[origin]} to {NAMES[target]}', 'kind': 'heliocentric', 'dv_kms': leg['dv_equivalent_kms'],
                 'days': leg['tof_days'], 'fuel_kg': leg['fuel_kg']},
                {'phase': f'capture {mode3} into {end_state(target)["label"]}' if dv3 else 'rendezvous at ' + NAMES[target],
                 'kind': 'planet-side', 'dv_kms': dv3, 'days': days3, 'fuel_kg': fuel3}],
            'leg': leg, 'depart_mjd2000': leg['t0_mjd2000'] - days1, 'arrive_mjd2000': leg['t0_mjd2000'] + leg['tof_days'] + days3,
            'fuel_kg': fuel_total, 'payload_delivered_kg': ship.payload_kg if feasible else 0.0, 'final_mass_kg': m3}
