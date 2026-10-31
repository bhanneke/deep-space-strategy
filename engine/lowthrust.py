"""Low-thrust transfers with continuous thrust and changing mass (pykep Sims-Flanagan, independently checked).

A low-thrust journey has three phases, kept apart in the results:
  1. planet-side departure from a parking orbit to escape (conventions.planet_side_plan). When the escape burn lasts at
     most a tenth of the parking-orbit period, the ship makes one Oberth burn of sqrt(2 mu/r) - sqrt(mu/r); longer burns
     split over up to 50 perigee passes with the same velocity change. Beyond that the ship flies a slow tangential spiral, whose velocity change from a circular orbit to escape equals the circular
     speed (Edelbaum); the phase's time and fuel follow from thrust and mass flow. At small bodies the ship starts at rest relative to
     the body. Optionally a launcher or kick stage supplies a departure excess speed (vinf_dep), not paid by the ship.
  2. heliocentric leg, in two steps. Sims-Flanagan (pykep.trajopt.sf_pl2pl) finds the route quickly from several
     random starts; it replaces each segment's thrust with one impulse, which is coarse for long segments. Its throttle
     profile and dates then seed pykep's zero-order-hold model (pykep.trajopt.zoh_pl2pl): thrust constant in magnitude
     and direction within each segment and continuous in time, integrated with the heyoka Taylor integrator together
     with the falling mass (dm/dt = -|T|/ve). SLSQP maximises the arrival mass subject to position, velocity and mass
     continuity. The zero-order-hold solution is the model's result.
  3. planet-side capture: the reverse burn or spiral into the destination's parking orbit (none at small bodies).

Which ships fly this model: those below 4 mm/s2 at full mass (Ship.low_thrust). Stronger ships burn for hours to days
against flights of months and fly short-burn legs (journey.py); segments of weeks cannot describe their burns.

Independent check: the solution's thrust profile is integrated as truly continuous thrust (two-body gravity plus
thrust, dm/dt = -|T|/ve) with scipy's DOP853 integrator, separate from pykep's propagation, and compared with the
arrival planet's position and velocity and with the optimiser's final mass.
"""
import math, sys
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work')
import numpy as np

from engine.conventions import end_state, NAMES, planet_side_plan, pass_days
from engine.ship import G0

MU_SUN = 1.32712440018e20            # m^3/s^2
# The zero-order-hold problem needs a non-zero range for the excess speeds at both ends (a zero-width bound leaves their
# direction variables without effect and stalls the solvers). 1e-4 in units of Earth's orbital speed allows up to about
# 3 m/s of numerical slack; the result reports the excess speeds the solution uses.
VINF_SLACK = 1e-4
# a refined leg counts only when the independent integration confirms it: arrival within 100 km of the target planet and
# the final mass within 1 kg (a feasible-looking optimiser solution once missed Mars by 4 million km in the check)
CHECK_KM, CHECK_KG = 100.0, 1.0
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

def lambert_seeds(origin, target, ship, m0, t0_bounds, tof_bounds, nseg, count=4):
    """Sims-Flanagan starting points from two-impulse (Lambert) transfers with up to two revolutions: the cheapest few
    departure dates and flight times on a 24 x 24 grid inside the bounds, each turned into a thrust profile that applies the departure impulse at
    full throttle from the first segment on and the arrival impulse from the last segment backwards. Strong engines fly
    close to these profiles; random starts alone end in local optima above about 1 mm/s2."""
    import pykep as pk, bodies
    atlas = bodies.load_atlas()
    pls, plf = bodies.body(origin, atlas), bodies.body(target, atlas)
    cand = []
    for t0 in np.linspace(t0_bounds[0], t0_bounds[1], 24):
        r0, v0p = pls.eph(float(t0))
        for tof in np.linspace(tof_bounds[0], tof_bounds[1], 24):
            r1, v1p = plf.eph(float(t0 + tof))
            try:
                l = pk.lambert_problem(r0=r0, r1=r1, tof=float(tof) * DAY, mu=MU_SUN, multi_revs=2)
            except Exception:
                continue
            for v0, v1 in zip(l.v0, l.v1):                     # zero- and multi-revolution solutions
                d1 = np.array(v0) - np.array(v0p); d2 = np.array(v1p) - np.array(v1)
                cand.append((float(np.linalg.norm(d1) + np.linalg.norm(d2)), float(t0), float(tof), d1, d2))
    cand.sort(key=lambda c: c[0])
    seeds, used = [], []
    for dv, t0, tof, d1, d2 in cand:
        if any(abs(t0 - a) < 0.1 * (t0_bounds[1] - t0_bounds[0] + 1) and abs(tof - b) < 0.1 * (tof_bounds[1] - tof_bounds[0] + 1) for a, b in used):
            continue                                              # distinct seeds only
        used.append((t0, tof))
        cap = ship.thrust_n * tof * DAY / nseg / m0              # velocity change of one segment at full throttle, m/s
        thr = np.zeros((nseg, 3))
        for dvec, order in ((d1, range(nseg)), (d2, range(nseg - 1, -1, -1))):
            left = float(np.linalg.norm(dvec)); u = dvec / max(left, 1e-12)
            for i in order:
                if left <= 0:
                    break
                take = min(left, cap * (1 - np.linalg.norm(thr[i])))
                thr[i] += u * take / cap; left -= take
        mf = m0 * math.exp(-dv / ship.ve_kms / 1000)
        seeds.append([t0, max(mf, ship.dry_kg + ship.payload_kg)] + [0.0] * 6 + list(thr.ravel()) + [tof])
        if len(seeds) == count:
            break
    return seeds

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
    lb, ub = prob.get_bounds()
    firsts = [list(np.clip(x, lb, ub)) for x in lambert_seeds(origin, target, ship, m0, t0_bounds, tof_bounds, nseg)] \
        if vinf_dep == 0 and vinf_arr == 0 else []
    for k in range(len(firsts) + starts):                           # Lambert seeds first, then scaled random starts
        pop = pg.population(prob, 1, seed=seed + k)
        if k < len(firsts):
            pop.set_x(0, firsts[k])
        else:
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
    prob.c_tol = [1e-8] * prob.get_nc()              # 1e-8 AU = 1.5 km; the independent check (CHECK_KM) is the real gate
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
    # pykep applies throttle x max thrust x direction, with the direction vector allowed below unit length
    # (inequalities_for_tc); normalising it here added thrust and put the check off by up to 47 km at 5 kN
    thrust = [list(r[0] * ship.thrust_n * np.asarray(r[1:])) for r in seg]
    vdep = list(z[2] * V * np.asarray(z[3:6])); varr = list(z[6] * V * np.asarray(z[7:10]))
    leg = {'origin': origin, 'target': target, 'model': 'zero-order hold (continuous thrust)', 't0_mjd2000': t0, 'tof_days': tof,
           'm0_kg': m0, 'mf_kg': mf, 'fuel_kg': m0 - mf, 'nseg': nseg, 'thrust_n': thrust,
           'max_throttle': float(seg[:, 0].max()), 'dv_equivalent_kms': ship.ve_kms * math.log(m0 / mf),
           'vinf_dep_vec_ms': vdep, 'vinf_arr_vec_ms': varr,
           'vinf_dep_ms': float(np.linalg.norm(vdep)), 'vinf_arr_ms': float(np.linalg.norm(varr))}
    leg['check'] = independent_check(t0, tof, vdep, varr, thrust, mf, ship, m0, pls, plf)
    c = leg['check']
    if c['position_error_km'] > CHECK_KM or abs(c['mass_error_kg']) > CHECK_KG:
        return None                                   # the optimiser called it feasible; the independent integration disagrees
    return leg

A_REF = 0.25e-3          # m/s2: reference acceleration for continuation (the search is reliable here)

def heliocentric(origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep=0.0, vinf_arr=0.0, nseg=24, starts=8, seed=1):
    """best continuous-thrust leg: Sims-Flanagan search, then zero-order-hold refinement. None if infeasible.

    Continuation for stronger engines: a ship with more thrust can fly exactly what a weaker ship flies (the same force
    gives the same trajectory and, at the same exhaust speed, the same mass flow). Above 2 x A_REF the engine also solves
    a reference ship with A_REF and starts a refinement from its solution, throttles scaled down. That start is feasible,
    so the result is at least as good as the weaker ship's; random starts alone ended in local optima above ~1 mm/s2."""
    import dataclasses
    best = None
    def take(leg, label):
        nonlocal best
        if leg is not None and (best is None or leg['mf_kg'] > best['mf_kg']):
            leg['found_by'] = label; best = leg
    sfs = sims_flanagan(origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep, vinf_arr, nseg, starts, seed)
    for sf in sfs:                                    # every feasible Sims-Flanagan solution seeds a refinement
        leg = zero_order_hold(sf, origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep, vinf_arr)
        if leg is not None:
            leg['seed_sims_flanagan_fuel_kg'] = sf['fuel_kg']
        take(leg, 'search')
    if ship.thrust_n / m0 > 2 * A_REF:
        ref = dataclasses.replace(ship, thrust_n=A_REF * m0)
        r = heliocentric(origin, target, ref, m0, t0_bounds, tof_bounds, vinf_dep, vinf_arr, nseg, starts, seed)
        if r is not None:
            x = ([r['t0_mjd2000'], r['mf_kg']] + list(r['vinf_dep_vec_ms']) + list(r['vinf_arr_vec_ms'])
                 + [c / ship.thrust_n for T in r['thrust_n'] for c in T] + [r['tof_days']])
            leg = zero_order_hold({'x': x, 'nseg': nseg}, origin, target, ship, m0, t0_bounds, tof_bounds, vinf_dep, vinf_arr)
            if leg is not None:
                leg['continuation_reference_dv_kms'] = r['dv_equivalent_kms']
            take(leg, f'continuation from {1000 * A_REF:g} mm/s2')
    if best is not None and sfs:
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
    """between the node's parking orbit and escape (zero excess speed), as conventions.planet_side_plan decides:
    one Oberth burn, the same burn split over perigee passes (plus the time on the intermediate orbits), or the
    Edelbaum spiral. Returns (mode, dv km/s, days, fuel kg, mass after kg); mode 'burn', 'passes', 'spiral' or 'none'."""
    e = end_state(node)
    if e['kind'] != 'parking':
        return 'none', 0.0, 0.0, 0.0, mass_kg
    mode, n = planet_side_plan(node, mass_kg, ship.ve_kms, ship.mdot_kgs)
    if mode == 'spiral':
        return ('spiral',) + spiral(node, ship, mass_kg)
    dv = math.sqrt(2 * e['mu'] / e['r']) - math.sqrt(e['mu'] / e['r'])
    m1 = mass_kg * math.exp(-dv / ship.ve_kms)
    days = (mass_kg - m1) / ship.mdot_kgs / DAY + pass_days(node, n)
    return mode, dv, days, mass_kg - m1, m1

def segments(ship, m0):
    """24 segments for weak engines; 48 from 1 mm/s2 of initial acceleration, where burns are shorter"""
    return 48 if ship.thrust_n / m0 >= 1e-3 else 24

def transfer(origin, target, ship, t0_bounds, tof_bounds, vinf_dep=0.0, nseg=None, starts=8, seed=1):
    """full low-thrust transfer: departure spiral, heliocentric leg, capture spiral; mass carried through all three"""
    m0 = ship.mass(ship.fuel_capacity_kg)
    nseg = nseg or segments(ship, m0)
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
