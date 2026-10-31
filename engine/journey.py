"""Chronological multi-leg journeys for a ship with fuel, payload and refuelling.

A journey is a sequence of legs. Each leg leaves a node on a network date no earlier than the ship's arrival there,
flies one of the precomputed transfer alternatives, and costs two impulsive burns: the departure burn (applied to the
ship's mass at departure) and the arrival burn (applied to the lighter ship at arrival). Fuel, payload and dry mass
are carried through the whole journey. At nodes that supply propellant (and that the query allows) the ship refills to
capacity; the journey records how much fuel it loaded in total.

The search is a label-setting search over the time-expanded graph. A label holds node, arrival date, fuel on board and
fuel loaded so far; a label is dropped when another label at the same node arrived no later, carries at least as much
fuel and loaded no more. The result is the set of non-dominated arrivals at the destination (earliest, least fuel, and
the trade-offs between them), each with its full itinerary and mass history.

Nodes: the 12 network sites and Earth-Moon L1 ('eml1'). L1 shares Earth's heliocentric legs and prices its departure and
arrival with the L1 model; L1 and low Earth orbit connect directly for 3.77 km/s each way in 4 days.

Flight model per ship and transfer (conventions.planet_side_plan). At each planet the ship makes one Oberth burn when
that burn lasts at most a tenth of the parking-orbit period; longer burns split into perigee passes with the same
velocity change, plus the time on the intermediate orbits; beyond 50 passes the ship spirals between the parking orbit
and escape and then makes an impulse of v-infinity at the edge of the sphere of influence. Pass and spiral time count
in the chronology: the departure phase starts after the ship's arrival and inside the departure window, the arrival
phase adds to the flight.
A leg counts only when its heliocentric impulses (v-infinity impulses after a spiral, and rendezvous burns at small
bodies) take at most IMPULSIVE_SHARE of the flight time; longer burns need the continuous-thrust model (lowthrust.py).
"""
import heapq, math
import numpy as np
from engine.conventions import SITES, PLACES, FUEL_SOURCES, NAMES, planet_side_mode, planet_side_arr, planet_side_plan, pass_days, end_state, parking_period_s, helio, local_transfer
from engine.ship import burn

IMPULSIVE_SHARE = 0.05      # heliocentric impulses may last at most this share of the flight time
DAY = 86400.0

def node_mode(node, ship):
    """(mode, passes) for the ship at the node, decided at full mass (the heaviest, so the longest burns)"""
    return planet_side_plan(node, ship.mass(ship.fuel_capacity_kg), ship.ve_kms, ship.mdot_kgs)

def helio_end(node, mode):
    """does the burn at this end happen in heliocentric space (after a spiral, or at a small body)?"""
    return mode == 'spiral' and node != 'eml1' or end_state(node)['kind'] == 'rendezvous'

NODES = SITES + ['eml1']          # the 12-site default; plan() uses every place the loaded data covers

def nodes_for(alts):
    """places whose heliocentric body appears in the loaded alternatives"""
    return [p for p in PLACES if helio(p) in alts.bodies]

class Network:
    """priced transfer alternatives for every ordered node pair and planet-side modes: arrays [date, alternative] of
    flight time, the two planet-side costs and their spiral parts"""
    def __init__(self, alts):
        self.alts, self._cache = alts, {}

    def edge(self, a, b, ma=('burn', 1), mb=('burn', 1)):
        key = (a, b, ma, mb)
        if key in self._cache:
            return self._cache[key]
        ha, hb = helio(a), helio(b)
        if ha == hb:                                       # inside one planet system (L1, the Moon, moons), or no edge
            res = None
            lt = local_transfer(a, b)
            if lt is not None and 'spiral' not in (ma[0], mb[0]):   # a spiralling ship has no short local transfers
                n = self.alts.n_dates
                z = np.zeros((n, 1))
                res = (np.full((n, 1), lt[2], float), np.full((n, 1), lt[0]), np.full((n, 1), lt[1]), z, z, z, z)
        else:
            i = self.alts.index.get(f'{ha}>{hb}')
            if i is None:
                res = None
            else:
                arr = self.alts.route(i).astype(float)     # [date, k, (tof, vinf_dep, vinf_arr)]
                vdep, varr = arr[:, :, 1], arr[:, :, 2]
                d1, s1 = planet_side_arr(a, vdep, ma[0]); d2, s2 = planet_side_arr(b, varr, mb[0])
                t1 = np.full_like(d1, pass_days(a, ma[1]) if ma[0] == 'passes' else 0.0)
                t2 = np.full_like(d2, pass_days(b, mb[1]) if mb[0] == 'passes' else 0.0)
                res = (arr[:, :, 0], d1, d2, s1, s2, t1, t2)
        self._cache[key] = res
        return res

class Label:
    __slots__ = ('node', 'k', 'fuel', 'loaded', 'legs', 'parent')
    def __init__(self, node, k, fuel, loaded, legs):
        self.node, self.k, self.fuel, self.loaded, self.legs = node, k, fuel, loaded, legs
    def __lt__(self, other):
        return (self.k, -self.fuel) < (other.k, -other.fuel)

def _dominated(front, k, fuel, loaded):
    return any(f[0] <= k and f[1] >= fuel - 1e-6 and f[2] <= loaded + 1e-6 for f in front)

def plan(alts, ship, origin, destination, depart_from_jd, depart_until_jd, arrive_by_jd=None, refuel_at=None,
         max_legs=4, max_wait_days=1100, start_fuel_kg=None, net=None, max_results=8):
    """non-dominated journeys from origin to destination for an impulsive ship.

    depart_from_jd / depart_until_jd: window for the first departure; arrive_by_jd: deadline (default: 15 years after
    the window opens). refuel_at: nodes where the ship may refill (default: all fuel sources). start_fuel_kg: fuel at
    the origin (default: full tanks).
    """
    net = net or Network(alts)
    nodes = nodes_for(alts)
    mode = {n: node_mode(n, ship) for n in nodes}
    mdot = ship.mdot_kgs
    step = alts.step
    # the departure window covers the network dates inside it: the first on or after its start, the last on or before its end
    k0 = max(0, math.ceil((depart_from_jd - alts.jd0) / step - 1e-9))
    k_dep_end = min(alts.n_dates - 1, math.floor((depart_until_jd - alts.jd0) / step + 1e-9))
    x_dep_start = (depart_from_jd - alts.jd0) / step             # earliest parking-orbit departure, in steps
    x_dep_end = (depart_until_jd - alts.jd0) / step              # latest parking-orbit departure, in steps
    # the search runs on network dates; the deadline itself is checked against each leg's actual arrival time
    arrive_by_jd = arrive_by_jd or alts.jd(k0) + 15 * 365.25
    k_end = min(alts.n_dates - 1, alts.k_of(arrive_by_jd) + 1)
    refuel = set(FUEL_SOURCES if refuel_at is None else refuel_at)
    wait = int(max_wait_days / step)
    fuel0 = ship.fuel_capacity_kg if start_fuel_kg is None else start_fuel_kg
    start = Label(origin, k0, fuel0, fuel0 if origin in refuel else 0.0, [])
    fronts = {n: [] for n in nodes}
    heap = [start]; done = []
    ve, dry, pay = ship.ve_kms, ship.dry_kg, ship.payload_kg
    while heap:
        L = heapq.heappop(heap)
        if L.node == destination:
            continue
        if len(L.legs) >= max_legs:
            continue
        fuel, loaded = L.fuel, L.loaded
        if L.node in refuel and L.legs:                  # refill at a propellant source (the origin's fill counts at start)
            loaded += ship.fuel_capacity_kg - fuel; fuel = ship.fuel_capacity_kg
        last_dep = x_dep_end if not L.legs else min(k_end, L.k + wait)   # latest parking-orbit departure, in steps
        if L.k > last_dep:
            continue
        for b in nodes:
            if b == L.node or (b != destination and b not in refuel):   # stops only where the ship can refuel
                continue
            e = net.edge(L.node, b, mode[L.node], mode[b])
            if e is None:
                continue
            # mass accounting for the four stages on every date and alternative at once: departure spiral, departure
            # impulse, arrival impulse, arrival spiral (spirals are zero for ships that burn at the planet)
            m0 = dry + pay + fuel
            s1 = float(np.nanmax(e[3])) if np.isfinite(e[3]).any() else 0.0      # the spiral part is the same on every date
            days1s = m0 * (1 - math.exp(-s1 / ve)) / mdot / DAY
            hi = min(alts.n_dates - 1, int(math.floor(last_dep + (days1s + float(np.nanmax(e[5]))) / step + 1e-9)))
            tof, dv1, dv2, _, s2, t1, t2 = (x[L.k:hi + 1] for x in e)
            if tof.size == 0:
                continue
            ma = m0 * math.exp(-s1 / ve)
            mb = ma * np.exp(-(dv1 - s1) / ve)
            mc = mb * np.exp(-(dv2 - s2) / ve)
            m2 = mc * np.exp(-s2 / ve)
            fuel_after = fuel - (m0 - m2)
            days1 = days1s + t1                                          # departure spiral or perigee passes
            days2 = (mc - m2) / mdot / DAY + t2
            helio_kg = (ma - mb) * helio_end(L.node, mode[L.node][0]) + (mb - mc) * helio_end(b, mode[b][0])
            kd = np.arange(L.k, hi + 1)[:, None] + np.zeros_like(tof)
            x_park = kd - days1 / step                                   # parking-orbit departure, in steps
            k2 = kd + np.ceil((tof + days2) / step)
            arr_jd = alts.jd0 + kd * step + tof + days2
            ok = (~np.isnan(tof) & (fuel_after >= 0) & (k2 <= k_end) & (arr_jd <= arrive_by_jd)
                  & (x_park >= (x_dep_start if not L.legs else L.k) - 1e-9) & (x_park <= last_dep + 1e-9) & (helio_kg / mdot <= IMPULSIVE_SHARE * tof * DAY))
            if not ok.any():
                continue
            cand = sorted(zip(k2[ok], -fuel_after[ok], kd[ok], np.argwhere(ok)[:, 1]))
            best_fuel = -math.inf
            for k2v, nf, kdv, alt in cand:                   # earliest first; keep only those that also gain fuel
                fa = -nf
                if fa <= best_fuel + 1e-6:
                    continue
                best_fuel = fa
                k2v, kdv, alt = int(k2v), int(kdv), int(alt)
                if _dominated(fronts[b], k2v, fa, loaded):
                    continue
                fronts[b] = [f for f in fronts[b] if not (k2v <= f[0] and fa >= f[1] - 1e-6 and loaded <= f[2] + 1e-6)]
                fronts[b].append((k2v, fa, loaded))
                leg = {'from': L.node, 'to': b, 'k_dep': kdv, 'k_arr': k2v, 'alt': alt, 'modes': (mode[L.node][0], mode[b][0]),
                       'passes': (mode[L.node][1], mode[b][1]), 'pass_days': (float(e[5][kdv, alt]), float(e[6][kdv, alt])),
                       'tof_days': float(e[0][kdv, alt]), 'dv_dep': float(e[1][kdv, alt]), 'dv_arr': float(e[2][kdv, alt]),
                       'sp_dep': float(e[3][kdv, alt]), 'sp_arr': float(e[4][kdv, alt]),
                       'fuel_before': fuel, 'fuel_after': fa, 'refuelled_before': L.node in refuel and bool(L.legs)}
                nl = Label(b, k2v, fa, loaded, L.legs + [leg])
                if b == destination:
                    done.append(nl)
                heapq.heappush(heap, nl)
    # non-dominated arrivals at the destination
    # three criteria: arrive earlier, load less fuel, arrive with more fuel left
    res = []
    for L in sorted(done, key=lambda x: (x.k, x.loaded, -x.fuel)):
        if any(r.k <= L.k and r.loaded <= L.loaded + 1e-6 and r.fuel >= L.fuel - 1e-6 for r in res):
            continue
        res.append(L)
    return sorted((describe(alts, ship, L) for L in res[:max_results]), key=lambda j: j['arrive_jd'])

def stages(ship, fuel, g, payload_kg=None):
    """the leg's four stages with Ship.burn, in order: departure spiral, departure impulse, arrival impulse, arrival spiral.
    Returns [(name, BurnResult)]; zero stages are left out."""
    park = lambda n: parking_period_s('earth' if n == 'eml1' else n)
    def name(end, mode, n):
        return {'burn': f'{end} burn', 'passes': f'{end} burns over {n} perigee passes',
                'spiral': f'{end} impulse at the sphere of influence'}[mode]
    n = g.get('passes', (1, 1))
    plan_ = [('departure spiral', g['sp_dep'], None),
             (name('departure', g['modes'][0], n[0]), g['dv_dep'] - g['sp_dep'], park(g['from']) if g['modes'][0] == 'burn' else None),
             (name('arrival', g['modes'][1], n[1]), g['dv_arr'] - g['sp_arr'], park(g['to']) if g['modes'][1] == 'burn' else None),
             ('arrival spiral', g['sp_arr'], None)]
    out = []
    for name, dv, period in plan_:
        if dv <= 0:
            continue
        b = burn(ship, fuel, dv, payload_kg, parking_period_s=period)
        out.append((name, b)); fuel = b.fuel_after_kg
    return out

def helio_seconds(g, st):
    """seconds of heliocentric impulses in a leg (the burns IMPULSIVE_SHARE limits)"""
    sec = 0.0
    for name, b in st:
        if 'spiral' in name:
            continue
        end = g['from'] if name.startswith('departure') else g['to']
        if helio_end(end, g['modes'][0] if name.startswith('departure') else g['modes'][1]):
            sec += b.duration_s
    return sec

def describe(alts, ship, L):
    """itinerary with dates, burns and the mass history. Each stage is recomputed with Ship.burn as a check of the
    search's vectorised mass accounting, and fuel must carry over from leg to leg (or be refilled to capacity)."""
    legs, prev_after = [], None
    for g in L.legs:
        f0 = g['fuel_before']
        if prev_after is not None:
            expected = ship.fuel_capacity_kg if g['refuelled_before'] else prev_after
            assert abs(f0 - expected) < 1e-6 * ship.fuel_capacity_kg + 1e-6, 'fuel does not carry over between legs'
        st = stages(ship, f0, g)
        assert all(b.feasible for _, b in st), 'infeasible burn in a returned journey'
        after = st[-1][1].fuel_after_kg if st else f0
        assert abs(after - g['fuel_after']) < 1e-6 * ship.fuel_capacity_kg + 1e-6, 'mass accounting mismatch'
        assert helio_seconds(g, st) <= IMPULSIVE_SHARE * g['tof_days'] * DAY * (1 + 1e-9), 'heliocentric burns too long'
        d_dep = sum(b.duration_s for n, b in st if n == 'departure spiral') / DAY + g['pass_days'][0]
        d_arr = sum(b.duration_s for n, b in st if n == 'arrival spiral') / DAY + g['pass_days'][1]
        helio_jd = alts.jd(g['k_dep'])
        legs.append({'from': g['from'], 'to': g['to'], 'from_name': NAMES[g['from']], 'to_name': NAMES[g['to']],
                     'depart_jd': helio_jd - d_dep, 'arrive_jd': helio_jd + g['tof_days'] + d_arr,
                     'tof_days': round(g['tof_days'] + d_dep + d_arr, 1), 'heliocentric_days': round(g['tof_days'], 1),
                     'dv_departure': round(g['dv_dep'], 3), 'dv_arrival': round(g['dv_arr'], 3),
                     'spiral_departure_kms': round(g['sp_dep'], 3), 'spiral_arrival_kms': round(g['sp_arr'], 3), 'modes': list(g['modes']),
                     'passes': list(g['passes']), 'pass_days': [round(x, 2) for x in g['pass_days']],
                     'stages': [{'stage': n, 'dv_kms': round(b.dv_kms, 3), 'days': round(b.duration_s / DAY, 3),
                                 'fuel_kg': round(b.fuel_used_kg, 1)} for n, b in st],
                     'heliocentric_burn_share': round(helio_seconds(g, st) / (g['tof_days'] * DAY), 4) if g['tof_days'] else 0.0,
                     'fuel_before_kg': round(f0, 1), 'fuel_after_kg': round(after, 1), 'refuelled_before': g['refuelled_before'],
                     'warnings': [w for _, b in st for w in b.warnings]})
        prev_after = after
    first = legs[0]['depart_jd']; last = legs[-1]['arrive_jd']
    return {'arrive_jd': last, 'depart_jd': first, 'days': round(last - first, 1), 'legs': legs,
            'dv_total': round(sum(l['dv_departure'] + l['dv_arrival'] for l in legs), 3),
            'fuel_loaded_kg': round(L.loaded, 1), 'fuel_left_kg': round(prev_after, 1),
            'payload_delivered_kg': ship.payload_kg, 'stops': [l['to'] for l in legs[:-1]]}

def max_payload(ship, journey, hi_kg=None, tol_kg=1.0):
    """largest payload the ship can carry through the journey's legs and refuel stops (start with full tanks), with
    heliocentric impulses within IMPULSIVE_SHARE of the flight. Same stages as the search; bisection on the payload.
    Spirals of a heavier ship take longer; the dates in the journey refer to its own payload."""
    import dataclasses
    def feasible(pay):
        s = dataclasses.replace(ship, payload_kg=pay)
        fuel = s.fuel_capacity_kg
        for leg in journey['legs']:
            if leg['refuelled_before']:
                fuel = s.fuel_capacity_kg
            g = {'from': leg['from'], 'to': leg['to'], 'modes': tuple(leg['modes']), 'passes': tuple(leg['passes']), 'dv_dep': leg['dv_departure'],
                 'dv_arr': leg['dv_arrival'], 'sp_dep': leg['spiral_departure_kms'], 'sp_arr': leg['spiral_arrival_kms']}
            st = stages(s, fuel, g)
            if not all(b.feasible for _, b in st):
                return False
            if helio_seconds(g, st) > IMPULSIVE_SHARE * leg['heliocentric_days'] * DAY:
                return False
            fuel = st[-1][1].fuel_after_kg if st else fuel
        return True
    lo, hi = 0.0, hi_kg or 50 * (ship.dry_kg + ship.fuel_capacity_kg)
    if not feasible(lo):
        return 0.0
    while hi - lo > tol_kg:
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if feasible(mid) else (lo, mid)
    return lo
