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
"""
import heapq, math
import numpy as np
from engine.conventions import SITES, FUEL_SOURCES, NAMES, planet_side_dv, parking_period_s, L1_LEO
from engine.ship import burn

NODES = SITES + ['eml1']
L1_DAYS = 4

class Network:
    """priced transfer alternatives for every ordered node pair: arrays [date, alternative] of flight time and burns"""
    def __init__(self, alts):
        self.alts, self._cache = alts, {}

    def edge(self, a, b):
        key = (a, b)
        if key in self._cache:
            return self._cache[key]
        ha, hb = ('earth' if a == 'eml1' else a), ('earth' if b == 'eml1' else b)
        if ha == hb:                                       # L1 <-> low Earth orbit, or no edge
            res = None
            if {a, b} == {'earth', 'eml1'}:
                n = self.alts.n_dates
                res = (np.full((n, 1), L1_DAYS, float), np.full((n, 1), L1_LEO), np.zeros((n, 1)))
        else:
            i = self.alts.index.get(f'{ha}>{hb}')
            if i is None:
                res = None
            else:
                arr = self.alts.alt[i].astype(float)       # [date, k, (tof, vinf_dep, vinf_arr)]
                vdep, varr = arr[:, :, 1], arr[:, :, 2]
                f = np.vectorize(lambda node, v: planet_side_dv(node, v) if not math.isnan(v) else math.nan)
                res = (arr[:, :, 0], f(a, vdep), f(b, varr))
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
    if ship.low_thrust:
        raise ValueError('low-thrust ships fly continuous-thrust legs; use engine.lowthrust for them')
    net = net or Network(alts)
    step = alts.step
    # the departure window covers the network dates inside it: the first on or after its start, the last on or before its end
    k0 = max(0, math.ceil((depart_from_jd - alts.jd0) / step - 1e-9))
    k_dep_end = min(alts.n_dates - 1, math.floor((depart_until_jd - alts.jd0) / step + 1e-9))
    # the search runs on network dates; the deadline itself is checked against each leg's actual arrival time
    arrive_by_jd = arrive_by_jd or alts.jd(k0) + 15 * 365.25
    k_end = min(alts.n_dates - 1, alts.k_of(arrive_by_jd) + 1)
    refuel = set(FUEL_SOURCES if refuel_at is None else refuel_at)
    wait = int(max_wait_days / step)
    fuel0 = ship.fuel_capacity_kg if start_fuel_kg is None else start_fuel_kg
    start = Label(origin, k0, fuel0, fuel0 if origin in refuel else 0.0, [])
    fronts = {n: [] for n in NODES}
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
        last_dep = k_dep_end if not L.legs else min(k_end, L.k + wait)
        if L.k > last_dep:
            continue
        for b in NODES:
            if b == L.node or (b != destination and b not in refuel):   # stops only where the ship can refuel
                continue
            e = net.edge(L.node, b)
            if e is None:
                continue
            tof, dv1, dv2 = (x[L.k:last_dep + 1] for x in e)
            if tof.size == 0:
                continue
            # mass accounting for both burns on every date and alternative at once
            m0 = dry + pay + fuel
            m2 = m0 * np.exp(-(dv1 + dv2) / ve)
            fuel_after = fuel - (m0 - m2)
            kd = np.arange(L.k, last_dep + 1)[:, None] + np.zeros_like(tof)
            k2 = kd + np.ceil(tof / step)
            arr_jd = alts.jd0 + kd * step + tof
            ok = ~np.isnan(tof) & (fuel_after >= 0) & (k2 <= k_end) & (arr_jd <= arrive_by_jd)
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
                leg = {'from': L.node, 'to': b, 'k_dep': kdv, 'k_arr': k2v, 'alt': alt,
                       'tof_days': float(e[0][kdv, alt]), 'dv_dep': float(e[1][kdv, alt]), 'dv_arr': float(e[2][kdv, alt]),
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

def describe(alts, ship, L):
    """itinerary with dates, burns and the mass history. Each burn is recomputed with Ship.burn as a check of the
    search's vectorised mass accounting, and fuel must carry over from leg to leg (or be refilled to capacity)."""
    legs, prev_after = [], None
    for g in L.legs:
        f0 = g['fuel_before']
        if prev_after is not None:
            expected = ship.fuel_capacity_kg if g['refuelled_before'] else prev_after
            assert abs(f0 - expected) < 1e-6 * ship.fuel_capacity_kg + 1e-6, 'fuel does not carry over between legs'
        b1 = burn(ship, f0, g['dv_dep'], parking_period_s=parking_period_s('earth' if g['from'] == 'eml1' else g['from']))
        b2 = burn(ship, b1.fuel_after_kg, g['dv_arr'], parking_period_s=parking_period_s('earth' if g['to'] == 'eml1' else g['to']))
        assert b1.feasible and b2.feasible, 'infeasible burn in a returned journey'
        assert abs(b2.fuel_after_kg - g['fuel_after']) < 1e-6 * ship.fuel_capacity_kg + 1e-6, 'mass accounting mismatch'
        legs.append({'from': g['from'], 'to': g['to'], 'from_name': NAMES[g['from']], 'to_name': NAMES[g['to']],
                     'depart_jd': alts.jd(g['k_dep']), 'arrive_jd': alts.jd(g['k_dep']) + g['tof_days'],
                     'tof_days': round(g['tof_days'], 1), 'dv_departure': round(g['dv_dep'], 3), 'dv_arrival': round(g['dv_arr'], 3),
                     'fuel_before_kg': round(f0, 1), 'fuel_after_kg': round(b2.fuel_after_kg, 1), 'refuelled_before': g['refuelled_before'],
                     'warnings': b1.warnings + b2.warnings})
        prev_after = b2.fuel_after_kg
    first = legs[0]['depart_jd']; last = legs[-1]['arrive_jd']
    return {'arrive_jd': last, 'depart_jd': first, 'days': round(last - first, 1), 'legs': legs,
            'dv_total': round(sum(l['dv_departure'] + l['dv_arrival'] for l in legs), 3),
            'fuel_loaded_kg': round(L.loaded, 1), 'fuel_left_kg': round(prev_after, 1),
            'payload_delivered_kg': ship.payload_kg, 'stops': [l['to'] for l in legs[:-1]]}

def max_payload(ship, journey, hi_kg=None, tol_kg=1.0):
    """largest payload the ship can carry through the journey's legs and refuel stops (start with full tanks).
    Same mass accounting as the search; bisection on the payload."""
    import dataclasses
    def feasible(pay):
        s = dataclasses.replace(ship, payload_kg=pay)
        fuel = s.fuel_capacity_kg
        for leg in journey['legs']:
            if leg['refuelled_before']:
                fuel = s.fuel_capacity_kg
            for dv in (leg['dv_departure'], leg['dv_arrival']):
                b = burn(s, fuel, dv)
                if not b.feasible:
                    return False
                fuel = b.fuel_after_kg
        return True
    lo, hi = 0.0, hi_kg or 50 * (ship.dry_kg + ship.fuel_capacity_kg)
    if not feasible(lo):
        return 0.0
    while hi - lo > tol_kg:
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if feasible(mid) else (lo, mid)
    return lo
