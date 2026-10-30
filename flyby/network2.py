# Direct transfer network on DE440 positions, 2026-2099: the Python port of network/compute.js.
# Same nodes, parking orbits, Oberth costs and flight-time grid; output keeps the network.json schema.
import json, math, sys, time
from multiprocessing import Pool
sys.path.insert(0, '/work/flyby')
import numpy as np
import pykep as pk
import bodies

AUD = 1731.456837                       # km/s per AU/day (used only for documentation parity)
NODES = ['earth', 'venus', 'mercury', 'mars', 'ceres', 'vesta', 'psyche', '433', '101955', '162173', 'jupiter', 'saturn']
# GM in km^3/s^2 (NASA NSSDC) and parking radius in km, as in network/nodes.js
PARK = {'earth': (398600.4, 6378 + 300), 'mars': (42828.4, 3396 + 300), 'venus': (324858.6, 6052 + 300),
        'mercury': (22031.8, 2440 + 200), 'jupiter': (126686534, 1882700), 'saturn': (37931187, 1221900), 'ceres': (62.63, 470 + 50)}
PARK_LABEL = {'earth': 'low Earth orbit, 300 km', 'mars': 'low Mars orbit, 300 km', 'venus': 'low Venus orbit, 300 km',
              'mercury': 'low Mercury orbit, 200 km', 'jupiter': 'circular orbit at Callisto distance',
              'saturn': 'circular orbit at Titan distance', 'ceres': 'low Ceres orbit, 50 km'}
JD0 = 2461345.5                          # 2026-11-01
JD1 = 2488433.5                          # 2099-01-01
STEP, TOFS = 10, 30
K2 = 2.9591220828e-4

def well(vinf, node):
    if node not in PARK:
        return vinf
    mu, r = PARK[node]
    return math.sqrt(vinf * vinf + 2 * mu / r) - math.sqrt(mu / r)

def mean_a(node, atlas):
    b = {x['id']: x for x in atlas['bodies']}[node]
    return b['orbit']['semi_major_axis_au']

def tof_grid(a1, a2):
    at = (a1 + a2) / 2
    tH = math.pi * math.sqrt(at ** 3 / K2)
    lo, hi = max(25, 0.3 * tH), min(3200, 1.8 * tH)
    return [lo * (hi / lo) ** (i / (TOFS - 1)) for i in range(TOFS)]

def route(args):
    a, b = args
    atlas = bodies.load_atlas()
    A, B = bodies.body(a, atlas), bodies.body(b, atlas)
    grid = tof_grid(mean_a(a, atlas), mean_a(b, atlas))
    dep, arr, tof = [], [], []
    for jd in np.arange(JD0, JD1 + 0.5, STEP):
        t0 = jd - 2451544.5
        ra, va = A.eph(t0)
        best = None
        for T in grid:
            rb, vb = B.eph(t0 + T)
            try:
                l = pk.lambert_problem(ra, rb, T * pk.DAY2SEC, pk.MU_SUN)
            except Exception:
                continue
            vd = np.linalg.norm(np.array(l.v0[0]) - va) / 1000
            vr = np.linalg.norm(np.array(l.v1[0]) - vb) / 1000
            d, r = well(vd, a), well(vr, b)
            if best is None or d + r < best[0] + best[1]:
                best = (d, r, T)
        dep.append(round(best[0], 2) if best else None)
        arr.append(round(best[1], 2) if best else None)
        tof.append(round(best[2]) if best else None)
    return {'from': a, 'to': b, 'dep': dep, 'arr': arr, 'tof': tof}

def main():
    atlas = bodies.load_atlas()
    pairs = [(a, b) for a in NODES for b in NODES if a != b]
    t = time.time()
    with Pool(8) as p:
        routes = p.map(route, pairs)
    names = {x['id']: x['name'] for x in atlas['bodies']}
    n = len(routes[0]['dep'])
    out = {'meta': {'built': time.strftime('%Y-%m-%d'), 'jd0': JD0, 'step_days': STEP, 'n': n, 'tof_steps': TOFS,
                    'ephemeris': 'JPL DE440 (pykep de440s, planet-system barycentres) for planets; JPL SBDB elements with two-body motion for small bodies',
                    'method': 'Single-revolution Lambert transfers (pykep) between DE440 positions; best of a log-spaced flight-time grid (0.3 to 1.8 x Hohmann time, 25 to 3200 days). Departure and arrival burns from the stated parking orbits via the Oberth relation sqrt(vinf^2 + 2mu/r) - sqrt(mu/r); small bodies use rendezvous vinf. Impulsive burns, no gravity assists, no low-thrust.'},
           'nodes': [{'id': k, 'name': names[k], 'park': PARK_LABEL.get(k, 'rendezvous (gravity negligible)')} for k in NODES],
           'routes': routes}
    json.dump(out, open('/work/network/network.json', 'w'), separators=(',', ':'))
    print('routes', len(routes), 'dates', n, 'seconds', round(time.time() - t))

if __name__ == '__main__':
    main()
