"""Daily transfer alternatives for every ordered pair of the heliocentric bodies (precompute).

For each pair of conventions.HELIO (30 bodies, 870 ordered pairs) and each day from 2026-11-01 to 2100-12-31, the module
solves Lambert's problem (single revolution, prograde, pykep) over 100 flight times, log-spaced from 0.3 to 2.5 Hohmann
times and capped at 100 years, and keeps the Pareto front of flight time against cost, thinned to at most K. The cost
for the front uses each body's own end state (conventions.end_state); moons and Earth-Moon L1 share their planet's
front and are priced at query time. Arrivals after the end of a body's ephemeris (Bennu: 2135-09-30) are skipped.

Output, one file per origin (resumable): data/daily/<origin>.npz with
  targets [n]               the other bodies
  tof     [n, date, K]      flight time, days (float32; NaN where fewer than K alternatives)
  vinf    [n, date, K, 2]   excess speed at departure and arrival, km/s (float16: 0.004 km/s resolution at 5 km/s)
The cheapest priced alternative per date is the network entry; the 10-day page files come from derive_pages().

Run: python -m engine.daily [procs]   (in the solmaps-env container, from /work)
"""
import json, math, os, sys, time
from multiprocessing import Pool
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work')
import numpy as np

from engine.conventions import HELIO, SITES, NAMES, PARKING, planet_side_dv, end_state

JD0, JD1 = 2461345.5, 2488433.5                    # 2026-11-01 to 2100-12-31, every day
K, NTOF, TOF_MAX = 6, 100, 100 * 365.25
K2 = 2.9591220828e-4
OUT = '/work/data/daily'
SMALL = json.load(open('/work/engine/small_bodies.json'))['bodies']

def semi_major_axis(node):
    if node in SMALL:
        return SMALL[node]['a_au']
    import bodies
    return {b['id']: b for b in bodies.load_atlas()['bodies']}[node]['orbit']['semi_major_axis_au']

def tof_grid(a1, a2):
    """100 flight times, log-spaced from 0.3 to 2.5 Hohmann times (at least 10 days, at most 100 years)"""
    tH = math.pi * math.sqrt(((a1 + a2) / 2) ** 3 / K2)
    lo, hi = max(10.0, 0.3 * tH), min(TOF_MAX, 2.5 * tH)
    return [lo * (hi / lo) ** (i / (NTOF - 1)) for i in range(NTOF)]

def pareto_thin(T, cost, vd, va):
    """indices of the Pareto front (faster and cheaper), increasing flight time, thinned to K keeping both ends"""
    order = np.lexsort((cost, T)); out, best = [], math.inf
    for i in order:
        if cost[i] < best - 1e-9:
            out.append(i); best = cost[i]
    if len(out) > K:
        idx = sorted(set([0, len(out) - 1] + [round(j * (len(out) - 1) / (K - 1)) for j in range(K)]))
        out = [out[j] for j in idx][:K]
    return out

def route(args):
    import pykep as pk, bodies
    a, b = args
    atlas = bodies.load_atlas()
    A, B = bodies.body(a, atlas), bodies.body(b, atlas)
    end_b = bodies.coverage_end(b)
    grid = np.array(tof_grid(semi_major_axis(a), semi_major_axis(b)))
    dates = np.arange(JD0, JD1 + 0.5, 1.0)
    tof = np.full((len(dates), K), np.nan, np.float32)
    vinf = np.full((len(dates), K, 2), np.nan, np.float16)
    for j, jd in enumerate(dates):
        t0 = jd - 2451544.5
        ra, va = A.eph(t0)
        T, cost, vd, vr = [], [], [], []
        for t in grid:
            if t0 + t > end_b:
                break
            rb, vb = B.eph(t0 + t)
            try:
                l = pk.lambert_problem(r0=ra, r1=rb, tof=t * 86400.0, mu=pk.MU_SUN)
            except Exception:
                continue
            d = math.dist(l.v0[0], va) / 1000; r = math.dist(l.v1[0], vb) / 1000
            T.append(t); vd.append(d); vr.append(r); cost.append(planet_side_dv(a, d) + planet_side_dv(b, r))
        if not T:
            continue
        # refine the three cheapest local minima of the grid, each between its grid neighbours (two rounds of halving
        # in log time). Cost valleys next to the 180-degree transfer singularity can be narrower than the grid spacing:
        # Psyche to Mercury has one about 10 days wide at 378 days of flight.
        def solve(t):
            if t0 + t > end_b or t in T:
                return
            rb, vb = B.eph(t0 + t)
            try:
                l = pk.lambert_problem(r0=ra, r1=rb, tof=t * 86400.0, mu=pk.MU_SUN)
            except Exception:
                return
            d = math.dist(l.v0[0], va) / 1000; r = math.dist(l.v1[0], vb) / 1000
            T.append(t); vd.append(d); vr.append(r); cost.append(planet_side_dv(a, d) + planet_side_dv(b, r))
        n = len(T)
        minima = [q for q in range(n) if (q == 0 or cost[q] <= cost[q - 1]) and (q == n - 1 or cost[q] <= cost[q + 1])]
        for q in sorted(minima, key=lambda q: cost[q])[:3]:
            lo, mid, hi = T[max(q - 1, 0)], T[q], T[min(q + 1, n - 1)]
            for _ in range(2):
                for t in (math.sqrt(lo * mid), math.sqrt(mid * hi)):
                    solve(t)
                near = sorted(t for t in T if lo <= t <= hi)
                mid = min(near, key=lambda t: cost[T.index(t)]); m = near.index(mid)
                lo, hi = near[max(m - 1, 0)], near[min(m + 1, len(near) - 1)]
        idx = pareto_thin(np.array(T), np.array(cost), vd, vr)
        for i, p in enumerate(idx):
            tof[j, i] = T[p]; vinf[j, i, 0] = vd[p]; vinf[j, i, 1] = vr[p]
    return a, b, tof, vinf

def main(procs=8):
    os.makedirs(OUT, exist_ok=True)
    todo = [a for a in HELIO if not os.path.exists(f'{OUT}/{a}.npz')]
    pairs = [(a, b) for a in todo for b in HELIO if a != b]
    print('origins to do', len(todo), 'pairs', len(pairs), flush=True)
    t = time.time(); got = {}
    with Pool(procs) as p:
        for a, b, tof, vinf in p.imap_unordered(route, pairs):
            got.setdefault(a, {})[b] = (tof, vinf)
            if len(got[a]) == len(HELIO) - 1:            # origin complete: write it and free the memory
                targets = [x for x in HELIO if x != a]
                np.savez_compressed(f'{OUT}/{a}.tmp.npz', targets=np.array(targets), jd0=JD0, step=1.0, k=K,
                                    tof=np.stack([got[a][x][0] for x in targets]), vinf=np.stack([got[a][x][1] for x in targets]))
                os.replace(f'{OUT}/{a}.tmp.npz', f'{OUT}/{a}.npz')
                del got[a]
                print(f'{a} done, {time.time() - t:.0f} s', flush=True)
    print('all origins done', round(time.time() - t), 's')

def derive_pages(step=10):
    """the 12-site, 10-day files the pages and the 10-day engine data use, from the daily data:
    network/network.json (cheapest priced transfer per route and date) and data/alternatives.npz (the alternatives)"""
    import bodies
    atlas = bodies.load_atlas(); names = {x['id']: x['name'] for x in atlas['bodies']}
    routes, alt, pairs = [], [], []
    for a in SITES:
        d = np.load(f'{OUT}/{a}.npz'); tg = list(d['targets'])
        for b in SITES:
            if a == b:
                continue
            i = tg.index(b)
            T = d['tof'][i, ::step].astype(np.float64); V = d['vinf'][i, ::step].astype(np.float64)
            dep = np.vectorize(lambda v: planet_side_dv(a, v) if v == v else np.nan)(V[..., 0])
            arr = np.vectorize(lambda v: planet_side_dv(b, v) if v == v else np.nan)(V[..., 1])
            tot = np.where(np.isnan(T), np.inf, dep + arr); k = tot.argmin(axis=1); ok = np.isfinite(tot.min(axis=1))
            n = np.arange(len(k))
            routes.append({'from': a, 'to': b,
                           'dep': [round(float(x), 2) if o else None for x, o in zip(dep[n, k], ok)],
                           'arr': [round(float(x), 2) if o else None for x, o in zip(arr[n, k], ok)],
                           'tof': [round(float(x)) if o else None for x, o in zip(T[n, k], ok)]})
            alt.append(np.concatenate([T[..., None], V], axis=-1).astype(np.float32)); pairs.append(f'{a}>{b}')
    n = len(routes[0]['dep'])
    out = {'meta': {'built': time.strftime('%Y-%m-%d'), 'jd0': JD0, 'step_days': step, 'n': n, 'tof_steps': NTOF,
                    'ephemeris': 'JPL DE440 for planets (Horizons DE441 barycentres after 2149); daily JPL Horizons states for small bodies',
                    'method': 'Cheapest of the daily transfer alternatives (engine/daily.py): single-revolution Lambert transfers (pykep), '
                              '100 flight times from 0.3 to 2.5 Hohmann times (at most 100 years). Departure and arrival burns from the '
                              'stated parking orbits via the Oberth relation sqrt(vinf^2 + 2mu/r) - sqrt(mu/r); small bodies use '
                              'rendezvous vinf. Impulsive burns, no gravity assists, no low-thrust.'},
           'nodes': [{'id': k, 'name': names[k], 'park': end_state(k)['label']} for k in SITES], 'routes': routes}
    json.dump(out, open('/work/network/network.json', 'w'), separators=(',', ':'))
    np.savez_compressed('/work/data/alternatives.npz', alt=np.stack(alt), pairs=np.array(pairs), jd0=JD0, step=step, k=K)
    print('pages: routes', len(routes), 'dates', n)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'pages':
        derive_pages()
    else:
        main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
