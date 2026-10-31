"""Transfer alternatives per departure date (precompute).

For every ordered pair of the 12 sites and every network date (every 10 days, 2026-11-01 to 2100), the module solves
Lambert's problem over a grid of flight times and keeps the alternatives that no other flight time beats on both
arrival time and total cost (a Pareto front), thinned to at most K. A deadline or an onward connection can then pick a
faster, more expensive transfer; the network keeps only the cheapest.

Stored per alternative: flight time and the heliocentric excess speeds at both ends. Planet-side burns are not stored;
the engine adds them at runtime from engine/conventions.py for the end state the journey uses (parking orbit,
rendezvous or Earth-Moon L1). The flight-time grid contains the network's 30 values, so the cheapest alternative
reproduces network.json; it extends to 2.5 Hohmann times for slow, cheap transfers.

Run: python -m engine.alternatives  (in the solmaps-env container, from /work) -> data/alternatives.npz
"""
import math, sys, time
from multiprocessing import Pool
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work')
import numpy as np

from engine.conventions import SITES, planet_side_dv

JD0, JD1, STEP = 2461345.5, 2488433.5, 10          # same date grid as network2.py
K = 6                                               # alternatives kept per route and date
K2 = 2.9591220828e-4

def tof_grid(a1, a2):
    """the network's 30 flight times (0.3-1.8 Hohmann) plus 10 longer ones up to 2.5 Hohmann"""
    tH = math.pi * math.sqrt(((a1 + a2) / 2) ** 3 / K2)
    lo, hi = max(25, 0.3 * tH), min(3200, 1.8 * tH)
    base = [lo * (hi / lo) ** (i / 29) for i in range(30)]
    hi2 = min(4500, 2.5 * tH)
    extra = [hi * (hi2 / hi) ** (i / 10) for i in range(1, 11)] if hi2 > hi else []
    return base + extra

def pareto(points):
    """points: (tof, cost, vdep, varr); keep those with no other point both faster and cheaper"""
    pts = sorted(points, key=lambda p: (p[0], p[1]))
    out, best = [], math.inf
    for p in pts:
        if p[1] < best - 1e-9:
            out.append(p); best = p[1]
    return out                                           # increasing flight time, decreasing cost

def thin(front, k=K):
    """keep the cheapest and the fastest point and spread the rest evenly in flight time"""
    if len(front) <= k:
        return front
    idx = sorted(set([0, len(front) - 1] + [round(i * (len(front) - 1) / (k - 1)) for i in range(k)]))
    return [front[i] for i in idx][:k]

def route(args):
    import pykep as pk
    import bodies, network2
    a, b = args
    atlas = bodies.load_atlas()
    A, B = bodies.body(a, atlas), bodies.body(b, atlas)
    grid = tof_grid(network2.mean_a(a, atlas), network2.mean_a(b, atlas))
    dates = np.arange(JD0, JD1 + 0.5, STEP)
    out = np.full((len(dates), K, 3), np.nan, np.float32)
    for j, jd in enumerate(dates):
        t0 = jd - 2451544.5
        ra, va = A.eph(t0)
        pts = []
        for T in grid:
            rb, vb = B.eph(t0 + T)
            try:
                l = pk.lambert_problem(ra, rb, T * pk.DAY2SEC, pk.MU_SUN)
            except Exception:
                continue
            vd = float(np.linalg.norm(np.array(l.v0[0]) - va) / 1000)
            vr = float(np.linalg.norm(np.array(l.v1[0]) - vb) / 1000)
            pts.append((T, planet_side_dv(a, vd) + planet_side_dv(b, vr), vd, vr))
        for i, p in enumerate(thin(pareto(pts))):
            out[j, i] = (p[0], p[2], p[3])
    return a, b, out

def main(path='/work/data/alternatives.npz', procs=8):
    import os
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pairs = [(a, b) for a in SITES for b in SITES if a != b]
    t = time.time()
    with Pool(procs) as p:
        res = p.map(route, pairs)
    arr = np.stack([r[2] for r in res])                  # [route, date, alternative, (tof, vinf_dep, vinf_arr)]
    np.savez_compressed(path, alt=arr, pairs=np.array([f'{a}>{b}' for a, b, _ in res]), jd0=JD0, step=STEP, k=K)
    print('routes', len(res), 'dates', arr.shape[1], 'alternatives per date', K, 'seconds', round(time.time() - t))

if __name__ == '__main__':
    main(procs=int(sys.argv[1]) if len(sys.argv) > 1 else 8)
