# Route dataset for the learned surrogate (phase 3).
# Each sample: two heliocentric bodies at a departure epoch -> cheapest rendezvous transfer over a
# flight-time grid (single-revolution Lambert, pykep). Features are the two states at departure;
# targets are departure vinf, arrival vinf, their sum and the best flight time.
# Orbits are sampled at random so the model can generalise to depot orbits that do not exist yet.
import sys, time
from multiprocessing import Pool
import numpy as np
import pykep as pk

K2 = 2.9591220828e-4
AU, DAY = pk.AU, pk.DAY2SEC
VSCALE = AU / DAY                        # m/s per AU/day

def random_elements(rng):
    a = float(np.exp(rng.uniform(np.log(0.3), np.log(10.0))))           # AU, log-uniform
    e = float(rng.uniform(0.0, 0.6) * min(1.0, (a - 0.25) / a))        # keep perihelion above about 0.25 AU
    i = float(np.radians(rng.uniform(0.0, 20.0)))
    W, w, M = (float(rng.uniform(0.0, 2 * np.pi)) for _ in range(3))
    return [a * AU, e, i, W, w, M]

def tof_grid(a1, a2, n=24):
    at = (a1 + a2) / 2
    tH = np.pi * np.sqrt(at ** 3 / K2)
    lo, hi = max(25.0, 0.3 * tH), min(4000.0, 1.8 * tH)
    return lo * (hi / lo) ** (np.arange(n) / (n - 1))

def sample(seed):
    rng = np.random.default_rng(seed)
    ep = pk.epoch(0.0)
    e1, e2 = random_elements(rng), random_elements(rng)
    b1 = pk.planet(pk.udpla.keplerian(ep, e1, pk.MU_SUN, 'a', [-1, -1, -1], pk.el_type.KEP_M))
    b2 = pk.planet(pk.udpla.keplerian(ep, e2, pk.MU_SUN, 'b', [-1, -1, -1], pk.el_type.KEP_M))
    r1, v1 = b1.eph(0.0); r2, v2 = b2.eph(0.0)
    best = None
    for T in tof_grid(e1[0] / AU, e2[0] / AU):
        rb, vb = b2.eph(T)
        try:
            l = pk.lambert_problem(r1, rb, T * DAY, pk.MU_SUN)
        except Exception:
            continue
        vd = np.linalg.norm(np.array(l.v0[0]) - v1) / 1000
        va = np.linalg.norm(np.array(l.v1[0]) - vb) / 1000
        if best is None or vd + va < best[0] + best[1]:
            best = (vd, va, T)
    if best is None:
        return None
    feats = list(np.array(r1) / AU) + list(np.array(v1) / VSCALE) + list(np.array(r2) / AU) + list(np.array(v2) / VSCALE)
    return feats + [best[0], best[1], best[0] + best[1], best[2]]

def main():
    n, seed0, procs = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    t = time.time()
    with Pool(procs) as p:
        rows = [r for r in p.map(sample, range(seed0, seed0 + n), chunksize=200) if r is not None]
    cols = ['r1x', 'r1y', 'r1z', 'v1x', 'v1y', 'v1z', 'r2x', 'r2y', 'r2z', 'v2x', 'v2y', 'v2z', 'vinf_dep', 'vinf_arr', 'dv_total', 'tof_days']
    np.savez_compressed(f'/work/surrogate/dataset_{seed0}_{n}.npz', X=np.array(rows, dtype=np.float32), columns=np.array(cols))
    print(f'{len(rows)} samples in {time.time() - t:.0f} s')

if __name__ == '__main__':
    main()
