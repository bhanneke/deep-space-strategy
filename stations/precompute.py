# Station siting, step 1: transfer costs between candidate station orbits and the 12 network sites.
# For every candidate solar orbit and every site, the cheapest Lambert transfer in each direction, every 30 days
# from 2027-01-01 for 20 years, on the same DE440 positions, parking orbits and Oberth costs as network2.py.
# A station has negligible gravity, so the station end of a transfer costs its excess speed (rendezvous).
# Output: stations/costs.npz (arrays [candidate, site, date] in km/s) and stations/candidates.json.
import json, math, sys, time
from multiprocessing import Pool
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work')
import numpy as np
import pykep as pk
import bodies
from network2 import NODES, PARK, well

EPOCH_JD = 2461406.5                       # 2027-01-01
STEP, YEARS, NTOF = 30, 20, 12
DATES = np.arange(0, YEARS * 365.25, STEP)
# Earth-Moon L1 departure, as in the growth map: the cheaper of a burn at L1 and a drop to a 300 km perigee
# followed by a perigee burn (Oberth). Patched-conic estimate; lunar flybys would do slightly better.
from engine.conventions import GM_E, planet_side_dv
def from_eml1(vinf):
    return planet_side_dv('eml1', vinf)

def candidates():
    out = []
    for a in np.round(np.arange(0.6, 3.31, 0.1), 2):
        for e in (0.0, 0.1, 0.2, 0.3, 0.45):
            if a * (1 - e) < 0.45:                 # keep perihelion outside 0.45 AU
                continue
            for i in (0.0, 4.0, 8.0):
                for w in ((0.0,) if e == 0 else (0.0, 90.0, 180.0, 270.0)):
                    for M0 in (0.0, 120.0, 240.0):
                        out.append({'a': float(a), 'e': e, 'i': i, 'node': 0.0, 'w': w, 'M0': M0})
    return out

def station_body(c):
    ep = pk.epoch(EPOCH_JD - 2451544.5)
    elem = [c['a'] * pk.AU, c['e'], math.radians(c['i']), math.radians(c['node']), math.radians(c['w']), math.radians(c['M0'])]
    return pk.planet(pk.udpla.keplerian(ep, elem, pk.MU_SUN, 'stn', [-1, -1, -1], pk.el_type.KEP_M))

def site_a(s, atlas):
    return {x['id']: x for x in atlas['bodies']}[s]['orbit']['semi_major_axis_au']

def tofs(a1, a2):
    tH = math.pi * math.sqrt(((a1 + a2) / 2) ** 3 / 2.9591220828e-4)
    lo, hi = max(25, 0.3 * tH), min(3200, 1.8 * tH)
    return [lo * (hi / lo) ** (k / (NTOF - 1)) for k in range(NTOF)]

def best(A, B, t0, grid, cost):
    """cheapest transfer over the flight-time grid: (cost, flight time in days)"""
    ra, va = A.eph(t0)
    b, bt = np.inf, np.nan
    for T in grid:
        rb, vb = B.eph(t0 + T)
        try:
            l = pk.lambert_problem(ra, rb, T * pk.DAY2SEC, pk.MU_SUN)
        except Exception:
            continue
        vd = np.linalg.norm(np.array(l.v0[0]) - va) / 1000
        vr = np.linalg.norm(np.array(l.v1[0]) - vb) / 1000
        c = cost(vd, vr)
        if c < b:
            b, bt = c, T
    return b, bt

def one(args):
    k, c = args
    atlas = bodies.load_atlas()
    C = station_body(c)
    out = np.full((len(NODES), len(DATES)), np.nan, np.float32)
    inn, out_t, inn_t = np.full_like(out, np.nan), np.full_like(out, np.nan), np.full_like(out, np.nan)
    in_eml1, in_eml1_t = np.full(len(DATES), np.nan, np.float32), np.full(len(DATES), np.nan, np.float32)
    t0s = EPOCH_JD - 2451544.5 + DATES
    for s, name in enumerate(NODES):
        S = bodies.body(name, atlas)
        grid = tofs(c['a'], site_a(name, atlas))
        for j, t0 in enumerate(t0s):
            out[s, j], out_t[s, j] = best(C, S, t0, grid, lambda vd, vr: vd + well(vr, name))
            inn[s, j], inn_t[s, j] = best(S, C, t0, grid, lambda vd, vr: well(vd, name) + vr)
            if name == 'earth':
                in_eml1[j], in_eml1_t[j] = best(S, C, t0, grid, lambda vd, vr: from_eml1(vd) + vr)
    return k, out, inn, in_eml1, out_t, inn_t, in_eml1_t

def main():
    cands = candidates()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else len(cands)
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    cands = cands[:n]
    t = time.time()
    OUT = np.zeros((n, len(NODES), len(DATES)), np.float32); IN = np.zeros_like(OUT); EML = np.zeros((n, len(DATES)), np.float32)
    OUTt, INt, EMLt = np.zeros_like(OUT), np.zeros_like(OUT), np.zeros_like(EML)
    with Pool(procs) as p:
        for done, (k, o, i, e, ot, it, et) in enumerate(p.imap_unordered(one, list(enumerate(cands)), chunksize=4), 1):
            OUT[k], IN[k], EML[k], OUTt[k], INt[k], EMLt[k] = o, i, e, ot, it, et
            if done % 200 == 0:
                print(done, '/', n, round(time.time() - t), 's', flush=True)
    np.savez_compressed('/work/stations/costs.npz', out=OUT, inn=IN, in_eml1=EML, out_tof=OUTt, inn_tof=INt, in_eml1_tof=EMLt,
                        dates=DATES, epoch_jd=EPOCH_JD)
    json.dump({'meta': {'epoch_jd': EPOCH_JD, 'step_days': STEP, 'years': YEARS, 'sites': NODES, 'tof_steps': NTOF,
                        'eml1_departure': 'cheaper of a burn at Earth-Moon L1 and a perigee drop to 300 km with a perigee burn',
                        'method': 'Cheapest single-revolution Lambert transfer (pykep, DE440 sites) per 30-day departure date; '
                                  'parking-orbit Oberth costs at planets as in the network; rendezvous at the station.'},
               'candidates': cands}, open('/work/stations/candidates.json', 'w'))
    print('candidates', n, 'seconds', round(time.time() - t))

if __name__ == '__main__':
    main()
