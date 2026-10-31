# Station siting, step 3: transfers between stations, and between stations and the 12 sites, with flight times,
# so journeys can stop at stations. Pool: the 50 best single relays of stations.json, every station of the best sets,
# and Earth-Moon L1 and low Venus, Mars and Ceres orbit. Each node has an end cost: a free station pays its excess
# speed (rendezvous), a low planetary orbit the Oberth burn of the planet's gravity well, Earth-Moon L1 the growth map's
# L1 model. Same dates as precompute.py (every 30 days from 2027-01-01 for 20 years).
# Output: stations/pairs.json in km/s and days.
import json, math, sys, time
from multiprocessing import Pool
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work/stations')
import numpy as np
import pykep as pk
import bodies
from network2 import NODES, well
from precompute import EPOCH_JD, DATES, station_body, from_eml1, tofs

FIXED_SITE = {'eml1': 'earth', 'venus': 'venus', 'mars': 'mars', 'ceres': 'ceres'}

def node_body(c, atlas):
    if c.get('kind') == 'site':
        return bodies.body(c['id'], atlas)
    return bodies.body(FIXED_SITE[c['site']], atlas) if c.get('site') else station_body(c)

def end_cost(c, vinf):
    if c.get('kind') == 'site':
        return well(vinf, c['id'])                 # parking orbit at planets, rendezvous at small bodies
    if not c.get('site'):
        return vinf
    return from_eml1(vinf) if c['site'] == 'eml1' else well(vinf, c['site'])

SAME = {'venus': 'venus', 'mars': 'mars', 'ceres': 'ceres'}     # a low planetary orbit and its site are the same place
L1_LEO, L1_DAYS = 3.77, 4                                          # low Earth orbit <-> Earth-Moon L1, as in the growth map

def same_place(ci, cj):
    a, b = (ci, cj) if ci.get('kind') == 'site' else (cj, ci)
    if a.get('kind') != 'site' or not b.get('site'):
        return None
    if b['site'] == 'eml1' and a['id'] == 'earth':
        return L1_LEO, L1_DAYS
    if SAME.get(b['site']) == a['id']:
        return 0.0, 0
    return None

def edge(args):
    key, ci, cj = args
    fixed = same_place(ci, cj)
    if fixed:
        return key, np.full(len(DATES), fixed[0], np.float32), np.full(len(DATES), fixed[1], np.float32)
    atlas = bodies.load_atlas()
    A, B = node_body(ci, atlas), node_body(cj, atlas)
    grid = tofs(ci['a'], cj['a'])
    cost = np.full(len(DATES), 99.0, np.float32); tof = np.zeros(len(DATES), np.float32)
    for k, t0 in enumerate(EPOCH_JD - 2451544.5 + DATES):
        ra, va = A.eph(t0)
        for T in grid:
            rb, vb = B.eph(t0 + T)
            try:
                l = pk.lambert_problem(ra, rb, T * pk.DAY2SEC, pk.MU_SUN)
            except Exception:
                continue
            c = end_cost(ci, np.linalg.norm(np.array(l.v0[0]) - va) / 1000) + end_cost(cj, np.linalg.norm(np.array(l.v1[0]) - vb) / 1000)
            if c < cost[k]:
                cost[k], tof[k] = c, T
    return key, cost, tof

def main():
    s = json.load(open('/work/stations/stations.json'))
    cands = s['candidates']
    sets = [c for st in s['network']['sets'] for c in st['stations']]
    pool = list(dict.fromkeys(s['network']['qubo']['candidates'] + sets + [k for k, c in enumerate(cands) if c.get('site')]))
    nodes = [cands[k] for k in pool]
    atlas = bodies.load_atlas()
    sa = {x['id']: x['orbit']['semi_major_axis_au'] for x in atlas['bodies'] if x['id'] in NODES}
    sites = [{'kind': 'site', 'id': n, 'a': sa[n]} for n in NODES]
    jobs = [(('nn', i, j), nodes[i], nodes[j]) for i in range(len(nodes)) for j in range(len(nodes)) if i != j]
    jobs += [(('ns', i, j), nodes[i], sites[j]) for i in range(len(nodes)) for j in range(len(sites))]
    jobs += [(('sn', j, i), sites[j], nodes[i]) for i in range(len(nodes)) for j in range(len(sites))]
    n, m, t = len(nodes), len(sites), len(DATES)
    NN = np.full((n, n, t), 99.0, np.float32); NNt = np.zeros_like(NN)
    NS = np.full((n, m, t), 99.0, np.float32); NSt = np.zeros_like(NS)
    SN = np.full((m, n, t), 99.0, np.float32); SNt = np.zeros_like(SN)
    t0 = time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as p:
        for (kind, a, b), c, f in p.imap_unordered(edge, jobs, chunksize=8):
            tgt = {'nn': (NN, NNt), 'ns': (NS, NSt), 'sn': (SN, SNt)}[kind]
            tgt[0][a, b], tgt[1][a, b] = c, f
    r = lambda x: np.round(x, 2).tolist()
    json.dump({'epoch_jd': EPOCH_JD, 'step_days': float(DATES[1] - DATES[0]), 'dates': len(DATES), 'sites': NODES,
               'nodes': [{**c, 'index': k} for k, c in zip(pool, nodes)],
               'node_node': {'cost': r(NN), 'tof': np.rint(NNt).astype(int).tolist()},
               'node_site': {'cost': r(NS), 'tof': np.rint(NSt).astype(int).tolist()},
               'site_node': {'cost': r(SN), 'tof': np.rint(SNt).astype(int).tolist()}},
              open('/work/stations/pairs.json', 'w'), separators=(',', ':'))
    print('nodes', n, 'edges', len(jobs), 'seconds', round(time.time() - t0))

if __name__ == '__main__':
    main()
