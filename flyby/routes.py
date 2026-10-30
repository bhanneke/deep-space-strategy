# Flyby routes from low Earth orbit or low Mars orbit to outer planets, the main belt and Mercury.
# Each route: one deep-space manoeuvre per leg (pykep mga_1dsm), optimised with pygmo.
# Costs match the direct network: departure burn from the origin's parking orbit via the Oberth
# relation, capture into the same parking orbits at planets, rendezvous at small bodies.
# Planet positions: JPL DE440 (bodies.py), valid to about 2150.
import json, math, sys, time
sys.path.insert(0, '/work/flyby')
import pykep as pk, pygmo as pg
import bodies

ATLAS = bodies.load_atlas()
OUT_OF_RANGE = 999.0                         # penalty, km/s: routes outside the ephemeris range
# GM (km^3/s^2) and parking radius (km), same as network/nodes.js
PARK = {'earth': (398600.4, 6378 + 300), 'mars': (42828.4, 3396 + 300), 'venus': (324858.6, 6052 + 300),
        'mercury': (22031.8, 2440 + 200), 'jupiter': (126686534, 1882700), 'saturn': (37931187, 1221900)}

def park_burn(vinf, node):                   # km/s, burn between a parking orbit and a hyperbola
    mu, r = PARK[node]
    return math.sqrt(vinf * vinf + 2 * mu / r) - math.sqrt(mu / r)

A = {'mercury': 0.387, 'venus': 0.723, 'earth': 1.0, 'mars': 1.524, 'ceres': 2.77, 'vesta': 2.36, 'psyche': 2.92, 'jupiter': 5.20, 'saturn': 9.54}
def hohmann_days(a1, a2):
    at = (a1 + a2) / 2
    return math.pi * math.sqrt(at ** 3 / 2.9591220828e-4)
def leg_bounds(p, q):
    if p == q:                               # resonant return to the same planet
        per = 365.25 * A[p] ** 1.5
        return [0.8 * per, 2.2 * per]
    t = hohmann_days(A[p], A[q])
    return [max(30, 0.3 * t), min(4000, 1.9 * t)]

SEQS = {
    'earth': {
        'jupiter': [['earth', 'jupiter'], ['earth', 'venus', 'jupiter'], ['earth', 'earth', 'jupiter'], ['earth', 'mars', 'jupiter'],
                    ['earth', 'venus', 'earth', 'jupiter'], ['earth', 'venus', 'venus', 'jupiter']],
        'saturn':  [['earth', 'saturn'], ['earth', 'jupiter', 'saturn'], ['earth', 'venus', 'saturn'], ['earth', 'earth', 'saturn'],
                    ['earth', 'venus', 'jupiter', 'saturn'], ['earth', 'earth', 'jupiter', 'saturn']],
        'ceres':   [['earth', 'ceres'], ['earth', 'mars', 'ceres'], ['earth', 'earth', 'ceres'], ['earth', 'venus', 'ceres'], ['earth', 'earth', 'mars', 'ceres']],
        'vesta':   [['earth', 'vesta'], ['earth', 'mars', 'vesta'], ['earth', 'earth', 'vesta'], ['earth', 'earth', 'mars', 'vesta']],
        'psyche':  [['earth', 'psyche'], ['earth', 'mars', 'psyche'], ['earth', 'earth', 'psyche'], ['earth', 'earth', 'mars', 'psyche']],
        'mercury': [['earth', 'mercury'], ['earth', 'venus', 'mercury'], ['earth', 'venus', 'venus', 'mercury'], ['earth', 'earth', 'venus', 'mercury']],
    },
    'mars': {
        'jupiter': [['mars', 'jupiter'], ['mars', 'earth', 'jupiter'], ['mars', 'mars', 'jupiter'], ['mars', 'earth', 'earth', 'jupiter']],
        'saturn':  [['mars', 'saturn'], ['mars', 'jupiter', 'saturn'], ['mars', 'earth', 'jupiter', 'saturn'], ['mars', 'mars', 'jupiter', 'saturn']],
        'ceres':   [['mars', 'ceres'], ['mars', 'mars', 'ceres'], ['mars', 'earth', 'ceres']],
        'vesta':   [['mars', 'vesta'], ['mars', 'mars', 'vesta'], ['mars', 'earth', 'vesta']],
        'psyche':  [['mars', 'psyche'], ['mars', 'mars', 'psyche'], ['mars', 'earth', 'psyche']],
    },
}

class Route:
    """pygmo problem: departure burn from the origin's parking orbit + deep-space manoeuvres + arrival burn, km/s."""
    def __init__(self, seq, t0):
        self.origin, target = seq[0], seq[-1]
        kw = dict(seq=[bodies.body(n, ATLAS) for n in seq], t0=t0, tof=[leg_bounds(a, b) for a, b in zip(seq, seq[1:])],
                  vinf=[0.5, 9.0], add_vinf_dep=False, add_vinf_arr=True, tof_encoding='direct')
        if target in PARK:
            # with orbit insertion the arrival term is the insertion burn into the parking orbit
            kw.update(orbit_insertion=True, e_target=0.0, rp_target=PARK[target][1] * 1000.0)
        self.inner = pk.trajopt.mga_1dsm(**kw)
    def fitness(self, x):
        try:
            return [self.inner.fitness(x)[0] / 1000.0 + park_burn(x[3] / 1000.0, self.origin)]
        except Exception as e:
            # only an epoch outside the ephemeris range counts as infeasible; anything else must stop the run
            if 'INSUFFDATA' in str(e) or 'only valid' in str(e):
                return [OUT_OF_RANGE]
            raise
    def get_bounds(self):
        return self.inner.get_bounds()

def solve(seq, t0, islands, gens, pop, evolves, seed):
    udp = Route(seq, t0)
    archi = pg.archipelago(n=islands, algo=pg.algorithm(pg.sade(gen=gens, variant=6, variant_adptv=2, memory=True)),
                           prob=pg.problem(udp), pop_size=pop, seed=seed)
    for _ in range(evolves):
        archi.evolve(); archi.wait_check()
    best = min((isl.get_population() for isl in archi), key=lambda p: p.champion_f[0])
    x, f = best.champion_x, best.champion_f[0]
    n = len(seq) - 1
    tofs = [x[5]] + [x[9 + 4 * (i - 1)] for i in range(1, n)]
    # encounter epochs (MJD2000) let the pages draw the flyby path leg by leg; x is the full MGA-1DSM decision vector
    enc = [float(x[0])]
    for d in tofs:
        enc.append(enc[-1] + float(d))
    return {'seq': seq, 'dv': round(float(f), 3), 'launch_mjd2000': round(float(x[0]), 1),
            'encounters_mjd2000': [round(e, 2) for e in enc], 'x': [float(v) for v in x],
            'vinf_dep': round(float(x[3]) / 1000, 3), 'dep_burn': round(park_burn(x[3] / 1000, seq[0]), 3),
            'tof_days': round(float(sum(tofs)), 1), 'arrival_mjd2000': round(float(x[0] + sum(tofs)), 1)}

def main():
    origin, targets = sys.argv[1], sys.argv[2].split(',')
    first, last = int(sys.argv[3]), int(sys.argv[4])
    budget = dict(islands=int(sys.argv[5]), gens=int(sys.argv[6]), pop=int(sys.argv[7]), evolves=int(sys.argv[8]))
    width, out = 2, []
    tag = f'{origin}_{"_".join(targets)}_{first}_{last}'
    for tgt in targets:
        for y in range(first, last + 1, width):
            t0 = [pk.epoch(f'{y}-01-01').mjd2000, pk.epoch(f'{y + width}-01-01').mjd2000]
            for seq in SEQS[origin][tgt]:
                t = time.time()
                r = min((solve(seq, t0, seed=y * 10 + k, **budget) for k in range(2)), key=lambda s: s['dv'])
                r.update(origin=origin, target=tgt, bin=f'{y}-{y + width - 1}', feasible=bool(r['dv'] < OUT_OF_RANGE), seconds=round(time.time() - t, 1))
                out.append(r)
                print(json.dumps(r), flush=True)
            json.dump(out, open(f'/work/flyby/routes2_{tag}.json', 'w'), indent=1)   # save progress after each window

if __name__ == '__main__':      # pygmo starts worker processes that re-import this module
    main()
