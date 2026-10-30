# Flyby routes from low Earth orbit to outer planets, the main belt and Mercury.
# Each route: one deep-space manoeuvre per leg (pykep mga_1dsm), optimised with pygmo.
# Costs match the direct network: departure burn from a 300 km Earth orbit via the Oberth
# relation, capture into the same parking orbits at planets, rendezvous at small bodies.
import json, math, sys, time
import pykep as pk, pygmo as pg

ATLAS = json.load(open('/work/atlas/atlas_dataset.json'))
BY = {b['id']: b for b in ATLAS['bodies']}
MU_E, R_LEO = 398600.4, 6378 + 300          # km^3/s^2, km
# capture orbits, same as network/nodes.js (radius in km)
PARK = {'mars': 3396 + 300, 'venus': 6052 + 300, 'mercury': 2440 + 200, 'jupiter': 1882700, 'saturn': 1221900}

OUT_OF_RANGE = 999.0                         # penalty, km/s: routes that need positions after 2050

def leo_burn(vinf):                          # km/s
    return math.sqrt(vinf * vinf + 2 * MU_E / R_LEO) - math.sqrt(MU_E / R_LEO)

def small_body(bid):
    el = BY[bid]['orbit']['elements']
    ep = pk.epoch(el['epoch_jd_tdb'] - 2451544.5)          # MJD2000
    elem = [el['a_au'] * pk.AU, el['e'], math.radians(el['i_deg']), math.radians(el['node_deg']),
            math.radians(el['arg_peri_deg']), math.radians(el['M_deg'])]
    return pk.planet(pk.udpla.keplerian(ep, elem, pk.MU_SUN, bid, [-1, -1, -1], pk.el_type.KEP_M))

def body(n):
    return pk.planet(pk.udpla.jpl_lp(n)) if n in ('mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn') else small_body(n)

A = {'mercury': 0.387, 'venus': 0.723, 'earth': 1.0, 'mars': 1.524, 'ceres': 2.77, 'vesta': 2.36, 'psyche': 2.92, 'jupiter': 5.20, 'saturn': 9.54}
def hohmann_days(a1, a2):
    at = (a1 + a2) / 2
    return math.pi * math.sqrt(at ** 3 / 2.9591220828e-4)
def leg_bounds(p, q):
    if p == q:                                              # resonant return to the same planet
        per = 365.25 * A[p] ** 1.5
        return [0.8 * per, 2.2 * per]
    t = hohmann_days(A[p], A[q])
    return [max(30, 0.3 * t), min(4000, 1.9 * t)]

SEQS = {
    'jupiter': [['earth', 'jupiter'], ['earth', 'venus', 'jupiter'], ['earth', 'earth', 'jupiter'], ['earth', 'mars', 'jupiter'],
                ['earth', 'venus', 'earth', 'jupiter'], ['earth', 'venus', 'venus', 'jupiter']],
    'saturn':  [['earth', 'saturn'], ['earth', 'jupiter', 'saturn'], ['earth', 'venus', 'saturn'], ['earth', 'earth', 'saturn'],
                ['earth', 'venus', 'jupiter', 'saturn'], ['earth', 'earth', 'jupiter', 'saturn']],
    'ceres':   [['earth', 'ceres'], ['earth', 'mars', 'ceres'], ['earth', 'earth', 'ceres'], ['earth', 'venus', 'ceres'], ['earth', 'earth', 'mars', 'ceres']],
    'vesta':   [['earth', 'vesta'], ['earth', 'mars', 'vesta'], ['earth', 'earth', 'vesta'], ['earth', 'earth', 'mars', 'vesta']],
    'psyche':  [['earth', 'psyche'], ['earth', 'mars', 'psyche'], ['earth', 'earth', 'psyche'], ['earth', 'earth', 'mars', 'psyche']],
    'mercury': [['earth', 'mercury'], ['earth', 'venus', 'mercury'], ['earth', 'venus', 'venus', 'mercury'], ['earth', 'earth', 'venus', 'mercury']],
}

class Route:
    """pygmo problem: LEO departure burn + deep-space manoeuvres + arrival burn, km/s."""
    def __init__(self, seq, t0):
        target = seq[-1]
        kw = dict(seq=[body(n) for n in seq], t0=t0, tof=[leg_bounds(a, b) for a, b in zip(seq, seq[1:])],
                  vinf=[0.5, 9.0], add_vinf_dep=False, tof_encoding='direct')
        if target in PARK:
            # pykep requires add_vinf_arr with orbit insertion; the arrival term is then the insertion burn
            kw.update(add_vinf_arr=True, orbit_insertion=True, e_target=0.0, rp_target=PARK[target] * 1000.0)
        else:
            kw.update(add_vinf_arr=True)
        self.inner = pk.trajopt.mga_1dsm(**kw)
    def fitness(self, x):
        try:
            return [self.inner.fitness(x)[0] / 1000.0 + leo_burn(x[3] / 1000.0)]
        except ValueError:                                  # an encounter after 2050, outside the ephemeris range
            return [OUT_OF_RANGE]
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
    t_launch = x[0]
    return {'seq': seq, 'dv': round(f, 3), 'launch_mjd2000': round(t_launch, 1),
            'vinf_dep': round(x[3] / 1000, 3), 'leo_burn': round(leo_burn(x[3] / 1000), 3), 'tof_days': round(sum(tofs), 1),
            'arrival_mjd2000': round(t_launch + sum(tofs), 1)}

if __name__ == '__main__':
    targets = sys.argv[1].split(',')
    budget = dict(islands=int(sys.argv[2]), gens=int(sys.argv[3]), pop=int(sys.argv[4]), evolves=int(sys.argv[5]))
    first, last, width = 2027, 2047, 2                   # two-year launch bins, launches up to end of 2048
    out = []
    for tgt in targets:
        for y in range(first, last + 1, width):
            t0 = [pk.epoch(f'{y}-01-01').mjd2000, pk.epoch(f'{y + width}-01-01').mjd2000]
            for seq in SEQS[tgt]:
                t = time.time()
                # two independent starts, keep the better one
                r = min((solve(seq, t0, seed=y * 10 + k, **budget) for k in range(2)), key=lambda s: s['dv'])
                r['feasible'] = bool(r['dv'] < OUT_OF_RANGE)
                r.update(target=tgt, bin=f'{y}-{y + width - 1}', seconds=round(time.time() - t, 1))
                out.append(r)
                print(json.dumps(r), flush=True)
    json.dump(out, open(f'/work/flyby/routes_{"_".join(targets)}.json', 'w'), indent=1)
