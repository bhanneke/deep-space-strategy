# Which optimiser setting finds the best flyby route reliably? Four hard problems from convergence.log, four settings,
# six independent repeats each. Hit rate = share of repeats within 0.05 km/s of the best value any run found.
import sys, json, time
sys.path.insert(0, '/work/flyby')
import pygmo as pg, pykep as pk
import routes3 as R

PROBLEMS = [(['earth', 'venus', 'venus', 'mercury'], 2035), (['earth', 'jupiter', 'saturn'], 2045),
            (['earth', 'earth', 'jupiter'], 2067), (['mars', 'jupiter', 'saturn'], 2071)]

def run(udp, setting, seed):
    if setting == 'A unconnected sade':
        a = pg.archipelago(n=8, algo=pg.algorithm(pg.sade(gen=2000, variant=6, variant_adptv=2, memory=True)), prob=pg.problem(udp), pop_size=30, seed=seed)
        rounds = 6
    elif setting == 'B ring sade':
        a = pg.archipelago(n=8, t=pg.ring(), algo=pg.algorithm(pg.sade(gen=2000, variant=6, variant_adptv=2, memory=True)), prob=pg.problem(udp), pop_size=30, seed=seed)
        rounds = 6
    elif setting == 'C basin hopping':
        inner = pg.algorithm(pg.compass_search(max_fevals=20000, start_range=0.1, stop_range=1e-6))
        a = pg.archipelago(n=8, algo=pg.algorithm(pg.mbh(inner, stop=10, perturb=0.1, seed=seed)), prob=pg.problem(udp), pop_size=20, seed=seed)
        rounds = 6
    else:                                   # 'D ring sade, double population'
        a = pg.archipelago(n=8, t=pg.ring(), algo=pg.algorithm(pg.sade(gen=1000, variant=6, variant_adptv=2, memory=True)), prob=pg.problem(udp), pop_size=60, seed=seed)
        rounds = 6
    for _ in range(rounds):
        a.evolve(); a.wait_check()
    return min(isl.get_population().champion_f[0] for isl in a)

if __name__ == '__main__':
    settings = ['A unconnected sade', 'B ring sade', 'C basin hopping', 'D ring sade, double population']
    for seq, y in PROBLEMS:
        t0 = [pk.epoch(f'{y}-01-01').mjd2000, pk.epoch(f'{y + 2}-01-01').mjd2000]
        udp = R.Route(seq, t0, 'capture')
        for st in settings:
            t = time.time()
            vals = [round(run(udp, st, seed=y * 1000 + k), 3) for k in range(6)]
            print(json.dumps({'seq': '-'.join(seq), 'year': y, 'setting': st, 'values': vals, 'seconds_per_run': round((time.time() - t) / 6)}), flush=True)
