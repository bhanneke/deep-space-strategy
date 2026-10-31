# Optimiser study, part 2: two-stage hybrids. Stage 1 is a population search (fast, finds the region); stage 2 refines
# the five best routes of every island with monotonic basin hopping and compass search. Same problems and repeats as
# optimizer_study.py; best known values come from both studies.
import sys, json, time
sys.path.insert(0, '/work/flyby')
import pygmo as pg, pykep as pk
import routes3 as R
from optimizer_study import PROBLEMS

def stage1(udp, ring, pop, gens, seed):
    kw = dict(t=pg.ring()) if ring else {}
    a = pg.archipelago(n=8, algo=pg.algorithm(pg.sade(gen=gens, variant=6, variant_adptv=2, memory=True)), prob=pg.problem(udp), pop_size=pop, seed=seed, **kw)
    for _ in range(6):
        a.evolve(); a.wait_check()
    return a

def polish(udp, a, seed, rounds=2):
    prob = pg.problem(udp)
    b = pg.archipelago()
    for i, isl in enumerate(a):
        src = isl.get_population()
        best = sorted(range(len(src)), key=lambda j: src.get_f()[j][0])[:5]
        pop = pg.population(prob, 0, seed=seed + i)
        for j in best:
            pop.push_back(src.get_x()[j], src.get_f()[j])
        inner = pg.algorithm(pg.compass_search(max_fevals=10000, start_range=0.05, stop_range=1e-6))
        b.push_back(algo=pg.algorithm(pg.mbh(inner, stop=5, perturb=0.05, seed=seed + i)), pop=pop)
    for _ in range(rounds):
        b.evolve(); b.wait_check()
    return min(isl.get_population().champion_f[0] for isl in b)

if __name__ == '__main__':
    for seq, y in PROBLEMS:
        t0 = [pk.epoch(f'{y}-01-01').mjd2000, pk.epoch(f'{y + 2}-01-01').mjd2000]
        udp = R.Route(seq, t0, 'capture')
        for name, ring, pop, gens in [('E ring double population + polish', True, 60, 1000), ('F unconnected + polish', False, 30, 2000)]:
            t = time.time(); vals = []
            for k in range(6):
                a = stage1(udp, ring, pop, gens, seed=y * 1000 + k)
                vals.append(round(polish(udp, a, seed=y * 1000 + 100 + k), 3))
            print(json.dumps({'seq': '-'.join(seq), 'year': y, 'setting': name, 'values': vals, 'seconds_per_run': round((time.time() - t) / 6)}), flush=True)
