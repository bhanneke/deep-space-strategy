# Validate the pykep/pygmo flyby toolchain on pykep's built-in benchmark problems.
import sys, time
import pykep as pk, pygmo as pg

def solve(udp, gens, islands, pop, evolves):
    prob = pg.problem(udp)
    archi = pg.archipelago(n=islands, algo=pg.algorithm(pg.sade(gen=gens, variant=6, variant_adptv=2, memory=True)), prob=prob, pop_size=pop, seed=7)
    for _ in range(evolves):
        archi.evolve(); archi.wait_check()
    return min(isl.get_population().champion_f[0] for isl in archi)

def main():
    gens, islands, pop, evolves = (int(a) for a in sys.argv[1:5])
    for name in sys.argv[5].split(','):
        t0 = time.time()
        f = solve(getattr(pk.trajopt.gym, name), gens, islands, pop, evolves)
        print(f'{name:10s} best found {f:10.4f}   ({time.time() - t0:5.0f} s)', flush=True)

if __name__ == '__main__':      # pygmo starts worker processes that re-import this module
    main()
