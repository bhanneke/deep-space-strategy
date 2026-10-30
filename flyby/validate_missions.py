# Replay real flyby missions: give the optimiser only the planet sequence and a launch year,
# then compare its encounter dates with the flown ones. Planet positions: JPL DE440 (bodies.py).
import sys
sys.path.insert(0, '/work/flyby')
import pykep as pk, pygmo as pg
import bodies

ATLAS = bodies.load_atlas()
MISSIONS = {
    'Juno': (['earth', 'earth', 'jupiter'], ('2011-01-01', '2012-01-01'), [[600, 900], [700, 1300]],
             ['2011-08-05', '2013-10-09', '2016-07-04']),
    'Galileo': (['earth', 'venus', 'earth', 'earth', 'jupiter'], ('1989-06-01', '1990-03-01'), [[60, 200], [200, 400], [600, 800], [900, 1500]],
                ['1989-10-18', '1990-02-10', '1990-12-08', '1992-12-08', '1995-12-07']),
    'Dawn': (['earth', 'mars', 'vesta'], ('2007-06-01', '2008-03-01'), [[300, 700], [600, 1200]],
             ['2007-09-27', '2009-02-17', '2011-07-16']),
}

def run(name):
    seq, (a, b), tof, real = MISSIONS[name]
    udp = pk.trajopt.mga_1dsm(seq=[bodies.body(n, ATLAS) for n in seq], t0=[pk.epoch(a).mjd2000, pk.epoch(b).mjd2000],
                              tof=tof, vinf=[0.5, 7.0], add_vinf_dep=True, add_vinf_arr=True, tof_encoding='direct')
    archi = pg.archipelago(n=8, algo=pg.algorithm(pg.sade(gen=3000, variant=6, variant_adptv=2, memory=True)),
                           prob=pg.problem(udp), pop_size=30, seed=3)
    for _ in range(4):
        archi.evolve(); archi.wait_check()
    pop = min((isl.get_population() for isl in archi), key=lambda p: p.champion_f[0])
    x, n = pop.champion_x, len(seq) - 1
    tofs = [x[5]] + [x[9 + 4 * (i - 1)] for i in range(1, n)]
    t = [x[0]]
    for d in tofs:
        t.append(t[-1] + d)
    print(f'{name}: objective {pop.champion_f[0] / 1000:.2f} km/s (launch vinf + manoeuvres + arrival vinf)')
    for body, r, m in zip(seq, real, t):
        print(f'   {body:8s} real {r}   model {str(pk.epoch(m))[:10]}   difference {m - pk.epoch(r).mjd2000:+6.0f} days')

def main():
    for name in sys.argv[1].split(','):
        run(name)

if __name__ == '__main__':      # pygmo starts worker processes that re-import this module
    main()
