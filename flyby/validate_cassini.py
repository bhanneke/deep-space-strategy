# Compare the optimised E-VVEJ-S trajectory (pykep gym cassini1) with the real Cassini timeline.
import pykep as pk, pygmo as pg

REAL = [('launch', '1997-10-15'), ('Venus flyby 1', '1998-04-26'), ('Venus flyby 2', '1999-06-24'),
        ('Earth flyby', '1999-08-18'), ('Jupiter flyby', '2000-12-30'), ('Saturn arrival', '2004-07-01')]

def main():
    udp = pk.trajopt.gym.cassini1
    archi = pg.archipelago(n=8, algo=pg.algorithm(pg.sade(gen=1500, variant=6, variant_adptv=2, memory=True)),
                           prob=pg.problem(udp), pop_size=20, seed=11)
    for _ in range(4):
        archi.evolve(); archi.wait_check()
    pop = min((isl.get_population() for isl in archi), key=lambda p: p.champion_f[0])
    x = pop.champion_x
    t = [x[0]]
    for tof in x[1:]:
        t.append(t[-1] + tof)
    print(f'objective {pop.champion_f[0]:.1f} m/s (first 3 km/s of launch vinf free)')
    for (label, real), mjd in zip(REAL, t):
        model = str(pk.epoch(mjd))[:10]
        diff = mjd - pk.epoch(real).mjd2000
        print(f'{label:15s} real {real}   model {model}   difference {diff:+7.0f} days')

if __name__ == '__main__':
    main()
