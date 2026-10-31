# Smoke test for the return-route origins: three short searches, quick settings.
import sys
sys.path.insert(0, '/work/flyby'); sys.argv = ['x']
import routes, pykep as pk

if __name__ == '__main__':
    for seq in [['ceres', 'mars', 'earth'], ['jupiter', 'earth'], ['psyche', 'mars', 'mars']]:
        r = routes.solve(seq, [pk.epoch('2031-01-01').mjd2000, pk.epoch('2033-01-01').mjd2000], islands=2, gens=200, pop=20, evolves=1, seed=1)
        print(seq, r['dv'], 'dep', r['dep_burn'], 'tof', r['tof_days'])
