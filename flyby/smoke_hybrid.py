# One hybrid solve (routes3 setting) on a hard problem: checks the code path and lets the caller measure memory.
import sys, time
sys.path.insert(0, '/work/flyby')
import routes3 as R, pykep as pk

if __name__ == '__main__':
    R.init_pool()
    t = time.time()
    r = R.solve(['earth', 'earth', 'jupiter'], [pk.epoch('2067-01-01').mjd2000, pk.epoch('2069-01-01').mjd2000], 'capture',
                islands=8, gens=1000, pop=60, evolves=6, seed=2067000)
    print('dv', r['dv'], 'vinf_arr', r['vinf_arr'], 'seconds', round(time.time() - t))
    r = R.solve(['vesta', 'earth'], [pk.epoch('2051-01-01').mjd2000, pk.epoch('2053-01-01').mjd2000], 'aero',
                islands=8, gens=1000, pop=60, evolves=6, seed=7)
    print('vesta-earth aero dv', r['dv'], 'vinf_arr', r['vinf_arr'])
