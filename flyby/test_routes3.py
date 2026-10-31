# Checks for routes3: the generated sequence lists, the screen cost of known mission sequences, and one quick solve
# per arrival variant (the aerocapture variant must not cost more than capture for the same route).
import sys
sys.path.insert(0, '/work/flyby')
import routes3 as R, pykep as pk

if __name__ == '__main__':
    for o, t in [('earth', 'jupiter'), ('earth', 'saturn'), ('earth', 'mercury'), ('ceres', 'earth'), ('jupiter', 'earth'), ('psyche', 'mars'), ('mars', 'saturn')]:
        print(o, '->', t, ['-'.join(x[:2] for x in s) for s in R.sequences(o, t)])
    for s in [['earth', 'venus', 'earth', 'earth', 'jupiter'], ['earth', 'venus', 'venus', 'earth', 'jupiter']]:
        print('screen cost', '-'.join(s), round(R.screen_cost(s), 2), '| direct', round(R.screen_cost([s[0], s[-1]]), 2))
    t0 = [pk.epoch('2031-01-01').mjd2000, pk.epoch('2033-01-01').mjd2000]
    for seq in [['ceres', 'mars', 'earth'], ['vesta', 'earth']]:
        a = R.solve(seq, t0, 'capture', islands=2, gens=300, pop=20, evolves=1, seed=7)
        b = R.solve(seq, t0, 'aero', islands=2, gens=300, pop=20, evolves=1, seed=7)
        print('-'.join(seq), 'capture', a['dv'], 'vinf_arr', a['vinf_arr'], '| aero', b['dv'], 'vinf_arr', b['vinf_arr'])
