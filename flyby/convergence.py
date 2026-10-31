# Convergence check for routes3: ten varied problems at the full budget, six independent starts each.
# Reports how far the best of the first three starts lies above the best of all six.
import sys, json, time
sys.path.insert(0, '/work/flyby')
import routes3 as R, pykep as pk

PROBLEMS = [(['earth', 'mars', 'vesta'], 2031, 'capture'), (['earth', 'venus', 'venus', 'mercury'], 2035, 'capture'),
            (['ceres', 'mars', 'earth'], 2033, 'aero'), (['jupiter', 'earth'], 2041, 'capture'),
            (['mars', 'mars', 'psyche'], 2029, 'capture'), (['earth', 'jupiter', 'saturn'], 2045, 'capture'),
            (['vesta', 'earth'], 2051, 'aero'), (['psyche', 'earth', 'mars'], 2061, 'capture'),
            (['earth', 'earth', 'jupiter'], 2067, 'capture'), (['mars', 'jupiter', 'saturn'], 2071, 'capture')]

if __name__ == '__main__':
    gaps = []
    for seq, y, arr in PROBLEMS:
        t0 = [pk.epoch(f'{y}-01-01').mjd2000, pk.epoch(f'{y + 2}-01-01').mjd2000]
        t = time.time()
        dv = [R.solve(seq, t0, arr, islands=8, gens=2000, pop=30, evolves=6, seed=y * 100 + k)['dv'] for k in range(6)]
        gap = min(dv[:3]) - min(dv)
        gaps.append(gap)
        print(json.dumps({'seq': '-'.join(seq), 'year': y, 'arrival': arr, 'starts': dv, 'best3': min(dv[:3]), 'best6': min(dv),
                          'gap': round(gap, 3), 'spread3': round(max(dv[:3]) - min(dv[:3]), 3), 'seconds': round(time.time() - t)}), flush=True)
    print('max gap of 3 starts against 6:', round(max(gaps), 3), 'km/s; problems with gap > 0.1:', sum(g > 0.1 for g in gaps))
