# Compare the best flyby route with the direct route per target and launch window.
import json, sys
from collections import defaultdict

AB = {'earth': 'E', 'venus': 'V', 'mars': 'Ma', 'jupiter': 'J', 'saturn': 'S', 'mercury': 'Me', 'ceres': 'C', 'vesta': 'Ve', 'psyche': 'P'}

def main(paths):
    rows = [r for p in paths for r in json.load(open(p))]
    tab = defaultdict(dict)
    for r in rows:
        tab[(r['target'], r['bin'])]['-'.join(AB[s] for s in r['seq'])] = r
    for (t, b), d in tab.items():
        direct = next(v for k, v in d.items() if k.count('-') == 1)
        name, best = min(d.items(), key=lambda kv: kv[1]['dv'])
        dtxt = f"{direct['dv']:.2f} ({direct['tof_days']:.0f} d)" if direct['feasible'] else 'beyond 2050'
        btxt = f"{best['dv']:.2f} ({best['tof_days']:.0f} d)" if best['feasible'] else 'beyond 2050'
        save = f"{direct['dv'] - best['dv']:+.2f}" if direct['feasible'] and best['feasible'] else ''
        print(f'{t:8s} {b}  direct {dtxt:16s} best {name:10s} {btxt:16s} saving {save}')
    print('optimiser seconds', round(sum(r['seconds'] for r in rows)))

if __name__ == '__main__':
    main(sys.argv[1:])
