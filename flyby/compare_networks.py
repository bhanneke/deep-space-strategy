# Compare two network.json files on their overlapping departure dates.
import json, sys, statistics

def totals(r):
    return [None if d is None else d + a for d, a in zip(r['dep'], r['arr'])]

def main(old_path, new_path):
    old, new = json.load(open(old_path)), json.load(open(new_path))
    assert abs(old['meta']['jd0'] - new['meta']['jd0']) < 1e-6 and old['meta']['step_days'] == new['meta']['step_days']
    n = min(old['meta']['n'], new['meta']['n'])
    O = {(r['from'], r['to']): totals(r) for r in old['routes']}
    diffs, worst = [], []
    for r in new['routes']:
        key = (r['from'], r['to'])
        a, b = O[key][:n], totals(r)[:n]
        d = [abs(x - y) for x, y in zip(a, b) if x is not None and y is not None]
        diffs += d
        worst.append((max(d), statistics.median(d), key))
    worst.sort(reverse=True)
    print(f'overlap {n} dates x {len(new["routes"])} routes; median abs difference {statistics.median(diffs):.3f} km/s; '
          f'95th percentile {sorted(diffs)[int(0.95 * len(diffs))]:.3f}; share under 0.1 km/s {sum(x < 0.1 for x in diffs) / len(diffs):.1%}')
    for mx, med, key in worst[:6]:
        print(f'  {key[0]:>8s} -> {key[1]:<8s} max {mx:6.2f}  median {med:.3f}')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
