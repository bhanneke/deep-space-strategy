# How much the Horizons orbits change the direct network: routes that touch a small body, old (two-body) against new.
import json, statistics as st, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
A = {(r['from'], r['to']): r for r in a['routes']}; B = {(r['from'], r['to']): r for r in b['routes']}
small = {'ceres', 'vesta', 'psyche', '433', '101955', '162173'}
n = b['meta']['n']
def cost(r, k): return None if r['dep'][k] is None else r['dep'][k] + r['arr'][k]
for period, (k0, k1) in {'2026-2040': (0, 511), '2041-2070': (511, 1607), '2071-2100': (1607, n)}.items():
    d = sorted(abs(cost(A[key], k) - cost(B[key], k)) for key in B if key[0] in small or key[1] in small
               for k in range(k0, k1, 5) if cost(A[key], k) is not None and cost(B[key], k) is not None)
    print(period, 'routes touching small bodies: median change %.2f km/s, 90th percentile %.2f, max %.2f' % (st.median(d), d[int(.9 * len(d))], d[-1]))
for key in [('earth', 'ceres'), ('earth', 'psyche'), ('ceres', 'vesta'), ('earth', '101955')]:
    for k0, k1, lab in [(0, 365, '2026-36'), (1607, 1972, '2071-81')]:
        o = min(cost(A[key], k) for k in range(k0, k1) if cost(A[key], k) is not None)
        m = min(cost(B[key], k) for k in range(k0, k1) if cost(B[key], k) is not None)
        print(key, lab, 'best window old %.2f new %.2f' % (o, m))
