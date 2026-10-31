"""External check against NASA NHATS (JPL SSD API): for each published NHATS round trip to one of our near-Earth asteroids,
solve the same two legs (Earth to the asteroid, and back after the stay) with solmaps' ephemerides and Lambert solver,
and compare the four excess speeds NHATS reports. NHATS uses its own trajectory software and JPL ephemerides, so
agreement tests positions, timing and solver against an independent source.

Step 1 (host, network): python3 tests/nhats_check.py fetch   -> tests/nhats_cases.json
Step 2 (container):     python tests/nhats_check.py           -> table of differences
"""
import json, math, pathlib, sys, urllib.parse, urllib.request
HERE = pathlib.Path(__file__).resolve().parent
NEAS = ['99942', '65803', '25143', '341843', '175706', '4660', '101955', '162173', '433']

def fetch():
    cases = []
    for des in NEAS:
        url = 'https://ssd-api.jpl.nasa.gov/nhats.api?' + urllib.parse.urlencode({'des': des, 'dv': 12, 'dur': 450, 'stay': 8, 'launch': '2020-2045'})
        try:
            d = json.load(urllib.request.urlopen(url, timeout=120))
        except Exception as e:
            print(des, 'not in NHATS for these limits:', e); continue
        if 'min_dv_traj' not in d:
            print(des, 'no trajectories:', str(d)[:120]); continue
        for kind in ('min_dv_traj', 'min_dur_traj'):
            t = d[kind]
            cases.append({'des': des, 'kind': kind, 'source': url, **{k: t[k] for k in t}})
        print(des, d.get('fullname'), 'min dv', d['min_dv']['dv'], 'km/s')
    json.dump({'fetched': __import__('time').strftime('%Y-%m-%d'), 'cases': cases}, open(HERE / 'nhats_cases.json', 'w'), indent=1)

def check():
    sys.path.insert(0, str(HERE.parent / 'flyby'))
    import pykep as pk, bodies
    atlas = bodies.load_atlas()
    cases = json.load(open(HERE / 'nhats_cases.json'))['cases']
    E = bodies.body('earth', atlas)
    worst = 0.0
    print(f"{'object':8s} {'trajectory':12s} {'launch':10s} {'out':>4s} {'stay':>4s} {'back':>4s} | v_dep_earth  vrel_arr_neo  vrel_dep_neo  entry speed (NHATS / solmaps, km/s)")
    for c in cases:
        A = bodies.body(c['des'], atlas)
        y, m, d = map(int, c['launch'].split('-'))
        jd = 367 * y - (7 * (y + (m + 9) // 12)) // 4 + (275 * m) // 9 + d + 1721013.5
        t0 = jd - 2451544.5; t1, ts, t2 = float(c['dur_out']), float(c['dur_at']), float(c['dur_ret'])
        if jd < 2461041.5:                                   # before 2026-01-01: outside the Horizons files
            print(f"{c['des']:8s} {c['kind'][4:]:12s} {c['launch']:10s} launch before 2026, skipped"); continue
        r0, v0 = E.eph(t0); r1, v1 = A.eph(t0 + t1)
        l = pk.lambert_problem(r0=r0, r1=r1, tof=t1 * 86400, mu=pk.MU_SUN)
        a1, a2 = math.dist(l.v0[0], v0) / 1000, math.dist(l.v1[0], v1) / 1000
        r2, v2 = A.eph(t0 + t1 + ts); r3, v3 = E.eph(t0 + t1 + ts + t2)
        l2 = pk.lambert_problem(r0=r2, r1=r3, tof=t2 * 86400, mu=pk.MU_SUN)
        b1, b2 = math.dist(l2.v0[0], v2) / 1000, math.dist(l2.v1[0], v3) / 1000
        # NHATS's Earth return dV is a burn only when the entry speed would exceed its limit; compare the entry speed at
        # 125 km instead: sqrt(vinf^2 + 2 GM / (R + 125 km))
        entry = math.sqrt(b2 * b2 + 2 * 398600.4 / (6378.1 + 125))
        ref = [float(c['v_dep_earth']), float(c['vrel_arr_neo']), float(c['vrel_dep_neo']), float(c['v_arr_earth'])]
        ours = [a1, a2, b1, entry]
        worst = max(worst, max(abs(x - y) for x, y in zip(ref, ours)))
        print(f"{c['des']:8s} {c['kind'][4:]:12s} {c['launch']:10s} {t1:4.0f} {ts:4.0f} {t2:4.0f} | " +
              '  '.join(f'{x:5.3f}/{y:5.3f}' for x, y in zip(ref, ours)))
    print('largest difference', round(worst, 4), 'km/s over', len(cases), 'trajectories')
    return worst

if __name__ == '__main__':
    fetch() if len(sys.argv) > 1 and sys.argv[1] == 'fetch' else check()
