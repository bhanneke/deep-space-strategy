"""Build site/orbits.json: real orbital elements for the landing page figure.
Planets: JPL Table 1 elements with rates, from the solmaps atlas dataset (valid 1800-2050).
The other solmaps transfer bodies: JPL Small-Body Database (sbdb.api).
Background populations: numbered objects from the JPL SBDB query API (sbdb_query.api).
Usage: python3 site/fetch_orbits.py /path/to/solmaps   (network access; run when the data should refresh)"""
import json
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request

solmaps = pathlib.Path(sys.argv[1])
OUT = pathlib.Path(__file__).with_name('orbits.json')
UA = {'User-Agent': 'deep-space-strategy-site/0.1'}
EPOCH = 2461000.5     # common epoch (JD) for the stored mean anomalies
TWO_BODY = 0.9856076686  # Gaussian mean motion, degrees per day for a = 1 AU


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
        return json.loads(r.read())


def to_epoch(ma, ep, a):
    return round((float(ma) + TWO_BODY / float(a) ** 1.5 * (EPOCH - float(ep))) % 360, 2)


src = (solmaps / 'engine/conventions.py').read_text()
sites = re.findall(r"'([^']+)'", re.search(r'^SITES = \[(.*?)\]', src, re.M).group(1))
new = re.findall(r"'([^']+)'", re.search(r'^NEW_SMALL\s*=\s*\[(.*?)\]', src, re.M | re.S).group(1))
helio = sites + ['uranus', 'neptune', 'pluto'] + new

atlas = {b['id']: b for b in json.loads((solmaps / 'atlas/atlas_dataset.json').read_text())['bodies']}
small = json.loads((solmaps / 'engine/small_bodies.json').read_text())['bodies']

planets, others = [], []
for h in helio:
    orbit = (atlas.get(h) or {}).get('orbit') or {}
    if orbit.get('j2000_elements') and orbit.get('rates_per_julian_century') and h != 'pluto':
        e, r = orbit['j2000_elements'], orbit['rates_per_julian_century']
        planets.append({'id': h, 'name': atlas[h]['name'],
                        'el': [e['a_au'], e['e'], e['I_deg'], e['L_deg'], e['long_peri_deg'], e['long_node_deg']],
                        'rate': [r['a_au_per_cy'], r['e_per_cy'], r['I_deg_per_cy'], r['L_deg_per_cy'], r['long_peri_deg_per_cy'], r['long_node_deg_per_cy']]})
        continue
    # small bodies and Pluto: SBDB elements. Ids are designation numbers; Pluto by name
    key = (small.get(h) or {}).get('short', h).split()[0] if h in small else h
    data = get('https://ssd-api.jpl.nasa.gov/sbdb.api?' + urllib.parse.urlencode({'sstr': '134340' if h == 'pluto' else key}))
    el = {x['name']: x['value'] for x in data['orbit']['elements']}
    name = atlas[h]['name'] if h in atlas else small[h]['short']
    others.append({'id': h, 'name': name, 'a': round(float(el['a']), 4), 'e': round(float(el['e']), 4), 'i': round(float(el['i']), 2),
                   'om': round(float(el['om']), 2), 'w': round(float(el['w']), 2),
                   'm': to_epoch(el['ma'], data['orbit']['epoch'], el['a'])})
    time.sleep(0.3)

populations = {}
for label, cls, limit in [('main belt', 'MBA', 2600), ('Jupiter Trojans', 'TJN', 500), ('near-Earth', 'APO', 250), ('Kuiper belt', 'TNO', 400)]:
    q = urllib.parse.urlencode({'fields': 'a,e,i,om,w,ma,epoch', 'sb-kind': 'a', 'sb-ns': 'n', 'sb-class': cls, 'limit': limit})
    rows = get('https://ssd-api.jpl.nasa.gov/sbdb_query.api?' + q)['data']
    populations[label] = [[round(float(a), 3), round(float(e), 3), round(float(i), 1), round(float(om), 1), round(float(w), 1), to_epoch(ma, ep, a)]
                          for a, e, i, om, w, ma, ep in rows if a and e and i and om and w and ma and ep]
    time.sleep(0.3)

OUT.write_text(json.dumps({'epoch': EPOCH, 'planets': planets, 'others': others, 'populations': populations,
                           'sources': ['JPL Approximate Positions of the Planets, Table 1 (valid 1800-2050)', 'JPL Small-Body Database'],
                           'fetched': time.strftime('%Y-%m-%d')}, separators=(',', ':'), ensure_ascii=False))
print('=== Summary ===')
print(f'planets {len(planets)}, other bodies {len(others)}, populations', {k: len(v) for k, v in populations.items()}, f'{OUT.stat().st_size // 1024} kB')
print('names', [b['name'] for b in others][:6], 'pluto', next(b for b in others if b['id'] == 'pluto'))
