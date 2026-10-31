# The JPL Small-Body Database catalogue: every asteroid and comet with orbital elements (full precision) and the physical
# data JPL lists (absolute magnitude, diameter, albedo, rotation period, GM, density, spectral class), plus near-Earth
# and potentially-hazardous flags and the orbit condition code. Basis of the catalogue screening (tier B) and of the
# work units for contributors. Output (not in git): data/sbdb/catalog.npz (numbers) and data/sbdb/names.json.gz.
# Run on the host (network access): python3 flyby/fetch_sbdb.py
import gzip, json, pathlib, time, urllib.parse, urllib.request

OUT = pathlib.Path(__file__).resolve().parent.parent / 'data' / 'sbdb'
FIELDS = ['spkid', 'full_name', 'pdes', 'kind', 'neo', 'pha', 'H', 'diameter', 'albedo', 'rot_per', 'GM', 'density', 'spec_B', 'spec_T',
          'class', 'a', 'e', 'i', 'om', 'w', 'ma', 'q', 'tp', 'epoch', 'condition_code']
NUMERIC = ['H', 'diameter', 'albedo', 'rot_per', 'GM', 'density', 'a', 'e', 'i', 'om', 'w', 'ma', 'q', 'tp', 'epoch', 'condition_code']
PAGE = 100000

def page(start):
    q = {'fields': ','.join(FIELDS), 'full-prec': 'true', 'limit': str(PAGE), 'limit-from': str(start)}
    url = 'https://ssd-api.jpl.nasa.gov/sbdb_query.api?' + urllib.parse.urlencode(q)
    return json.load(urllib.request.urlopen(url, timeout=900))

if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    rows, start, total = [], 0, None
    while total is None or start < total:
        d = page(start); total = int(d['count']); rows += d['data']; start += PAGE
        print(f'{len(rows):,} of {total:,}', flush=True)
        time.sleep(1)
    col = {f: [r[i] for r in rows] for i, f in enumerate(FIELDS)}
    num = lambda v: float(v) if v not in (None, '') else float('nan')
    meta = {'fetched': time.strftime('%Y-%m-%d'), 'source': 'https://ssd-api.jpl.nasa.gov/sbdb_query.api (full-prec)', 'count': total,
            'units': {'a': 'au', 'q': 'au', 'i': 'deg', 'om': 'deg', 'w': 'deg', 'ma': 'deg', 'tp': 'JD TDB', 'epoch': 'JD TDB',
                      'diameter': 'km', 'rot_per': 'h', 'GM': 'km^3/s^2', 'density': 'g/cm^3', 'H': 'mag'}}
    with gzip.open(OUT / 'names.json.gz', 'wt') as f:
        json.dump({'meta': meta, 'spkid': col['spkid'], 'full_name': [s.strip() for s in col['full_name']], 'pdes': col['pdes'],
                   'kind': col['kind'], 'neo': col['neo'], 'pha': col['pha'], 'spec_B': col['spec_B'], 'spec_T': col['spec_T'],
                   'class': col['class']}, f)
    with open(OUT / 'catalog_numeric.json', 'w') as f:              # converted to npz in the container (numpy)
        json.dump({k: [num(v) for v in col[k]] for k in NUMERIC}, f)
    print('saved', total, 'objects')
