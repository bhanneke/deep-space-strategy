# Perturbed heliocentric state vectors for the small-body sites from JPL Horizons (ecliptic J2000, AU and AU/day),
# every 10 days from 2026 to 2110. bodies.py interpolates them; this replaces two-body motion from SBDB elements,
# which ignores the pull of Jupiter and the other planets and drifts over decades.
import json, pathlib, urllib.parse, urllib.request

SITES = {'ceres': '1;', 'vesta': '4;', 'psyche': '16;', '433': '433;', '101955': '101955;', '162173': '162173;'}
OUT = pathlib.Path(__file__).parent / 'horizons'

def fetch(cmd):
    q = {'format': 'json', 'COMMAND': f"'{cmd}'", 'OBJ_DATA': "'NO'", 'MAKE_EPHEM': "'YES'", 'EPHEM_TYPE': "'VECTORS'",
         'CENTER': "'500@10'", 'REF_PLANE': "'ECLIPTIC'", 'REF_SYSTEM': "'ICRF'", 'OUT_UNITS': "'AU-D'", 'VEC_TABLE': "'2'",
         'CSV_FORMAT': "'YES'", 'START_TIME': "'2026-01-01'", 'STOP_TIME': "'2110-12-31'", 'STEP_SIZE': "'10 d'"}
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urllib.parse.urlencode(q)
    r = json.load(urllib.request.urlopen(url, timeout=120))['result']
    source = r[r.index('Target body name'):].split('\n')[0]
    rows = []
    for line in r[r.index('$$SOE') + 5:r.index('$$EOE')].strip().split('\n'):
        f = [x.strip() for x in line.split(',')]
        rows.append([float(f[0])] + [float(x) for x in f[2:8]])
    return source, rows

if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    for site, cmd in SITES.items():
        source, rows = fetch(cmd)
        json.dump({'site': site, 'source': source.strip(), 'frame': 'heliocentric ecliptic J2000, AU and AU/day, TDB',
                   'columns': ['jd_tdb', 'x', 'y', 'z', 'vx', 'vy', 'vz'], 'rows': rows}, open(OUT / f'{site}.json', 'w'))
        print(site, len(rows), 'rows', rows[0][0], '-', rows[-1][0], '|', source.strip()[:60])
