# Perturbed heliocentric state vectors from JPL Horizons (ecliptic J2000, AU and AU/day), every day.
#   - small-body sites, 2026 to 2205: bodies.py interpolates them; this replaces two-body motion from SBDB elements,
#     which ignores the pull of Jupiter and the other planets and drifts over decades. Daily states keep close planetary
#     encounters sharp (Apophis passes Earth at 38,000 km on 2029-04-13), which a 10-day step would smooth over.
#   - planet-system barycentres, 2140 to 2205: DE440 in pykep ends in 2150; flights of up to 100 years from departures
#     before 2100 arrive until about 2200, so bodies.py switches to these Horizons states (DE441) after 2149.
# Also writes engine/small_bodies.json: name, orbit size and spectral type from the JPL Small-Body Database (SBDB).
# The files are not in git; each data release archives the exact copies it used. Run on the host (network access):
#   python3 flyby/fetch_horizons.py
import json, pathlib, re, time, urllib.parse, urllib.request

SITES = {
    # existing sites
    'ceres': '1', 'vesta': '4', 'psyche': '16', '433': '433', '101955': '101955', '162173': '162173',
    # main belt, water-rich or large
    '2': '2', '10': '10', '24': '24', '3': '3',
    # near-Earth asteroids, easy to reach or hydrated
    '99942': '99942', '65803': '65803', '25143': '25143', '341843': '341843', '175706': '175706', '4660': '4660',
    # Jupiter Trojans visited by Lucy
    '3548': '3548', '15094': '15094', '11351': '11351', '21900': '21900', '617': '617',
}
# planet-system barycentres as in pykep's de440s (Earth: the Earth-Moon barycentre)
PLANETS = {'mercury': '1', 'venus': '2', 'earth': '3', 'mars': '4', 'jupiter': '5', 'saturn': '6', 'uranus': '7',
           'neptune': '8', 'pluto': '9'}
ROOT = pathlib.Path(__file__).parent
OUT = ROOT / 'horizons'

def fetch(cmd, start, stop, small=True):
    q = {'format': 'json', 'COMMAND': f"'{cmd};'" if small else f"'{cmd}'", 'OBJ_DATA': "'NO'", 'MAKE_EPHEM': "'YES'",
         'EPHEM_TYPE': "'VECTORS'", 'CENTER': "'500@10'", 'REF_PLANE': "'ECLIPTIC'", 'REF_SYSTEM': "'ICRF'",
         'OUT_UNITS': "'AU-D'", 'VEC_TABLE': "'2'", 'CSV_FORMAT': "'YES'", 'START_TIME': f"'{start}'",
         'STOP_TIME': f"'{stop}'", 'STEP_SIZE': "'1 d'"}
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api?' + urllib.parse.urlencode(q)
    r = json.load(urllib.request.urlopen(url, timeout=600))['result']
    m = re.search(r'No ephemeris for target .* after A\.D\. (\d{4})-([A-Z]{3})-(\d{2})', r)
    if m and '$$SOE' not in r:
        # Horizons ends some solutions early (Bennu: 2135-09-30, after its close Earth approach); stop the day before
        mon = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'].index(m.group(2)) + 1
        end = time.strftime('%Y-%m-%d', time.gmtime(time.mktime((int(m.group(1)), mon, int(m.group(3)), 12, 0, 0, 0, 0, 0)) - 86400))
        print('  Horizons ends this solution on', m.group(0)[-30:].strip(), '-> fetching to', end, flush=True)
        return fetch(cmd, start, end, small)
    source = r[r.index('Target body name'):].split('\n')[0]
    rows = []
    for line in r[r.index('$$SOE') + 5:r.index('$$EOE')].strip().split('\n'):
        f = [x.strip() for x in line.split(',')]
        rows.append([float(f[0])] + [float(x) for x in f[2:8]])
    return source, rows

def save(name, source, rows):
    json.dump({'site': name, 'source': source.strip(), 'frame': 'heliocentric ecliptic J2000, AU and AU/day, TDB',
               'fetched': time.strftime('%Y-%m-%d'), 'columns': ['jd_tdb', 'x', 'y', 'z', 'vx', 'vy', 'vz'], 'rows': rows},
              open(OUT / f'{name}.json', 'w'))

def sbdb(num):
    url = 'https://ssd-api.jpl.nasa.gov/sbdb.api?' + urllib.parse.urlencode({'sstr': num, 'phys-par': 'true'})
    d = json.load(urllib.request.urlopen(url, timeout=120))
    el = {e['name']: e['value'] for e in d['orbit']['elements']}
    phys = {p['name']: p['value'] for p in d.get('phys_par', [])}
    return {'name': d['object']['fullname'], 'short': d['object'].get('shortname') or d['object']['fullname'],
            'a_au': float(el['a']), 'e': float(el['e']), 'i_deg': float(el['i']),
            'class': d['object'].get('orbit_class', {}).get('name'), 'spectral': phys.get('spec_T') or phys.get('spec_B'),
            'source': url}

if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    reg = {}
    for site, cmd in SITES.items():
        source, rows = fetch(cmd, '2026-01-01', '2205-12-31')
        save(site, source, rows)
        reg[site] = sbdb(cmd)
        print(site, len(rows), 'rows', rows[0][0], '-', rows[-1][0], '|', reg[site]['short'], flush=True)
        time.sleep(1)                                    # be gentle with the JPL APIs
    json.dump({'fetched': time.strftime('%Y-%m-%d'), 'bodies': reg}, open(ROOT.parent / 'engine' / 'small_bodies.json', 'w'), indent=1)
    for planet, cmd in PLANETS.items():
        source, rows = fetch(cmd, '2140-01-01', '2205-12-31', small=False)
        save('planet_' + planet, source, rows)
        print('planet', planet, len(rows), 'rows |', source.strip()[:60], flush=True)
        time.sleep(1)
