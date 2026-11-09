"""Build site/library.json: the tagged library behind Deep Space Strategy.
Every item gets a permanent four-character tag, kept in site/tags.json (append only: a tag never changes or returns).
Items come only from what solmaps publishes: places and resources from the resource atlas dataset, routes between the
30 transfer bodies, and the computation jobs of the status page (statements verbatim, progress from progress.json).
Usage: python3 site/build_library.py /path/to/solmaps"""
import json
import pathlib
import re
import sys
import urllib.request

here = pathlib.Path(__file__).parent
solmaps = pathlib.Path(sys.argv[1])
ALPHABET = '0123456789ABCDEFGHJKMNPQRSTVWXYZ'   # no I, L, O, U: easy to read and to type


def tag_of(n):
    s = ''
    for _ in range(4):
        s = ALPHABET[n % 32] + s
        n //= 32
    return s


reg_path = here / 'tags.json'
registry = json.loads(reg_path.read_text()) if reg_path.exists() else {}
next_n = max((sum(ALPHABET.index(c) * 32 ** (3 - i) for i, c in enumerate(t)) for t in registry.values()), default=0) + 1


def tag(key):
    global next_n
    if key not in registry:
        registry[key] = tag_of(next_n)
        next_n += 1
    return registry[key]


SOL = 'https://bhanneke.github.io/solmaps/'
with urllib.request.urlopen(SOL + 'progress.json', timeout=30) as r:
    progress = json.load(r)
atlas = json.loads((solmaps / 'atlas/atlas_dataset.json').read_text())
src = (solmaps / 'engine/conventions.py').read_text()
sites = re.findall(r"'([^']+)'", re.search(r'^SITES = \[(.*?)\]', src, re.M).group(1))
new = re.findall(r"'([^']+)'", re.search(r'^NEW_SMALL\s*=\s*\[(.*?)\]', src, re.M | re.S).group(1))
helio = sites + ['uranus', 'neptune', 'pluto'] + new
small = json.loads((solmaps / 'engine/small_bodies.json').read_text())['bodies']
by_id = {b['id']: b for b in atlas['bodies']}
name_of = lambda h: by_id[h]['name'] if h in by_id else small[h]['short']

CHECK = {'page-checked': 'sourced', 'derived': 'derived', 'snippet-only': 'unchecked', 'unverified': 'unchecked'}
items, used_by = {}, {}


def add(key, kind, title, statement, status, rests_on=(), link=None, source=None, extra=None):
    t = tag(key)
    items[t] = {'tag': t, 'kind': kind, 'title': title, 'statement': statement, 'status': status,
                'rests_on': [tag(k) for k in rests_on]}
    if link:
        items[t]['link'] = link
    if source:
        items[t]['source'] = source
    if extra:
        items[t]['extra'] = extra
    return t


# computations: statements verbatim from the status page
JOBS = [
    ('transfers', 'Daily transfer options between 30 bodies, 2026–2100, every number of revolutions, flights up to 100 years', 30),
    ('screening', 'Every asteroid and comet of the JPL Small-Body Database against Earth, Mars and Ceres, both directions, every day to 2100', 791),
    ('flyby-returns', 'Flyby search: returns to Earth and Mars from Ceres, Vesta, Psyche and Jupiter', 288),
    ('flyby-earth', 'Flyby search: from Earth to the planets, the belt and the near-Earth asteroids', 396),
    ('flyby-mars', 'Flyby search: from Mars to the planets and the belt', 288),
    ('flyby-saturn', 'Flyby search: from Saturn to Earth and Mars', 72),
    ('flyby-other', 'Flyby search: every other pair of the 30 bodies (841 pairs)', 30276),
    ('cyclers', 'Earth–Mars cyclers, 2027–2097', 1),
    ('lowthrust', 'Low-thrust atlas: an electric ship on every site pair and two-year window', 1),
]
for key, text, total in JOBS:
    done = len(progress.get(key, {}).get('done', []))
    add(f'job:{key}', 'Computation', text, text, 'computed' if done >= total else 'running' if done else 'open',
        link=SOL + 'status.html', extra=f'{done} of {total} units')

# places and their resources (resource atlas); moons, populations and Lagrange points included
for b in atlas['bodies']:
    o = b.get('orbit') or {}
    p = b.get('parent')
    add(f"place:{b['id']}", 'Place', b['name'],
        f"{b['name']}: {b.get('type', 'body')}" + (f" of {by_id[p]['name']}" if p in by_id else '') +
        (f", semi-major axis {o['semi_major_axis_au']:g} AU" if isinstance(o.get('semi_major_axis_au'), (int, float)) else ''),
        CHECK.get(o.get('verification'), 'unchecked'), rests_on=[f'place:{p}'] if p in by_id else [],
        link=SOL + 'atlas.html', source=o.get('source'))
for b in atlas['bodies']:
    for i, r in enumerate(b.get('resources') or []):
        add(f"res:{b['id']}:{r['resource']}:{i}", 'Resource', f"{r['resource']} on {b['name']}",
            f"{b['name']}, {r['resource']}: {r.get('form_location', '')}", CHECK.get(r.get('verification'), 'unchecked'),
            rests_on=[f"place:{b['id']}"], link=SOL + 'atlas.html', source=r.get('source'),
            extra=f"evidence: {r.get('evidence', 'not graded')}")

# routes between the 30 transfer bodies: both directions
for h in helio:
    if h not in by_id:
        add(f'place:{h}', 'Place', name_of(h), f'{name_of(h)}: {small[h].get("class") or "small body"}, semi-major axis {small[h]["a_au"]:g} AU',
            'sourced', link=SOL + f'solar-system-3d.html#body={h}', source=small[h].get('source'))
done = set(progress['transfers']['done'])
for a in helio:
    for b in helio:
        if a != b:
            add(f'route:{a}:{b}', 'Route', f'{name_of(a)} to {name_of(b)}',
                f'Transfer options from {name_of(a)} to {name_of(b)} for every day from 2026 to 2100.',
                'computed' if a in done else 'open', rests_on=[f'place:{a}', f'place:{b}', 'job:transfers'],
                link=SOL + f'solar-system-3d.html#route={a}-{b}')

for t, it in items.items():
    for r in it['rests_on']:
        used_by.setdefault(r, []).append(t)
for t, u in used_by.items():
    if t in items:
        items[t]['used_by'] = u
reg_path.write_text(json.dumps(registry, indent=0, ensure_ascii=False))
(here / 'library.json').write_text(json.dumps(list(items.values()), separators=(',', ':'), ensure_ascii=False))
kinds, stats = {}, {}
for it in items.values():
    kinds[it['kind']] = kinds.get(it['kind'], 0) + 1
    stats[it['status']] = stats.get(it['status'], 0) + 1
print('=== Summary ===')
print(len(items), 'items', kinds, stats, f"{(here / 'library.json').stat().st_size // 1024} kB")
