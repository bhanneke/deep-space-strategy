import json, pathlib
here = pathlib.Path(__file__).parent
root = here.parent
j = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
data = json.load(open(root / 'atlas' / 'atlas_dataset.json'))
net = json.load(open(root / 'network' / 'network.json'))
flyby = []
for p in sorted((root / 'flyby').glob('routes2_*.json')):
    flyby += json.load(open(p))
paths = {}
for p in sorted((root / 'flyby').glob('paths_*.json')):
    paths.update(json.load(open(p)))
for r in flyby:                                   # the optimiser's sampled trajectory, where trajectories.py made one
    key = f"{r.get('origin', 'earth')}|{r['target']}|{r['launch_mjd2000']}|{'-'.join(r['seq'])}"
    if key in paths:
        r['path'] = paths[key]
    r.pop('x', None)
html = (here / 'solar3d_template.html').read_text()
html = html.replace('/*ATLAS_DATA*/', j(data)).replace('/*NET*/', j(net)).replace('/*FLYBY*/', j(flyby)).replace('/*ENGINE*/', (root / 'network' / 'engine.js').read_text())
(here / 'solar-system-3d.html').write_text(html)
print('built', len(html), 'bytes;', len(flyby), 'flyby routes;', sum('path' in r for r in flyby), 'with sampled trajectories')
