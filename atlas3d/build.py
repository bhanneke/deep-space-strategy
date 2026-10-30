import json, pathlib
here = pathlib.Path(__file__).parent
root = here.parent
j = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
data = json.load(open(root / 'atlas' / 'atlas_dataset.json'))
net = json.load(open(root / 'network' / 'network.json'))
flyby = []
for p in sorted((root / 'flyby').glob('routes2_*.json')):
    flyby += json.load(open(p))
html = (here / 'solar3d_template.html').read_text()
html = html.replace('/*ATLAS_DATA*/', j(data)).replace('/*NET*/', j(net)).replace('/*FLYBY*/', j(flyby)).replace('/*ENGINE*/', (root / 'network' / 'engine.js').read_text())
(here / 'solar-system-3d.html').write_text(html)
print('built', len(html), 'bytes;', len(flyby), 'flyby routes')
