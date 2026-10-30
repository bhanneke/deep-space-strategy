import json, pathlib
here = pathlib.Path(__file__).parent
atlas = json.load(open(here.parent / 'atlas' / 'atlas_dataset.json'))
net = json.load(open(here / 'network.json'))
stn = json.load(open(here / 'stations.json'))
keep = set(n['id'] for n in net['nodes'])
atlas_small = {'meta': atlas['meta'], 'bodies': [b for b in atlas['bodies'] if b['id'] in keep]}
for b in atlas_small['bodies']: b.pop('resources', None)
j = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
html = (here / 'growth_template.html').read_text()
flyby = []
for p in sorted((here.parent / 'flyby').glob('routes_*.json')):
    flyby += json.load(open(p))
html = html.replace('/*ENGINE*/', (here / 'engine.js').read_text()).replace('/*ATLAS*/', j(atlas_small)).replace('/*NET*/', j(net)).replace('/*STN*/', j(stn)).replace('/*FLYBY*/', j(flyby))
(here / 'solar-system-growth-map.html').write_text(html)
print('built', len(html), 'bytes')
