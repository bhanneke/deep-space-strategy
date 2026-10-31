import json, pathlib
here = pathlib.Path(__file__).parent
stn = json.load(open(here / 'stations.json'))
stn['network'].pop('qubo', None)              # the QUBO coefficients stay in stations.json for the solver work
stn.pop('control', None)
j = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
html = (here / 'stations_template.html').read_text().replace('/*STN*/', j(stn))
(here / 'station-siting.html').write_text(html)
print('built', len(html), 'bytes')
