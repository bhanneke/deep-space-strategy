import json, pathlib
here = pathlib.Path(__file__).parent
data = json.load(open(here.parent / 'atlas' / 'atlas_dataset.json'))
blob = json.dumps(data, ensure_ascii=False).replace('</', '<\\/')
html = (here / 'solar3d_template.html').read_text().replace('/*ATLAS_DATA*/', blob)
(here / 'solar-system-3d.html').write_text(html)
print('built', len(html), 'bytes')
