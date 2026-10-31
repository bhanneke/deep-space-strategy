# Packs the fetched SBDB catalogue (flyby/fetch_sbdb.py) into data/sbdb/catalog.npz and prints what it holds.
# Run in the solmaps-env container: python flyby/sbdb_pack.py
import gzip, json, pathlib
import numpy as np

D = pathlib.Path('/work/data/sbdb')
num = json.load(open(D / 'catalog_numeric.json'))
names = json.load(gzip.open(D / 'names.json.gz', 'rt'))
arr = {k: np.array(v, dtype=np.float64) for k, v in num.items()}
flags = {k: np.array([x == 'Y' for x in names[k]]) for k in ('neo', 'pha')}
kind = np.array([str(x) for x in names['kind']])
np.savez_compressed(D / 'catalog.npz', spkid=np.array(names['spkid'], dtype=np.int64), kind=kind, **flags, **arr)
n = len(kind)
print(f"{n:,} objects; asteroids {np.sum(np.char.startswith(kind, 'a')):,}, comets {np.sum(np.char.startswith(kind, 'c')):,}")
print(f"near-Earth {flags['neo'].sum():,}, potentially hazardous {flags['pha'].sum():,}")
for k in ('diameter', 'albedo', 'rot_per', 'GM', 'density'):
    print(f'with {k}: {np.isfinite(arr[k]).sum():,}')
spec = sum(1 for b, t in zip(names['spec_B'], names['spec_T']) if b or t)
print(f'with a spectral class: {spec:,}')
bound = np.isfinite(arr['a']) & (arr['e'] < 1)
print(f"bound orbits {bound.sum():,}; orbit condition code 0-5 (well determined) {np.sum(arr['condition_code'] <= 5):,}")
