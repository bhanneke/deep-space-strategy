# Queue duration on the Mac Mini from measured seconds per problem (all starts included), by sequence length.
import sys
sys.path.insert(0, '/work/flyby')
import routes3 as R
SEC = {2: 43, 3: 117, 4: 307}                      # measured on the Mini, 2026-10-31 (ceres:earth,mars job)
JOBS = [('ceres', 'earth,mars'), ('vesta', 'earth,mars'), ('psyche', 'earth,mars'), ('jupiter', 'earth,mars'),
        ('earth', 'jupiter,saturn,mercury'), ('earth', 'ceres,vesta,psyche'), ('mars', 'jupiter,saturn'), ('mars', 'ceres,vesta,psyche')]
total = 0
for o, ts in JOBS:
    h = 0
    for t in ts.split(','):
        for s in R.sequences(o, t):
            h += 36 * (2 if t in R.AERO else 1) * SEC[len(s)] / 3600
    total += h
    print(f'{o}:{ts}  {h:.0f} h')
print(f'total {total:.0f} h = {total / 24:.1f} days')
