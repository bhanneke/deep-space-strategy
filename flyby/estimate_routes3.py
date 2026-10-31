# Compute estimate for the routes3 queue from the sequence lists and the timings in convergence.log.
import sys, json
sys.path.insert(0, '/work/flyby')
import routes3 as R

conv = [json.loads(l) for l in open('/work/flyby/convergence.log') if l.startswith('{')]
per_start = {}
for c in conv:
    n = c['seq'].count('-') + 1
    per_start.setdefault(n, []).append(c['seconds'] / 6)
avg = {n: sum(v) / len(v) for n, v in per_start.items()}
print('seconds per start by sequence length:', {n: round(v) for n, v in sorted(avg.items())})
JOBS = {'returns': [('ceres', 'earth'), ('ceres', 'mars'), ('vesta', 'earth'), ('vesta', 'mars'), ('psyche', 'earth'), ('psyche', 'mars'),
                    ('jupiter', 'earth'), ('jupiter', 'mars')],
        'earth': [('earth', t) for t in ('jupiter', 'saturn', 'mercury', 'ceres', 'vesta', 'psyche')],
        'mars': [('mars', t) for t in ('jupiter', 'saturn', 'ceres', 'vesta', 'psyche')]}
for name, jobs in JOBS.items():
    total = 0
    for o, t in jobs:
        for s in R.sequences(o, t):
            starts = 3 * (2 if len(s) > 3 else 1)
            variants = 2 if t in R.AERO else 1
            total += 36 * variants * starts * avg.get(len(s), max(avg.values()))
    print(f'{name}: {total / 3600:.0f} hours of one run (8 islands)')
