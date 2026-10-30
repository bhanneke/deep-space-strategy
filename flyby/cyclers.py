# Earth-Mars cycler metric: for a start window, find the encounter schedule E-M-E-M-... that keeps a
# spacecraft cycling with the least maintenance, and report what a taxi must match at each encounter.
# Maintenance = sum of powered-flyby burns (pykep fb_dv) at the intermediate encounters; launch and final
# arrival are free, because a cycler is built once and runs on. Planet positions: JPL DE440 (bodies.py).
import json, sys, time
sys.path.insert(0, '/work/flyby')
import numpy as np
import pykep as pk, pygmo as pg
import bodies

ATLAS = bodies.load_atlas()

class Cycler:
    def __init__(self, seq, t0, legs):
        self.seq = [bodies.body(n, ATLAS) for n in seq]
        self.names, self.t0, self.legs = seq, t0, legs
    def get_bounds(self):
        return [self.t0[0]] + [lo for lo, hi in self.legs], [self.t0[1]] + [hi for lo, hi in self.legs]
    def _arcs(self, x):
        t = [x[0]]
        for d in x[1:]:
            t.append(t[-1] + d)
        st = [p.eph(ti) for p, ti in zip(self.seq, t)]
        arcs = []
        for i in range(len(t) - 1):
            l = pk.lambert_problem(st[i][0], st[i + 1][0], (t[i + 1] - t[i]) * pk.DAY2SEC, pk.MU_SUN)
            arcs.append((np.array(l.v0[0]), np.array(l.v1[0])))
        return t, st, arcs
    def fitness(self, x):
        try:
            t, st, arcs = self._arcs(x)
        except Exception:
            return [99.0]
        dv = 0.0
        for i in range(1, len(t) - 1):                         # intermediate encounters: powered flybys
            vin = arcs[i - 1][1] - np.array(st[i][1])
            vout = arcs[i][0] - np.array(st[i][1])
            dv += pk.fb_dv(vin, vout, self.seq[i])
        return [dv / 1000.0]
    def report(self, x):
        t, st, arcs = self._arcs(x)
        enc = []
        for i in range(len(t)):
            vrel = (arcs[i][0] if i < len(arcs) else arcs[i - 1][1]) - np.array(st[i][1]) if i == 0 else arcs[i - 1][1] - np.array(st[i][1])
            enc.append({'body': self.names[i], 'date': str(pk.epoch(t[i]))[:10], 'vinf': round(float(np.linalg.norm(vrel)) / 1000, 2)})
        return enc

SEQ = ['earth', 'mars', 'earth', 'mars', 'earth', 'mars', 'earth']     # three Earth-Mars round trips
LEGS = [[100, 450], [300, 1100], [100, 450], [300, 1100], [100, 450], [300, 1100]]

def run(year, islands, gens, pop, evolves):
    t0 = [pk.epoch(f'{year}-01-01').mjd2000, pk.epoch(f'{year + 2}-01-01').mjd2000]
    udp = Cycler(SEQ, t0, LEGS)
    archi = pg.archipelago(n=islands, algo=pg.algorithm(pg.sade(gen=gens, variant=6, variant_adptv=2, memory=True)),
                           prob=pg.problem(udp), pop_size=pop, seed=year)
    for _ in range(evolves):
        archi.evolve(); archi.wait_check()
    p = min((isl.get_population() for isl in archi), key=lambda q: q.champion_f[0])
    enc = udp.report(p.champion_x)
    return {'window': f'{year}-{year + 1}', 'maintenance_kms': round(float(p.champion_f[0]), 3), 'encounters': enc,
            'mean_vinf_earth': round(float(np.mean([e['vinf'] for e in enc if e['body'] == 'earth'])), 2),
            'mean_vinf_mars': round(float(np.mean([e['vinf'] for e in enc if e['body'] == 'mars'])), 2)}

def main():
    first, last = int(sys.argv[1]), int(sys.argv[2])
    budget = dict(islands=int(sys.argv[3]), gens=int(sys.argv[4]), pop=int(sys.argv[5]), evolves=int(sys.argv[6]))
    out = []
    for y in range(first, last + 1, 2):
        t = time.time(); r = run(y, **budget); r['seconds'] = round(time.time() - t, 1)
        out.append(r); print(json.dumps(r), flush=True)
        json.dump(out, open(f'/work/flyby/cyclers_{first}_{last}.json', 'w'), indent=1)

if __name__ == '__main__':      # pygmo starts worker processes that re-import this module
    main()
