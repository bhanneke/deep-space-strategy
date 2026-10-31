# Sampled flyby trajectories for the pages: rebuilds each optimised route from its stored decision vector with
# pykep's own ephemeris function (deep-space manoeuvres included) and samples the spacecraft position.
# Check: the sampled position at arrival must sit at the target planet. Output: flyby/paths_<file>.json
import json, sys, glob, os
sys.path.insert(0, '/work/flyby')
import numpy as np
import pykep as pk
import routes, bodies

NS = 80

def main():
    for f in sorted(glob.glob('/work/flyby/routes2_*.json')):
        rows = [r for r in json.load(open(f)) if r.get('feasible') and len(r['seq']) > 2 and 'x' in r]
        if not rows:
            continue
        out, worst = {}, 0.0
        for r in rows:
            y = int(r['bin'][:4])
            t0 = [pk.epoch(f'{y}-01-01').mjd2000, pk.epoch(f'{y + 2}-01-01').mjd2000]
            prob = routes.Route(r['seq'], t0).inner
            eph = prob.get_eph_function(r['x'])
            x = r['x']; start = x[0]; end = start + x[5] + sum(x[9::4])          # exact epochs from the decision vector
            enc = [start]
            for T in [x[5]] + list(x[9::4]):
                enc.append(enc[-1] + T)
            ts = np.unique(np.concatenate([np.linspace(start + 1e-6, end - 1e-6, NS), np.array(enc[1:-1])]))   # flyby moments sampled exactly
            pts = [[round(float(t), 2)] + [round(float(v) / pk.AU, 5) for v in eph(float(t))[0]] for t in ts]
            tgt = bodies.body(r['target'], routes.ATLAS).eph(float(end))[0]
            miss = float(np.linalg.norm(np.array(pts[-1][1:]) - np.array(tgt) / pk.AU))
            worst = max(worst, miss)
            out[f"{r['origin']}|{r['target']}|{r['launch_mjd2000']}|{'-'.join(r['seq'])}"] = pts
        name = '/work/flyby/paths_' + os.path.basename(f)[len('routes2_'):]
        json.dump(out, open(name, 'w'), separators=(',', ':'))
        print(os.path.basename(name), len(out), 'paths; largest arrival miss', round(worst, 6), 'AU')

if __name__ == '__main__':
    main()
