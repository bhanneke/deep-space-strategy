"""Low-thrust atlas (precompute): the electric preset on every ordered site pair and every two-year departure window,
2027-2098, through the refinement cache. The comparison page then answers these queries from the cache, and
data/lowthrust_atlas.json summarises fuel, flight time and the independent check for every pair and window.

Flight-time bounds per pair: 0.5 to 3 Hohmann times (at most 4,000 days). Resumable: finished entries come from the
cache. Run: python -m engine.lt_atlas  (in the solmaps-env container, from /work)
"""
import json, math, os, sys, time
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work')

from engine.conventions import SITES
from engine.ship import PRESETS
from engine import refine

K2 = 2.9591220828e-4
A = {'earth': 1.0, 'venus': 0.723, 'mercury': 0.387, 'mars': 1.524, 'ceres': 2.77, 'vesta': 2.36, 'psyche': 2.92,
     '433': 1.458, '101955': 1.126, '162173': 1.191, 'jupiter': 5.20, 'saturn': 9.54}

def mjd(year):
    return (year - 2000) * 365.25 - 0.5

def main(out='/work/data/lowthrust_atlas.json'):
    if os.environ.get('SOLMAPS_POOL'):
        import pygmo as pg
        pg.mp_island.init_pool(int(os.environ['SOLMAPS_POOL']))
    ship = PRESETS['solar-electric cargo']
    rows = json.load(open(out)) if os.path.exists(out) else []
    done = {(r['origin'], r['target'], r['window']) for r in rows}
    for a in SITES:
        for b in SITES:
            if a == b:
                continue
            tH = math.pi * math.sqrt(((A[a] + A[b]) / 2) ** 3 / K2)
            tof = (max(60, 0.5 * tH), min(4000, 3 * tH))
            for y in range(2027, 2098, 2):
                key = (a, b, f'{y}-{y + 1}')
                if key in done:
                    continue
                t = time.time()
                r = refine.lowthrust_transfer(a, b, ship, (mjd(y), mjd(y + 2)), tof)
                leg = r.get('leg') or {}
                row = {'origin': a, 'target': b, 'window': key[2], 'feasible': bool(r.get('feasible')), 'reason': r.get('reason', ''),
                       'fuel_kg': r.get('fuel_kg'), 'days': (r['arrive_mjd2000'] - r['depart_mjd2000']) if leg else None,
                       'depart_mjd2000': r.get('depart_mjd2000'), 'check_km': leg.get('check', {}).get('position_error_km'),
                       'seconds': round(time.time() - t, 1)}
                rows.append(row)
                print(json.dumps(row), flush=True)
                json.dump(rows, open(out + '.tmp', 'w')); os.replace(out + '.tmp', out)

if __name__ == '__main__':
    main()
