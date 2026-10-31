"""On-demand trajectory refinement with a disk cache.

The precomputed alternatives sit on a 10-day date grid with a fixed set of flight times. When a query needs more (an
exact departure date, the best flight time between grid points, or a low-thrust transfer, which takes 10-60 s), the
engine computes it once and stores the result under data/cache/, keyed by a hash of every input. A repeated query
reads the stored result.
"""
import hashlib, json, math, pathlib, sys
sys.path.insert(0, '/work/flyby'); sys.path.insert(0, '/work')
import numpy as np

from engine.conventions import planet_side_dv

ROOT = pathlib.Path(__file__).resolve().parent.parent
CACHE = ROOT / 'data' / 'cache'

def _plain(o):
    if isinstance(o, dict):
        return {k: _plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_plain(v) for v in o]
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o

def cached(kind, params, compute):
    CACHE.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(json.dumps({'kind': kind, **params}, sort_keys=True, default=str).encode()).hexdigest()[:24]
    path = CACHE / f'{kind}-{key}.json'
    if path.exists():
        try:
            res = json.load(open(path)); res['cached'] = True
            return res
        except ValueError:                        # an unreadable entry counts as missing and gets recomputed
            path.unlink()
    res = _plain(compute()); res['params'] = params
    tmp = path.with_suffix('.tmp')                # write, then rename: a crash never leaves a half-written entry
    json.dump(res, open(tmp, 'w'))
    tmp.replace(path)
    res['cached'] = False
    return res

def lowthrust_transfer(origin, target, ship, t0_bounds, tof_bounds, vinf_dep=0.0, nseg=None, starts=8):
    from engine import lowthrust as LT
    params = {'origin': origin, 'target': target, 'ship': ship.to_dict(), 't0_bounds': list(t0_bounds), 'tof_bounds': list(tof_bounds),
              'vinf_dep': vinf_dep, 'nseg': nseg, 'starts': starts, 'model': 'sf+zoh v7, Lambert seeds, auto segments, perigee passes, raw thrust directions, checked legs only, burn or spiral, window at parking-orbit departure, 3 m/s slack'}
    return cached('lowthrust', params, lambda: LT.transfer(origin, target, ship, t0_bounds, tof_bounds, vinf_dep, nseg, starts))

def impulsive(a, b, depart_jd, tof_guess_days, dep_node=None, arr_node=None):
    """exact Lambert transfer at the given departure date, best flight time within 0.6-1.4 times the guess"""
    from scipy.optimize import minimize_scalar
    dep_node, arr_node = dep_node or a, arr_node or b
    def compute():
        import pykep as pk, bodies
        atlas = bodies.load_atlas()
        A, B = bodies.body('earth' if a == 'eml1' else a, atlas), bodies.body('earth' if b == 'eml1' else b, atlas)
        t0 = depart_jd - 2451544.5
        ra, va = A.eph(t0)
        def parts(T):
            rb, vb = B.eph(t0 + T)
            l = pk.lambert_problem(ra, rb, T * pk.DAY2SEC, pk.MU_SUN)
            vd = float(np.linalg.norm(np.array(l.v0[0]) - va) / 1000); vr = float(np.linalg.norm(np.array(l.v1[0]) - vb) / 1000)
            return vd, vr, planet_side_dv(dep_node, vd) + planet_side_dv(arr_node, vr)
        r = minimize_scalar(lambda T: parts(T)[2], bounds=(0.6 * tof_guess_days, 1.4 * tof_guess_days), method='bounded',
                            options={'xatol': 0.01})
        vd, vr, dv = parts(r.x)
        return {'tof_days': float(r.x), 'vinf_dep': vd, 'vinf_arr': vr, 'dv_departure': planet_side_dv(dep_node, vd),
                'dv_arrival': planet_side_dv(arr_node, vr), 'dv_total': dv}
    return cached('impulsive', {'a': a, 'b': b, 'jd': depart_jd, 'tof_guess': tof_guess_days, 'dep': dep_node, 'arr': arr_node}, compute)
