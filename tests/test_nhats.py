"""External check: solmaps legs against NASA NHATS round trips (cases fetched by tests/nhats_check.py, committed).

NHATS reports, per trajectory, the excess speed leaving Earth, the relative speeds arriving at and leaving the asteroid,
and the atmospheric entry speed, capped at 12 km/s by a braking burn. Earth departure differs by up to about 0.012 km/s
because NHATS starts from Earth's centre and solmaps from the Earth-Moon barycentre (4,700 km apart).
"""
import json, math, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'flyby'))
import pytest
pk = pytest.importorskip('pykep')
CASES = ROOT / 'tests' / 'nhats_cases.json'
pytestmark = pytest.mark.skipif(not (ROOT / 'flyby' / 'horizons' / '99942.json').exists(), reason='run make ephemeris first')

def test_legs_match_nhats():
    import bodies
    atlas = bodies.load_atlas(); E = bodies.body('earth', atlas); n = 0
    for c in json.load(open(CASES))['cases']:
        y, m, d = map(int, c['launch'].split('-'))
        jd = 367 * y - (7 * (y + (m + 9) // 12)) // 4 + (275 * m) // 9 + d + 1721013.5
        if jd < 2461041.5:
            continue
        A = bodies.body(c['des'], atlas); t0 = jd - 2451544.5
        t1, ts, t2 = float(c['dur_out']), float(c['dur_at']), float(c['dur_ret'])
        r0, v0 = E.eph(t0); r1, v1 = A.eph(t0 + t1)
        l = pk.lambert_problem(r0=r0, r1=r1, tof=t1 * 86400, mu=pk.MU_SUN)
        r2, v2 = A.eph(t0 + t1 + ts); r3, v3 = E.eph(t0 + t1 + ts + t2)
        l2 = pk.lambert_problem(r0=r2, r1=r3, tof=t2 * 86400, mu=pk.MU_SUN)
        assert math.dist(l.v0[0], v0) / 1000 == pytest.approx(float(c['v_dep_earth']), abs=0.02)
        assert math.dist(l.v1[0], v1) / 1000 == pytest.approx(float(c['vrel_arr_neo']), abs=0.005)
        assert math.dist(l2.v0[0], v2) / 1000 == pytest.approx(float(c['vrel_dep_neo']), abs=0.005)
        entry = math.sqrt((math.dist(l2.v1[0], v3) / 1000) ** 2 + 2 * 398600.4 / (6378.1 + 125))
        if float(c['v_arr_earth']) < 11.999:                 # NHATS caps entry at 12 km/s with a braking burn
            assert entry == pytest.approx(float(c['v_arr_earth']), abs=0.01)
        n += 1
    assert n >= 6
