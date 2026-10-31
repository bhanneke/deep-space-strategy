"""End-state conventions: one definition, used the same way by every component (Python and the pages' JavaScript)."""
import math, re, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'flyby')); sys.path.insert(0, str(ROOT / 'stations'))
import pytest
from engine import conventions as C

def test_hohmann_departure_burn_matches_textbook_derivation():
    # Earth-Mars Hohmann: v-infinity 2.95 km/s at Earth; burn from a 300 km orbit 3.59 km/s (technical docs, validation)
    assert C.planet_side_dv('earth', 2.95) == pytest.approx(3.59, abs=0.01)
    # Earth-Jupiter Hohmann: v-infinity 8.79 km/s, burn 6.30 km/s
    assert C.planet_side_dv('earth', 8.79) == pytest.approx(6.30, abs=0.01)

@pytest.mark.parametrize('node', ['earth', 'mars', 'ceres', 'jupiter'])
def test_inverse_round_trip(node):
    for v in (0.5, 3.0, 7.5):
        assert C.vinf_from_planet_side(node, C.planet_side_dv(node, v)) == pytest.approx(v, abs=1e-9)

def test_ceres_is_a_parking_orbit_and_cheaper_than_rendezvous():
    assert C.end_state('ceres')['kind'] == 'parking'
    for v in (1.0, 4.0, 8.0):
        assert C.planet_side_dv('ceres', v) < v          # the Oberth effect, small as it is, helps at Ceres too

def test_small_bodies_are_rendezvous():
    for n in ('vesta', 'psyche', '433', '101955', '162173'):
        assert C.end_state(n)['kind'] == 'rendezvous'
        assert C.planet_side_dv(n, 3.3) == 3.3

def test_l1_departure_is_cheaper_than_low_earth_orbit_at_high_speed_and_never_negative():
    for v in (0.5, 3.0, 6.0, 10.0):
        assert C.planet_side_dv('eml1', v) >= 0
    assert C.planet_side_dv('eml1', 3.0) < C.planet_side_dv('earth', 3.0)

def test_python_components_share_the_conventions():
    pytest.importorskip('pykep')
    import network2, routes
    for node in C.SITES:
        for v in (1.0, 5.0):
            ref = C.planet_side_dv(node, v)
            assert network2.well(v, node) == ref
            assert routes.park_burn(v, node) == ref

def test_javascript_engine_uses_the_same_parking_orbits():
    src = (ROOT / 'network' / 'nodes.js').read_text()
    for node, p in C.PARKING.items():
        m = re.search(node + r":\s*\{\s*mu:\s*([\d.e+]+),\s*r:\s*([\d\s+.]+),", src)
        assert m, f'{node} missing in nodes.js'
        assert float(m.group(1)) == p['mu']
        assert eval(m.group(2)) == p['r']           # radius written as a sum such as 470 + 50
