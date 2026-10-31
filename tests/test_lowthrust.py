"""Low thrust: continuous thrust with changing mass, checked by an independent integration; planet-side phases separate."""
import math, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'flyby'))
import pytest
pk = pytest.importorskip('pykep')
from engine.ship import PRESETS, Ship
from engine import lowthrust as LT

SHIP = PRESETS['solar-electric cargo']
WINDOW, TOF = (10600, 10900), (300, 1200)          # MJD2000: departures Jan-Nov 2029

@pytest.fixture(scope='module')
def leg():
    l = LT.heliocentric('earth', 'mars', SHIP, SHIP.mass(SHIP.fuel_capacity_kg), WINDOW, TOF, nseg=24, starts=8)
    assert l is not None, 'no feasible low-thrust leg'
    return l

def test_independent_integration_reaches_mars(leg):
    c = leg['check']                                # scipy DOP853 on the thrust profile, separate from pykep's integrator
    assert c['position_error_km'] < 50
    assert c['velocity_error_ms'] < 0.05
    assert abs(c['mass_error_kg']) < 1.0

def test_mass_and_throttle_consistency(leg):
    assert leg['fuel_kg'] == pytest.approx(leg['m0_kg'] - leg['mf_kg'])
    assert leg['dv_equivalent_kms'] == pytest.approx(SHIP.ve_kms * math.log(leg['m0_kg'] / leg['mf_kg']))
    assert leg['max_throttle'] <= 1 + 1e-6
    assert leg['vinf_dep_ms'] <= 31 and leg['vinf_arr_ms'] <= 31     # numerical slack only (VINF_SLACK)

def test_regression_earth_mars_2029(leg):
    # 40 t ship, 10 N, Isp 2,500 s, rendezvous with Mars: about 8.2 t of fuel, 5.6-5.7 km/s equivalent (2026-10-31 run)
    assert 7800 < leg['fuel_kg'] < 8700
    assert 5.4 < leg['dv_equivalent_kms'] < 5.9

def test_two_transcriptions_agree_on_fuel(leg):
    """Sims-Flanagan (impulses at segment midpoints) and the zero-order-hold model (continuous thrust) discretise the same
    optimal-control problem; their best fuel use must agree within a few percent."""
    sf = leg['best_sims_flanagan_fuel_kg']
    assert abs(leg['fuel_kg'] - sf) / sf < 0.05

def test_edelbaum_spiral_from_low_earth_orbit():
    dv, days, fuel, m1 = LT.spiral('earth', SHIP, 40000.0)
    assert dv == pytest.approx(math.sqrt(398600.4 / 6678), rel=1e-12)          # circular speed, 7.73 km/s
    assert fuel == pytest.approx(40000 * (1 - math.exp(-dv / SHIP.ve_kms)))
    assert days == pytest.approx(fuel / SHIP.mdot_kgs / 86400)
    assert LT.spiral('psyche', SHIP, 40000.0)[0] == 0.0                         # rendezvous at a small body

def test_full_transfer_keeps_phases_apart_and_carries_mass():
    small = Ship('sep', dry_kg=4000, payload_kg=2000, fuel_capacity_kg=10000, thrust_n=6.0, isp_s=3000, propulsion='electric')
    t = LT.transfer('earth', 'mars', small, WINDOW, TOF, nseg=24, starts=8)
    assert [p['kind'] for p in t['phases']] == ['planet-side', 'heliocentric', 'planet-side']
    assert t['fuel_kg'] == pytest.approx(sum(p['fuel_kg'] for p in t['phases']))
    m = small.mass(small.fuel_capacity_kg)
    assert t['final_mass_kg'] == pytest.approx(m - t['fuel_kg'], rel=1e-9)
    assert t['phases'][0]['dv_kms'] == pytest.approx(7.73, abs=0.01)          # spiral out of low Earth orbit
    if t['feasible']:
        assert t['fuel_kg'] <= small.fuel_capacity_kg

def test_refinement_cache_returns_the_stored_result():
    from engine import refine
    a = refine.impulsive('earth', 'mars', 2461345.5, 300)
    b = refine.impulsive('earth', 'mars', 2461345.5, 300)
    assert b['cached'] and b['dv_total'] == a['dv_total']
    assert a['dv_total'] <= 5.69                                                # at least as good as the 5.68 grid value
