"""More thrust never costs more: an acceleration sweep for the continuous-thrust search (Earth to Mars, 2029).

A stronger ship can fly any profile a weaker ship flies, so the equivalent delta-v must not rise with the acceleration;
and it can fly the two-impulse (Lambert) profile, so it must not exceed the best two-impulse transfer in the same
window. Before the continuation start (lowthrust.heliocentric), the search found 6.70 km/s at 2 mm/s2 against
5.71 km/s at 0.25 mm/s2 and nothing at 3.5 mm/s2.
"""
import pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / 'flyby'))
import numpy as np
import pytest
pk = pytest.importorskip('pykep')
from engine.ship import Ship
from engine import lowthrust as LT

M0, WINDOW, TOF = 40000.0, (10600, 10900), (300, 1200)

def lambert_best():
    import bodies
    atlas = bodies.load_atlas(); E, M = bodies.body('earth', atlas), bodies.body('mars', atlas)
    best = np.inf
    for t0 in np.linspace(*WINDOW, 61):
        r0, v0 = E.eph(float(t0))
        for tof in np.linspace(*TOF, 181):
            r1, v1 = M.eph(float(t0 + tof))
            l = pk.lambert_problem(r0=r0, r1=r1, tof=float(tof) * 86400, mu=LT.MU_SUN)
            best = min(best, np.linalg.norm(np.array(l.v0[0]) - v0) + np.linalg.norm(np.array(v1) - l.v1[0]))
    return best / 1000

@pytest.fixture(scope='module')
def sweep():
    out = {}
    for a0 in (0.25e-3, 1e-3, 3.5e-3):
        ship = Ship('sweep', dry_kg=10000, payload_kg=10000, fuel_capacity_kg=20000, thrust_n=a0 * M0, isp_s=4000)
        leg = LT.heliocentric('earth', 'mars', ship, M0, WINDOW, TOF, nseg=LT.segments(ship, M0), starts=8)
        assert leg is not None, f'no leg at {1000 * a0} mm/s2'
        out[a0] = leg
    return out

def test_more_thrust_never_costs_more(sweep):
    dv = [sweep[a]['dv_equivalent_kms'] for a in sorted(sweep)]
    for weak, strong in zip(dv, dv[1:]):
        assert strong <= weak * 1.01

def test_never_worse_than_two_impulses(sweep):
    bound = lambert_best()
    assert bound == pytest.approx(8.264, abs=0.01)           # regression value of the two-impulse optimum
    for leg in sweep.values():
        assert leg['dv_equivalent_kms'] <= bound * 1.02

def test_every_leg_passes_the_independent_check(sweep):
    for leg in sweep.values():
        assert leg['check']['position_error_km'] < LT.CHECK_KM and abs(leg['check']['mass_error_kg']) < LT.CHECK_KG
