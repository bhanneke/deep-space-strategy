"""Ship mass accounting: the rocket equation, fuel limits and the finite-burn warning."""
import math, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import pytest
from engine.ship import Ship, burn, G0, PRESETS

S = Ship('test', dry_kg=10000, payload_kg=5000, fuel_capacity_kg=30000, thrust_n=1.0e5, isp_s=400)

def test_rocket_equation_and_mass_conservation():
    r = burn(S, 30000, 3.0)
    m0 = 45000.0
    assert r.mass_before_kg == m0
    assert r.mass_after_kg == pytest.approx(m0 * math.exp(-3.0 / (400 * G0 / 1000)))
    assert r.fuel_used_kg + r.fuel_after_kg == pytest.approx(30000)
    assert r.mass_before_kg - r.mass_after_kg == pytest.approx(r.fuel_used_kg)
    assert r.feasible

def test_capacity_is_the_largest_feasible_burn():
    cap = S.dv_capacity_kms()
    assert burn(S, 30000, cap - 1e-6).feasible
    assert not burn(S, 30000, cap + 1e-3).feasible
    assert cap == pytest.approx(400 * G0 / 1000 * math.log(45000 / 15000))

def test_two_burns_equal_one_combined_burn():
    a = burn(S, 30000, 1.5); b = burn(S, a.fuel_after_kg, 1.5)
    c = burn(S, 30000, 3.0)
    assert b.fuel_after_kg == pytest.approx(c.fuel_after_kg)

def test_less_payload_leaves_more_fuel():
    assert burn(S, 30000, 2.0, payload_kg=0).fuel_after_kg > burn(S, 30000, 2.0).fuel_after_kg

def test_burn_duration_and_finite_burn_warning():
    r = burn(S, 30000, 3.0, parking_period_s=5400)
    assert r.duration_s == pytest.approx(r.fuel_used_kg / (1.0e5 / (400 * G0)))
    weak = Ship('weak', 10000, 5000, 30000, thrust_n=2.0e3, isp_s=400)
    assert burn(weak, 30000, 3.0, parking_period_s=5400).warnings          # hours-long burn in a 90-minute orbit
    assert not burn(S, 30000, 0.1, parking_period_s=5400).warnings

def test_presets_are_consistent():
    for p in PRESETS.values():
        assert p.dv_capacity_kms() > 0 and p.ve_kms > 0
    assert PRESETS['solar-electric cargo'].low_thrust and not PRESETS['methalox tug'].low_thrust
