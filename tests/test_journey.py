"""Chronological journeys: legs in order, fuel carried through, deadlines respected, refuelling only where allowed."""
import pathlib, sys, time
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import pytest
from engine.ship import Ship, PRESETS

ALT = ROOT / 'data' / 'alternatives.npz'
pytestmark = pytest.mark.skipif(not ALT.exists(), reason='run make data first')
JD_2027 = 2461406.5

@pytest.fixture(scope='module')
def ctx():
    from engine.data import Alternatives
    from engine.journey import Network
    alts = Alternatives(ALT)
    return alts, Network(alts)

def plan(ctx, ship, a, b, **kw):
    from engine.journey import plan as p
    alts, net = ctx
    return p(alts, ship, a, b, kw.pop('dep_from', JD_2027), kw.pop('dep_until', JD_2027 + 800), net=net, **kw)

def check(res, ship, deadline=None, refuel=None):
    for j in res:
        t = j['legs'][0]['depart_jd']
        for i, leg in enumerate(j['legs']):
            assert leg['depart_jd'] >= t - 1e-6                      # a leg leaves after the previous arrival
            t = leg['arrive_jd']
            assert 0 <= leg['fuel_after_kg'] <= ship.fuel_capacity_kg + 1e-6
            if i < len(j['legs']) - 1 and refuel is not None:
                assert leg['to'] in refuel
        if deadline:
            assert j['arrive_jd'] <= deadline + 1e-6
    arrivals = [j['arrive_jd'] for j in res]
    assert arrivals == sorted(arrivals)

def test_single_leg_matches_the_best_alternative(ctx):
    ship = Ship('big', 10000, 0, 500000, 1e6, 450)
    res = plan(ctx, ship, 'earth', 'mars', max_legs=1, dep_until=JD_2027 + 60)
    alts, _ = ctx
    k0, k1 = alts.k_of(JD_2027), alts.k_of(JD_2027 + 60)
    best = min(sum(alts.priced(o, 'earth', 'mars')) for k in range(k0, k1 + 1) for o in alts.options('earth', 'mars', k))
    assert min(j['dv_total'] for j in res) == pytest.approx(best, abs=1e-3)
    check(res, ship)

def test_fuel_limits_the_journey_and_refuelling_extends_it(ctx):
    small = Ship('small', dry_kg=10000, payload_kg=5000, fuel_capacity_kg=110000, thrust_n=5e5, isp_s=380)   # 7.9 km/s
    direct = plan(ctx, small, 'earth', 'ceres', refuel_at=[], max_legs=1)
    via = plan(ctx, small, 'earth', 'ceres', refuel_at=['mars'], max_legs=2)
    assert small.dv_capacity_kms() < 9.0
    assert not direct                                                 # one tank cannot reach Ceres
    assert via and all(j['stops'] == ['mars'] for j in via)            # a refuel at Mars can
    check(via, small, refuel={'mars'})

def test_deadline_is_respected_and_tightening_it_removes_late_journeys(ctx):
    ship = PRESETS['hydrolox tug']
    loose = plan(ctx, ship, 'earth', 'ceres', arrive_by_jd=JD_2027 + 15 * 365.25)
    assert loose
    deadline = loose[0]['arrive_jd'] + 1
    tight = plan(ctx, ship, 'earth', 'ceres', arrive_by_jd=deadline)
    check(tight, ship, deadline=deadline)
    assert tight and all(j['arrive_jd'] <= deadline + 1e-6 for j in tight)
    early = plan(ctx, ship, 'earth', 'ceres', arrive_by_jd=loose[0]['arrive_jd'] - 30)
    assert all(j['arrive_jd'] <= loose[0]['arrive_jd'] - 30 for j in early)

def test_earliest_and_cheapest_are_both_reported(ctx):
    ship = PRESETS['hydrolox tug']
    res = plan(ctx, ship, 'earth', 'mars', dep_until=JD_2027 + 1200)
    assert len(res) >= 2
    assert res[0]['arrive_jd'] == min(j['arrive_jd'] for j in res)
    assert min(j['fuel_loaded_kg'] for j in res) <= res[0]['fuel_loaded_kg']

def test_l1_departure_saves_fuel_against_low_earth_orbit(ctx):
    ship = PRESETS['hydrolox tug']
    leo = plan(ctx, ship, 'earth', 'mars', max_legs=1, refuel_at=[])
    l1 = plan(ctx, ship, 'eml1', 'mars', max_legs=1, refuel_at=[])
    assert min(j['dv_total'] for j in l1) < min(j['dv_total'] for j in leo)

def test_query_speed(ctx):
    t = time.time()
    plan(ctx, PRESETS['methalox tug'], 'earth', 'psyche', max_legs=3)
    assert time.time() - t < 60

def test_max_payload_uses_the_tightest_segment(ctx):
    from engine.journey import max_payload
    from engine.ship import burn
    import dataclasses
    ship = PRESETS['hydrolox tug']
    j = plan(ctx, ship, 'earth', 'mars', max_legs=1)[0]
    pmax = max_payload(ship, j)
    assert pmax >= ship.payload_kg
    big = dataclasses.replace(ship, payload_kg=pmax)
    leg = j['legs'][0]
    b1 = burn(big, big.fuel_capacity_kg, leg['dv_departure']); b2 = burn(big, b1.fuel_after_kg, leg['dv_arrival'])
    assert b2.feasible and b2.fuel_after_kg < 0.01 * big.fuel_capacity_kg       # the last kilograms of fuel
