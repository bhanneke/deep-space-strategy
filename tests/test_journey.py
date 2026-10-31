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
    ship = Ship('big', 10000, 0, 500000, 1e8, 450)        # thrust high enough for single burns: pure Lambert costs
    res = plan(ctx, ship, 'earth', 'mars', max_legs=1, dep_until=JD_2027 + 60)
    alts, _ = ctx
    import math
    k0 = math.ceil((JD_2027 - alts.jd0) / alts.step); k1 = math.floor((JD_2027 + 60 - alts.jd0) / alts.step)
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

def test_departures_stay_inside_the_window(ctx):
    for start in (JD_2027, JD_2027 + 3, JD_2027 + 7):
        res = plan(ctx, PRESETS['methalox tug'], 'earth', 'mars', dep_from=start, dep_until=start + 400)
        assert res and all(start - 1e-6 <= j['depart_jd'] <= start + 400 + 1e-6 for j in res)

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

def test_planet_side_modes_follow_burn_length(ctx):
    """one burn, perigee passes or a spiral, from the escape burn's length against the parking-orbit period"""
    from engine.journey import node_mode
    assert node_mode('earth', Ship('fast', 10000, 0, 500000, 1e8, 450)) == ('burn', 1)
    mode, n = node_mode('earth', Ship('1 MN', 10000, 0, 500000, 1e6, 450))     # 20 min, a fifth of the period
    assert mode == 'passes' and 2 <= n <= 3
    assert node_mode('earth', PRESETS['solar-electric cargo'])[0] == 'spiral'
    assert node_mode('101955', PRESETS['solar-electric cargo']) == ('burn', 1)        # rendezvous: no parking orbit

def test_perigee_passes_keep_the_oberth_cost_and_take_time(ctx):
    """passes cost the same delta-v as one burn and start inside the window: the 1 MN ship over the window matches the
    single-burn ship over the window opened later by the pass time"""
    from engine.conventions import pass_days
    fast, slow = Ship('fast', 10000, 0, 500000, 1e8, 450), Ship('1 MN', 10000, 0, 500000, 1e6, 450)
    pd = pass_days('earth', 3)
    assert 0.1 < pd < 3                                 # two intermediate orbits, about 9 hours
    b = plan(ctx, slow, 'earth', 'mars', max_legs=1, dep_until=JD_2027 + 400, max_results=500)
    a = plan(ctx, fast, 'earth', 'mars', max_legs=1, dep_from=JD_2027 + pd, dep_until=JD_2027 + 400 + pd, max_results=500)
    assert min(j['dv_total'] for j in b) == pytest.approx(min(j['dv_total'] for j in a), abs=1e-3)
    for j in b:
        leg = j['legs'][0]
        assert leg['modes'][0] == 'passes' and leg['pass_days'][0] == pytest.approx(pd, abs=0.01)
        assert leg['tof_days'] == pytest.approx(leg['heliocentric_days'] + leg['pass_days'][0] + leg['pass_days'][1], abs=0.15)
    check(b, slow)

def test_spiralling_ship_obeys_the_burn_share(ctx):
    """a medium-thrust ship spirals out, then makes v-infinity impulses that must stay within 5% of the flight"""
    from engine.journey import IMPULSIVE_SHARE
    ship = Ship('2 kN, Isp 10,000 s', dry_kg=40000, payload_kg=20000, fuel_capacity_kg=40000, thrust_n=2000.0, isp_s=10000)
    res = plan(ctx, ship, 'earth', 'mars', max_legs=1, dep_until=JD_2027 + 800)
    assert res
    for j in res:
        leg = j['legs'][0]
        assert leg['modes'] == ['spiral', 'spiral'] and leg['heliocentric_burn_share'] <= IMPULSIVE_SHARE + 1e-9
        assert leg['depart_jd'] >= JD_2027 - 1e-6
    check(res, ship)
    assert not plan(ctx, PRESETS['solar-electric cargo'], 'earth', 'mars', max_legs=1, dep_until=JD_2027 + 800)
