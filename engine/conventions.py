"""How a ship enters and leaves each node: the single source of truth for the engine, the network and the flyby search.

Every transfer has two parts that the engine keeps apart:
  - the heliocentric leg between two bodies, described by the excess speeds (v-infinity) at both ends; and
  - the planet-side manoeuvres: the burn that turns a parking orbit into the departure hyperbola, and the burn that
    captures from the arrival hyperbola into a parking orbit. At bodies with negligible gravity the ship matches the
    body's velocity (rendezvous), so the planet-side cost equals the excess speed.

Earlier versions priced Ceres two ways (the flyby search as a rendezvous, the network as capture into a 50 km orbit);
all code now takes the end state from here.
"""
import math

# GM in km^3/s^2 (NASA NSSDC fact sheets) and parking-orbit radius in km
PARKING = {
    'earth':   {'mu': 398600.4,  'r': 6378 + 300, 'label': 'low Earth orbit, 300 km'},
    'venus':   {'mu': 324858.6,  'r': 6052 + 300, 'label': 'low Venus orbit, 300 km'},
    'mercury': {'mu': 22031.8,   'r': 2440 + 200, 'label': 'low Mercury orbit, 200 km'},
    'mars':    {'mu': 42828.4,   'r': 3396 + 300, 'label': 'low Mars orbit, 300 km'},
    'jupiter': {'mu': 126686534, 'r': 1882700,    'label': 'circular orbit at Callisto distance'},
    'saturn':  {'mu': 37931187,  'r': 1221900,    'label': 'circular orbit at Titan distance'},
    'ceres':   {'mu': 62.63,     'r': 470 + 50,   'label': 'low Ceres orbit, 50 km'},
}
SITES = ['earth', 'venus', 'mercury', 'mars', 'ceres', 'vesta', 'psyche', '433', '101955', '162173', 'jupiter', 'saturn']
NAMES = {'earth': 'Earth', 'venus': 'Venus', 'mercury': 'Mercury', 'mars': 'Mars', 'ceres': 'Ceres', 'vesta': 'Vesta',
         'psyche': '16 Psyche', '433': '433 Eros', '101955': '101955 Bennu', '162173': '162173 Ryugu',
         'jupiter': 'Jupiter', 'saturn': 'Saturn', 'eml1': 'Earth–Moon L1'}
# sites that can supply propellant: Earth by launch, the others from water ice or hydrated minerals (atlas evidence)
FUEL_SOURCES = {'earth', 'mars', 'ceres', '101955', '162173', 'eml1'}

# Earth-Moon L1 (growth map model): low Earth orbit <-> L1 costs 3.77 km/s; leaving L1 for interplanetary space is the
# cheaper of a burn at L1 and a drop to a 300 km perigee followed by a perigee burn (Oberth).
GM_E, R_L1 = 398600.4, 0.84 * 384400
V_ESC_L1, V_L1 = math.sqrt(2 * GM_E / R_L1), R_L1 * 2 * math.pi / (27.3217 * 86400)
L1_LEO = 3.77

def end_state(node):
    """'parking', 'rendezvous' or 'l1', with the parameters the planet-side cost needs"""
    if node == 'eml1':
        return {'kind': 'l1', 'label': 'Earth–Moon L1'}
    if node in PARKING:
        return {'kind': 'parking', **PARKING[node]}
    return {'kind': 'rendezvous', 'label': 'rendezvous (gravity negligible)'}

def planet_side_dv(node, vinf):
    """burn between the node's end state and a hyperbola with excess speed vinf (km/s); same for departure and capture"""
    e = end_state(node)
    if e['kind'] == 'parking':
        return math.sqrt(vinf * vinf + 2 * e['mu'] / e['r']) - math.sqrt(e['mu'] / e['r'])
    if e['kind'] == 'l1':
        rp = 6378 + 300
        at_l1 = math.sqrt(vinf * vinf + V_ESC_L1 ** 2) - V_L1
        v_apo = math.sqrt(2 * GM_E * rp / (R_L1 * (R_L1 + rp)))
        via_perigee = (V_L1 - v_apo) + math.sqrt(vinf * vinf + 2 * GM_E / rp) - v_apo * R_L1 / rp
        return min(at_l1, via_perigee)
    return vinf

PASS_ARC = 0.1        # one burn may last this share of the parking-orbit period (otherwise split into perigee passes)
MAX_PASSES = 50       # beyond this many perigee passes the ship spirals instead

def planet_side_plan(node, mass_kg, ve_kms, mdot_kgs):
    """how a ship of this mass and engine leaves or reaches the node's parking orbit: ('burn', 1), ('passes', n) or
    ('spiral', 0). The escape burn's duration decides: up to PASS_ARC of the orbit period it is one Oberth burn;
    longer burns split into n perigee burns of at most that length (apogee raising), which keep the Oberth effect at the
    cost of the time on the intermediate orbits; beyond MAX_PASSES the ship spirals (Edelbaum)."""
    e = end_state('earth' if node == 'eml1' else node)
    if e['kind'] != 'parking':
        return 'burn', 1
    dv = math.sqrt(2 * e['mu'] / e['r']) - math.sqrt(e['mu'] / e['r'])
    seconds = mass_kg * (1 - math.exp(-dv / ve_kms)) / mdot_kgs
    n = max(1, math.ceil(seconds / (PASS_ARC * 2 * math.pi * math.sqrt(e['r'] ** 3 / e['mu'])) - 1e-12))
    return ('burn', 1) if n == 1 else ('passes', n) if n <= MAX_PASSES else ('spiral', 0)

def pass_days(node, n):
    """days on the intermediate orbits when a planet-side burn splits into n perigee passes. The first n-1 burns raise
    the perigee speed in equal steps from circular toward escape, each followed by one full orbit; the last burn turns
    the orbit into the hyperbola (and carries any excess speed). Every intermediate orbit is bound."""
    e = end_state('earth' if node == 'eml1' else node)
    if n <= 1 or e['kind'] != 'parking':
        return 0.0
    mu, r = e['mu'], e['r']
    vc, vesc = math.sqrt(mu / r), math.sqrt(2 * mu / r)
    total = 0.0
    for k in range(1, n):
        v = vc + k * (vesc - vc) / n
        a = 1.0 / (2.0 / r - v * v / mu)
        total += 2 * math.pi * math.sqrt(a ** 3 / mu)
    return total / 86400.0

def planet_side_mode(node, vinf, mode):
    """planet-side cost for a ship flying in the given mode: (total dv, spiral part) in km/s.
    'burn' or 'passes': Oberth burns (planet_side_dv), in one go or split over perigee passes. 'spiral': a slow tangential spiral between the parking orbit and
    escape (Edelbaum: the circular speed), then an impulse of vinf at the edge of the sphere of influence; at L1 only the
    burn at L1 itself, since the perigee drop needs a short burn near Earth."""
    e = end_state(node)
    if mode in ('burn', 'passes') or e['kind'] == 'rendezvous':
        return planet_side_dv(node, vinf), 0.0
    if e['kind'] == 'l1':
        return math.sqrt(vinf * vinf + V_ESC_L1 ** 2) - V_L1, 0.0
    vc = math.sqrt(e['mu'] / e['r'])
    return vc + vinf, vc

def vinf_from_planet_side(node, dv):
    """inverse of planet_side_dv for parking orbits and rendezvous"""
    e = end_state(node)
    if e['kind'] == 'parking':
        vc = math.sqrt(e['mu'] / e['r'])
        return math.sqrt(max((dv + vc) ** 2 - 2 * e['mu'] / e['r'], 0.0))
    if e['kind'] == 'rendezvous':
        return dv
    raise ValueError('no closed-form inverse for ' + e['kind'])

def parking_period_s(node):
    e = end_state(node)
    return 2 * math.pi * math.sqrt(e['r'] ** 3 / e['mu']) if e['kind'] == 'parking' else None

def export():
    """conventions as JSON for the pages, so JavaScript uses the same numbers"""
    return {'sites': SITES, 'names': NAMES, 'parking': PARKING, 'fuel_sources': sorted(FUEL_SOURCES), 'l1_leo': L1_LEO}
