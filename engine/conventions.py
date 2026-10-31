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
    'uranus':  {'mu': 5793939,   'r': 436298,     'label': 'circular orbit at Titania distance'},
    'neptune': {'mu': 6836529,   'r': 354800,     'label': 'circular orbit at Triton distance'},
    'pluto':   {'mu': 869.6,     'r': 1188 + 100, 'label': 'low Pluto orbit, 100 km'},
    'ceres':   {'mu': 62.63,     'r': 470 + 50,   'label': 'low Ceres orbit, 50 km'},
}
# SITES: the 12 sites of the pages and of the 10-day files (network.json, alternatives.npz). HELIO: every body with its own
# heliocentric legs (the daily data); MOONS share their planet's legs; PLACES: every node of the journey graph.
SITES = ['earth', 'venus', 'mercury', 'mars', 'ceres', 'vesta', 'psyche', '433', '101955', '162173', 'jupiter', 'saturn']
NEW_SMALL = ['2', '10', '24', '3', '99942', '65803', '25143', '341843', '175706', '4660', '3548', '15094', '11351', '21900', '617']
HELIO = SITES + ['uranus', 'neptune', 'pluto'] + NEW_SMALL
# Moons as end states inside their planet's system. GM (km^3/s^2) and mean radius (km): JPL SSD planetary satellite
# physical parameters; orbit radius a (km) and retrograde motion: JPL SSD satellite mean elements (both 2026-10-31).
# 'alt': parking altitude above the surface (km); None for bodies too small to orbit usefully (rendezvous).
MOONS = {
    'moon':      {'planet': 'earth',   'mu': 4902.800,   'R': 1737.4,  'a': 384400,  'alt': 100,  'retro': False},
    'phobos':    {'planet': 'mars',    'mu': 0.0007087,  'R': 11.08,   'a': 9375,    'alt': None, 'retro': False},
    'deimos':    {'planet': 'mars',    'mu': 0.0000962,  'R': 6.2,     'a': 23457,   'alt': None, 'retro': False},
    'europa':    {'planet': 'jupiter', 'mu': 3202.71210, 'R': 1560.80, 'a': 671100,  'alt': 100,  'retro': False},
    'ganymede':  {'planet': 'jupiter', 'mu': 9887.83275, 'R': 2631.20, 'a': 1070400, 'alt': 100,  'retro': False},
    'callisto':  {'planet': 'jupiter', 'mu': 7179.28340, 'R': 2410.30, 'a': 1882700, 'alt': 100,  'retro': False},
    'titan':     {'planet': 'saturn',  'mu': 8978.13710, 'R': 2574.76, 'a': 1221900, 'alt': 1000, 'retro': False},
    'enceladus': {'planet': 'saturn',  'mu': 7.21037,    'R': 252.10,  'a': 238400,  'alt': 50,   'retro': False},
    'triton':    {'planet': 'neptune', 'mu': 1428.49546, 'R': 1352.60, 'a': 354800,  'alt': 100,  'retro': True},
}
PLACES = HELIO + list(MOONS) + ['eml1']
NAMES = {'earth': 'Earth', 'venus': 'Venus', 'mercury': 'Mercury', 'mars': 'Mars', 'ceres': 'Ceres', 'vesta': 'Vesta',
         'psyche': '16 Psyche', '433': '433 Eros', '101955': '101955 Bennu', '162173': '162173 Ryugu',
         'jupiter': 'Jupiter', 'saturn': 'Saturn', 'eml1': 'Earth–Moon L1',
         'uranus': 'Uranus', 'neptune': 'Neptune', 'pluto': 'Pluto', '2': '2 Pallas', '10': '10 Hygiea', '24': '24 Themis',
         '3': '3 Juno', '99942': '99942 Apophis', '65803': '65803 Didymos', '25143': '25143 Itokawa', '341843': '341843 (2008 EV5)',
         '175706': '175706 (1996 FG3)', '4660': '4660 Nereus', '3548': '3548 Eurybates', '15094': '15094 Polymele',
         '11351': '11351 Leucus', '21900': '21900 Orus', '617': '617 Patroclus',
         'moon': 'Moon', 'phobos': 'Phobos', 'deimos': 'Deimos', 'europa': 'Europa', 'ganymede': 'Ganymede',
         'callisto': 'Callisto', 'titan': 'Titan', 'enceladus': 'Enceladus', 'triton': 'Triton'}
# sites that can supply propellant: Earth by launch, the others from water ice or hydrated minerals (atlas evidence)
# The Moon, Europa, Ganymede, Callisto, Titan and Enceladus: water ice confirmed by spacecraft (LCROSS; Galileo; Cassini).
FUEL_SOURCES = {'earth', 'mars', 'ceres', '101955', '162173', 'eml1', 'moon', 'europa', 'ganymede', 'callisto', 'titan', 'enceladus'}

# Earth-Moon L1 (growth map model): low Earth orbit <-> L1 costs 3.77 km/s; leaving L1 for interplanetary space is the
# cheaper of a burn at L1 and a drop to a 300 km perigee followed by a perigee burn (Oberth).
GM_E, R_L1 = 398600.4, 0.84 * 384400
V_ESC_L1, V_L1 = math.sqrt(2 * GM_E / R_L1), R_L1 * 2 * math.pi / (27.3217 * 86400)
L1_LEO = 3.77

def helio(node):
    """the body whose heliocentric legs a node uses"""
    return 'earth' if node == 'eml1' else MOONS[node]['planet'] if node in MOONS else node

def end_state(node):
    """'parking', 'rendezvous', 'l1' or 'moon', with the parameters the planet-side cost needs"""
    if node == 'eml1':
        return {'kind': 'l1', 'label': 'Earth–Moon L1'}
    if node in MOONS:
        m = MOONS[node]
        return {'kind': 'moon', **m, 'mu_p': PARKING[m['planet']]['mu'],
                'label': f"low {NAMES[node]} orbit, {m['alt']} km" if m['alt'] is not None else f"rendezvous with {NAMES[node]}"}
    if node in PARKING:
        return {'kind': 'parking', **PARKING[node]}
    return {'kind': 'rendezvous', 'label': 'rendezvous (gravity negligible)'}

def planet_side_dv(node, vinf):
    """burn between the node's end state and a hyperbola with excess speed vinf (km/s); same for departure and capture"""
    e = end_state(node)
    if e['kind'] == 'parking':
        return math.sqrt(vinf * vinf + 2 * e['mu'] / e['r']) - math.sqrt(e['mu'] / e['r'])
    if e['kind'] == 'moon':
        # patched conic: the hyperbola reaches the moon's orbit radius tangentially; the ship then meets the moon with
        # the speed difference (sum for a retrograde moon) and captures into low orbit in one Oberth burn
        v_hyp = math.sqrt(vinf * vinf + 2 * e['mu_p'] / e['a'])
        v_moon = math.sqrt(e['mu_p'] / e['a'])
        return moon_capture(e, v_hyp + v_moon if e['retro'] else v_hyp - v_moon)
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
    if mode in ('burn', 'passes') or e['kind'] in ('rendezvous', 'moon'):   # spirals only from parking orbits
        return planet_side_dv(node, vinf), 0.0
    if e['kind'] == 'l1':
        return math.sqrt(vinf * vinf + V_ESC_L1 ** 2) - V_L1, 0.0
    vc = math.sqrt(e['mu'] / e['r'])
    return vc + vinf, vc

def moon_capture(e, v_rel):
    """burn from relative speed v_rel at the moon into its low orbit (Oberth at the parking radius); rendezvous for tiny moons"""
    if e['alt'] is None:
        return abs(v_rel)
    r = e['R'] + e['alt']
    return math.sqrt(v_rel * v_rel + 2 * e['mu'] / r) - math.sqrt(e['mu'] / r)

def local_radius(node):
    """orbit radius around the planet for nodes inside one planet system"""
    return MOONS[node]['a'] if node in MOONS else PARKING[node]['r']

def local_transfer(a, b):
    """transfer between two nodes in the same planet system (moon to moon, parking orbit to moon, L1 to the Moon):
    a Hohmann ellipse around the planet between the two orbit radii, with an Oberth escape or capture at each moon.
    Returns (dv at departure, dv at arrival, days) in km/s and days, or None.
    Earth-Moon L1 to the Moon: the ship leaves L1 with L1's orbital speed and meets the Moon with the difference."""
    if helio(a) != helio(b) or a == b:
        return None
    if 'eml1' in (a, b):
        other = b if a == 'eml1' else a
        if other == 'earth':
            return L1_LEO, 0.0, 4.0
        if other != 'moon':
            return None
        e = end_state('moon')
        cap = moon_capture(e, math.sqrt(e['mu_p'] / e['a']) - V_L1)
        return (0.0, cap, 4.0) if a == 'eml1' else (cap, 0.0, 4.0)
    mu = PARKING[helio(a)]['mu']
    r1, r2 = local_radius(a), local_radius(b)
    at = (r1 + r2) / 2
    v1c, v2c = math.sqrt(mu / r1), math.sqrt(mu / r2)
    v1t, v2t = math.sqrt(mu * (2 / r1 - 1 / at)), math.sqrt(mu * (2 / r2 - 1 / at))
    def end(node, v_orbit, v_transfer):
        if node in MOONS:
            e = end_state(node)
            return moon_capture(e, v_transfer + v_orbit if e['retro'] else v_transfer - v_orbit)
        return abs(v_transfer - v_orbit)
    days = math.pi * math.sqrt(at ** 3 / mu) / 86400
    return end(a, v1c, v1t), end(b, v2c, v2t), days

def planet_side_arr(node, vinf, mode='burn'):
    """planet_side_mode for an array of excess speeds (NaN stays NaN): (total, spiral part) arrays, same formulas"""
    import numpy as np
    v = np.asarray(vinf, dtype=np.float64); e = end_state(node); z = np.zeros_like(v)
    if e['kind'] == 'rendezvous':
        return v.copy(), z
    if e['kind'] == 'moon':
        vh = np.sqrt(v * v + 2 * e['mu_p'] / e['a']); vm = math.sqrt(e['mu_p'] / e['a'])
        rel = vh + vm if e['retro'] else vh - vm
        if e['alt'] is None:
            return np.abs(rel), z
        r = e['R'] + e['alt']
        return np.sqrt(rel * rel + 2 * e['mu'] / r) - math.sqrt(e['mu'] / r), z
    if e['kind'] == 'l1':
        at_l1 = np.sqrt(v * v + V_ESC_L1 ** 2) - V_L1
        if mode == 'spiral':
            return at_l1, z
        rp = 6378 + 300
        v_apo = math.sqrt(2 * GM_E * rp / (R_L1 * (R_L1 + rp)))
        via = (V_L1 - v_apo) + np.sqrt(v * v + 2 * GM_E / rp) - v_apo * R_L1 / rp
        return np.minimum(at_l1, via), z
    vc = math.sqrt(e['mu'] / e['r'])
    if mode == 'spiral':
        return vc + v, np.where(np.isnan(v), np.nan, vc)
    return np.sqrt(v * v + 2 * e['mu'] / e['r']) - vc, z

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
    return {'sites': SITES, 'helio': HELIO, 'moons': MOONS, 'places': PLACES, 'names': NAMES, 'parking': PARKING,
            'fuel_sources': sorted(FUEL_SOURCES), 'l1_leo': L1_LEO}
