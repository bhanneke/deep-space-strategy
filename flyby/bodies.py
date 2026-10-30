# Bodies for solmaps: planets from JPL DE440 (to about 2150), small bodies from JPL SBDB elements.
import json, math
import pykep as pk

PLANETS = ('mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn')
NAIF = {'mercury': 'MERCURY BARYCENTER', 'venus': 'VENUS BARYCENTER', 'earth': 'EARTH BARYCENTER',
        'mars': 'MARS BARYCENTER', 'jupiter': 'JUPITER BARYCENTER', 'saturn': 'SATURN BARYCENTER'}

class DE440Planet:
    """Heliocentric DE440 positions with the physical constants that flybys need.

    The small DE440 kernel shipped with pykep holds planet-system barycentres only. Earth's entry is
    the Earth-Moon barycentre, which lies about 4,700 km from Earth's centre; for heliocentric legs
    that offset is negligible. Mass, radius and safe flyby radius come from pykep's jpl_lp planets.
    """
    def __init__(self, name):
        self.name = name
        self.de = None
        ref = pk.planet(pk.udpla.jpl_lp(name))
        self.mu_self, self.radius, self.safe = ref.get_mu_self(), ref.get_radius(), ref.get_safe_radius()
    def __getstate__(self):
        # pygmo copies problems into worker processes; the SPICE kernel must be loaded again there
        state = dict(self.__dict__); state['de'] = None
        return state
    def eph(self, mjd2000):
        if self.de is None:
            self.de = pk.udpla.de440s(NAIF[self.name], 'ECLIPJ2000', 'SUN')   # loads the kernel in this process
        return self.de.eph(mjd2000)
    def get_name(self):
        return self.name
    def get_mu_central_body(self):
        return pk.MU_SUN
    def get_mu_self(self):
        return self.mu_self
    def get_radius(self):
        return self.radius
    def get_safe_radius(self):
        return self.safe

def small_body(bid, atlas):
    el = {b['id']: b for b in atlas['bodies']}[bid]['orbit']['elements']
    ep = pk.epoch(el['epoch_jd_tdb'] - 2451544.5)
    elem = [el['a_au'] * pk.AU, el['e'], math.radians(el['i_deg']), math.radians(el['node_deg']),
            math.radians(el['arg_peri_deg']), math.radians(el['M_deg'])]
    return pk.planet(pk.udpla.keplerian(ep, elem, pk.MU_SUN, bid, [-1, -1, -1], pk.el_type.KEP_M))

_cache = {}
def body(name, atlas):
    if name not in _cache:
        _cache[name] = pk.planet(DE440Planet(name)) if name in PLANETS else small_body(name, atlas)
    return _cache[name]

def load_atlas(path='/work/atlas/atlas_dataset.json'):
    return json.load(open(path))
