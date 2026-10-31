"""Loads the precomputed transfer alternatives and prices them for a journey's end states."""
import pathlib
import numpy as np
from engine.conventions import planet_side_dv, helio

ROOT = pathlib.Path(__file__).resolve().parent.parent

class Alternatives:
    def __init__(self, path=ROOT / 'data' / 'alternatives.npz'):
        d = np.load(path)
        self.alt = d['alt']                               # [route, date, k, (tof days, vinf_dep, vinf_arr)] in km/s
        self.index = {p: i for i, p in enumerate(d['pairs'])}
        self.bodies = {p.split('>')[0] for p in self.index}
        self.jd0, self.step = float(d['jd0']), float(d['step'])
        self.n_dates = self.alt.shape[1]

    def jd(self, k):
        return self.jd0 + k * self.step

    def k_of(self, jd):
        return int(round((jd - self.jd0) / self.step))

    def options(self, a, b, k):
        """transfer alternatives from a to b leaving on date index k, cheapest last (increasing flight time)"""
        i = self.index.get(f'{helio(a)}>{helio(b)}')
        if i is None or not 0 <= k < self.n_dates:
            return []
        rows = self.alt[i, k]
        return [{'tof': float(r[0]), 'vinf_dep': float(r[1]), 'vinf_arr': float(r[2])} for r in rows if not np.isnan(r[0])]

    @staticmethod
    def priced(opt, dep_node, arr_node):
        """planet-side burns for the given end states (a site, or 'eml1' for a departure from or arrival at L1)"""
        return planet_side_dv(dep_node, opt['vinf_dep']), planet_side_dv(arr_node, opt['vinf_arr'])
