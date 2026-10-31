"""Loads the precomputed transfer alternatives and prices them for a journey's end states.

Alternatives: the 10-day, 12-site file (data/alternatives.npz) that the pages use.
DailyAlternatives: the daily data for every heliocentric body (data/daily/<origin>.npz, engine/daily.py), loaded one
origin at a time on first use. Both offer the same interface: index, bodies, jd, k_of, route, options, priced.
"""
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

    def route(self, i):
        """[date, k, (tof, vinf_dep, vinf_arr)] for route index i"""
        return self.alt[i]

    def options(self, a, b, k):
        """transfer alternatives from a to b leaving on date index k, cheapest last (increasing flight time)"""
        i = self.index.get(f'{helio(a)}>{helio(b)}')
        if i is None or not 0 <= k < self.n_dates:
            return []
        rows = self.route(i)[k]
        return [{'tof': float(r[0]), 'vinf_dep': float(r[1]), 'vinf_arr': float(r[2])} for r in rows if not np.isnan(r[0])]

    @staticmethod
    def priced(opt, dep_node, arr_node):
        """planet-side burns for the given end states (a site, or 'eml1' for a departure from or arrival at L1)"""
        return planet_side_dv(dep_node, opt['vinf_dep']), planet_side_dv(arr_node, opt['vinf_arr'])

class DailyAlternatives(Alternatives):
    def __init__(self, folder=ROOT / 'data' / 'daily'):
        self.folder = pathlib.Path(folder)
        self.files = sorted(p for p in self.folder.glob('*.npz') if not p.stem.endswith('.tmp'))
        self.index, self._where, self._loaded = {}, [], {}
        for p in self.files:
            with np.load(p) as d:
                if not hasattr(self, 'jd0'):
                    self.jd0, self.step, self.n_dates = float(d['jd0']), float(d['step']), d['tof'].shape[1]
                for j, b in enumerate(d['targets']):
                    self.index[f'{p.stem}>{b}'] = len(self._where); self._where.append((p.stem, j))
        self.bodies = {p.stem for p in self.files}

    def route(self, i):
        origin, j = self._where[i]
        if origin not in self._loaded:
            with np.load(self.folder / f'{origin}.npz') as d:
                self._loaded[origin] = (d['tof'], d['vinf'])
        tof, vinf = self._loaded[origin]
        return np.concatenate([tof[j][..., None], vinf[j].astype(np.float32)], axis=-1)
