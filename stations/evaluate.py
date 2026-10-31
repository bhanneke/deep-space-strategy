# Station siting, step 2: objectives per candidate orbit and technology state, from stations/costs.npz.
#   reach      mean cheapest cost from the station to the 12 sites (km/s, lower is better)
#   value      share of Earth-launched mass a station saves on trips from low Earth orbit to the sites
#              (hub leverage of the growth map; each trip uses the best station of the set)
#   control    share of (origin, destination, date) trips between sites where a refuelling stop at the
#              station splits the trip into two legs that each cost less than the direct trip
#   supply     cheapest cost of delivering propellant to the station from the sources a technology state allows
# Waiting: every cost is the minimum over departures within the next 120 days. Journeys are chronological: an onward leg
# from a station leaves at or after the ship's arrival there (flight times from precompute.py).
import json, sys
import numpy as np

W = 4                                      # 4 x 30 days of waiting for stations
VE = {'methane/oxygen': 3.7, 'hydrogen/oxygen': 4.4, 'nuclear thermal': 8.5}
SOURCES = [                                # cumulative propellant-source ladder
    ('Earth only', []),
    ('+ lunar ice', ['moon']),
    ('+ water asteroids', ['moon', 'neo']),
    ('+ Mars ice', ['moon', 'neo', 'mars']),
    ('+ Ceres', ['moon', 'neo', 'mars', 'ceres']),
]
LIFT = {'moon': 2.60, 'mars': 4.1, 'ceres': 0.52, 'neo': 0.0}   # surface to depot, km/s, as in the growth map

def wait_min(x, w, axis=-1):
    """minimum over the next w steps (including the current one) along the date axis"""
    out = x.copy()
    for s in range(1, w):
        sh = np.full_like(x, np.nan)
        sl = [slice(None)] * x.ndim
        sl[axis] = slice(0, x.shape[axis] - s); sh[tuple(sl)] = np.take(x, range(s, x.shape[axis]), axis=axis)
        out = np.fmin(out, sh)
    return out

def wait_arg(x, w, axis=-1):
    """minimum over the next w steps and the offset (0..w-1) of the step that attains it"""
    shifts = []
    for k in range(w):
        sh = np.full_like(x, np.nan)
        sl = [slice(None)] * x.ndim
        sl[axis] = slice(0, x.shape[axis] - k); sh[tuple(sl)] = np.take(x, range(k, x.shape[axis]), axis=axis)
        shifts.append(sh)
    st = np.stack(shifts)
    filled = np.where(np.isnan(st), np.inf, st)
    off = filled.argmin(0)
    mn = np.take_along_axis(filled, off[None], 0)[0]
    return np.where(np.isinf(mn), np.nan, mn), off

def main(root):
    d = np.load(f'{root}/stations/costs.npz')
    meta = json.load(open(f'{root}/stations/candidates.json'))
    sites = meta['meta']['sites']; cands = meta['candidates']
    OUT, IN, EML = (wait_min(d[k], W) for k in ('out', 'inn', 'in_eml1'))
    n, ns, nt = OUT.shape
    # arrival date (JD) of the inbound leg actually taken: the cheapest departure within the waiting window plus its flight
    _, off = wait_arg(d['inn'], W)
    t_idx = np.arange(nt)[None, None, :] + off
    tof_in = np.take_along_axis(d['inn_tof'], np.minimum(t_idx, nt - 1), axis=2)
    ARR_IN = d['epoch_jd'] + d['dates'][np.minimum(t_idx, nt - 1)] + tof_in
    net = json.load(open(f'{root}/network/network.json'))
    R = {(r['from'], r['to']): r for r in net['routes']}
    k_of = np.rint((d['epoch_jd'] + d['dates'] - net['meta']['jd0']) / net['meta']['step_days']).astype(int)
    def direct(a, b):                      # site-to-site cost at the station dates, with 120 days of waiting
        r = R[(a, b)]
        full = np.array([np.nan if x is None else x + y for x, y in zip(r['dep'], r['arr'])], np.float32)
        return wait_min(full, 12)[k_of]
    D = np.full((ns, ns, nt), np.nan, np.float32)
    for i, a in enumerate(sites):
        for j, b in enumerate(sites):
            if a != b:
                D[i, j] = direct(a, b)
    e = sites.index('earth')
    dest = [j for j in range(ns) if j != e]

    # fixed candidates: Earth-Moon L1 and low orbits around Venus, Mars and Ceres, on the same footing as the solar orbits
    sys.path.insert(0, f'{root}/stations')
    from precompute import from_eml1, GM_E
    rp = 6378 + 300; vc2, ve2 = GM_E / rp, 2 * GM_E / rp
    def vinf_leo(burn):                    # excess speed from a burn at 300 km (inverse of the Oberth relation)
        return np.sqrt(np.maximum((burn + np.sqrt(vc2)) ** 2 - ve2, 0))
    def series(a, b, f, arrival=False):   # f(dep, arr) per network date, then waiting and the station dates
        r = R[(a, b)]
        dep = np.array([np.nan if x is None else x for x in r['dep']], np.float64)
        arr = np.array([np.nan if x is None else x for x in r['arr']], np.float64)
        val, o = wait_arg(f(dep, arr).astype(np.float32), 12)
        if not arrival:
            return val[k_of]
        kk = np.minimum(k_of + o[k_of], len(r['tof']) - 1)
        tof = np.array([np.nan if r['tof'][k] is None else r['tof'][k] for k in kk], np.float64)
        return val[k_of], net['meta']['jd0'] + kk * net['meta']['step_days'] + tof
    eml_from = lambda dep, arr: np.vectorize(from_eml1)(vinf_leo(dep)) + arr        # L1 -> site
    eml_to = lambda dep, arr: dep + np.vectorize(from_eml1)(vinf_leo(arr))          # site -> L1
    FIX = [('eml1', 'Earth–Moon L1'), ('venus', 'Low Venus orbit'), ('mars', 'Low Mars orbit'), ('ceres', 'Low Ceres orbit')]
    L1_LEO = 3.77                                                                    # low Earth orbit to L1, as in the growth map
    fo = np.full((len(FIX), ns, nt), np.nan, np.float32); fi = np.full_like(fo, np.nan); fe = np.full((len(FIX), nt), np.nan, np.float32)
    station_jd = d['epoch_jd'] + d['dates']
    fa = np.repeat(station_jd[None, None, :], len(FIX), 0).repeat(ns, 1).astype(np.float64)   # inbound arrival JD
    for q, (h, _) in enumerate(FIX):
        for j, s in enumerate(sites):
            if h == 'eml1':
                fo[q, j] = L1_LEO if s == 'earth' else series('earth', s, eml_from)
                if s == 'earth':
                    fi[q, j], fa[q, j] = L1_LEO, station_jd + 4
                else:
                    fi[q, j], fa[q, j] = series(s, 'earth', eml_to, arrival=True)
            else:
                fo[q, j] = 0 if s == h else D[sites.index(h), j]
                if s != h:
                    fi[q, j], fa[q, j] = series(s, h, lambda dep, arr: dep + arr, arrival=True)
                else:
                    fi[q, j] = 0
        fe[q] = 0 if h == 'eml1' else series('earth', h, eml_from)
    OUT, IN, EML = np.concatenate([OUT, fo]), np.concatenate([IN, fi]), np.concatenate([EML, fe])
    ARR_IN = np.concatenate([ARR_IN, fa])
    # onward legs leave after the inbound arrival: index of the first station date at or after it (beyond the grid: none)
    A_IDX = np.ceil((ARR_IN - d['epoch_jd']) / float(d['dates'][1] - d['dates'][0])).astype(np.int64)
    A_OK = (A_IDX < nt) & ~np.isnan(ARR_IN)
    A_IDX = np.clip(np.nan_to_num(A_IDX, nan=nt - 1), 0, nt - 1).astype(np.int64)
    def out_after(i):
        """onward cost OUT[c, j, arrival of the inbound leg from site i at date t] for every candidate, site j, date t"""
        idx = np.broadcast_to(A_IDX[:, i, :][:, None, :], OUT.shape)
        v = np.take_along_axis(OUT, idx, axis=2)
        return np.where(A_OK[:, i, :][:, None, :], v, np.nan)
    cands = cands + [{'name': name, 'site': h, 'a': {'eml1': 1.0, 'venus': 0.723, 'mars': 1.524, 'ceres': 2.77}[h],
                      'e': 0.0, 'i': 0.0, 'node': 0.0, 'w': 0.0, 'M0': 0.0} for h, name in FIX]
    n = len(cands)

    reach = np.nanmean(OUT[:, dest, :], axis=(1, 2))

    OUT_AFTER = [out_after(i) for i in range(ns)]                    # [inbound site][candidate, onward site, date]
    # relay control: max of the two legs below the direct trip, over all ordered site pairs and dates
    ctrl = np.zeros(n)
    for i in range(ns):
        for j in range(ns):
            if i == j:
                continue
            legs = np.fmax(IN[:, i, :], OUT_AFTER[i][:, j, :])
            ctrl += np.nansum(legs < D[i, j][None, :], axis=1)
    ctrl /= ns * (ns - 1) * nt

    # supply cost per source and per ladder step
    src_cost = {
        'earth': np.nanmean(IN[:, e, :], axis=1),
        'moon': LIFT['moon'] + np.nanmean(EML, axis=1),
        'neo': np.nanmean(np.fmin(IN[:, sites.index('101955'), :], IN[:, sites.index('162173'), :]), axis=1),
        'mars': LIFT['mars'] + np.nanmean(IN[:, sites.index('mars'), :], axis=1),
        'ceres': LIFT['ceres'] + np.nanmean(IN[:, sites.index('ceres'), :], axis=1),
    }
    supply = {}
    for label, srcs in SOURCES:
        supply[label] = np.min([src_cost[s] for s in (srcs or ['earth'])], axis=0)

    # fuel saved for trips from low Earth orbit (hub leverage of the growth map), one station at a time
    din = IN[:, e, :][:, None, :]                       # Earth -> station
    Dd = D[e][None, :, :]                               # Earth -> destination, direct
    out_e = OUT_AFTER[e]                                # station -> destination, leaving after the ship has arrived
    ok = (din < Dd) & (out_e < Dd)
    single = {}
    for name, ve in VE.items():
        L = np.nan_to_num(np.where(ok, 1 - np.exp(-(Dd - din) / ve), 0.0)[:, dest, :])
        single[name] = np.round(L.mean(axis=(1, 2)), 4).tolist()

    # refuelling network: every trip between two sites may stop once at a station of the set; a trip costs the
    # longer of its two legs, or the direct trip if that is cheaper. A set of k stations minimises the mean over
    # all site pairs and dates. Pool: the 400 best single relays plus the fixed locations.
    pairs = [(i, j) for i in range(ns) for j in range(ns) if i != j]
    Dp = np.stack([D[i, j] for i, j in pairs])                       # [pair, date]
    valid = ~np.isnan(Dp); base = float(Dp[valid].mean())
    def relay(c):
        m = np.stack([np.fmax(IN[c, i], OUT_AFTER[i][c, j]) for i, j in pairs])
        return np.where(np.isnan(m), np.inf, m)
    gain1 = np.array([base - float(np.fmin(Dp, relay(c))[valid].mean()) for c in range(n)])
    pool = sorted(set(np.argsort(-gain1)[:400].tolist()) | set(range(n - len(FIX), n)))
    M = {c: relay(c) for c in pool}
    def cost(S):
        x = Dp.copy()
        for c in S:
            x = np.fmin(x, M[c])
        return float(x[valid].mean())
    ranked = sorted(pool, key=lambda c: -gain1[c])
    sets, S = [], []
    for k in range(1, 6):                               # greedy, then swaps until no swap improves the set
        S = S + [min((c for c in ranked if c not in S), key=lambda c: cost(S + [c]))]
        improved = True
        while improved:
            improved = False
            for pos in range(len(S)):
                for c in ranked[:150]:
                    if c in S:
                        continue
                    T = S[:pos] + [c] + S[pos + 1:]
                    if cost(T) < cost(S) - 1e-9:
                        S, improved = T, True
        sets.append({'k': k, 'stations': [int(c) for c in S], 'mean_leg': round(cost(S), 3)})
    q = ranked[:50]                                     # QUBO coefficients on the 50 best single relays
    g = {c: base - cost([c]) for c in q}
    qubo = {'candidates': [int(c) for c in q], 'linear': [round(-g[c], 4) for c in q],
            'pairwise_overlap': [[round(g[c] + g[c2] - (base - cost([c, c2])), 4) if c2 > c else 0.0 for c2 in q] for c in q],
            'form': "minimise sum_o a_o x_o + sum_{o<o'} b_{oo'} x_o x_{o'} + lambda (sum_o x_o - k)^2; a = -gain, b = overlap of gains"}
    out = {'meta': {**meta['meta'], 'wait_days': W * 30, 've': VE, 'sources': [s[0] for s in SOURCES], 'lift': LIFT,
                    'objectives': {'reach': 'mean cheapest cost from the station to the 11 sites other than Earth, km/s',
                                   'value': 'share of Earth-launched mass one station saves on trips from low Earth orbit to the other sites, when the station makes its own propellant',
                                   'network': 'mean cost of the longer leg when every trip between two sites may refuel once at a station of the set, km/s; without stations the direct trip',
                                   'control': 'share of site-to-site trips where a refuelling stop at the station gives two legs each cheaper than the direct trip',
                                   'supply': 'cheapest mean cost to deliver propellant to the station from the available sources, km/s'}},
           'candidates': cands,
           'reach': np.round(reach, 3).tolist(), 'control': np.round(ctrl, 4).tolist(),
           'supply': {k: np.round(v, 3).tolist() for k, v in supply.items()},
           'value': single, 'network': {'base_mean': round(base, 3), 'single_gain': np.round(gain1, 4).tolist(), 'sets': sets, 'qubo': qubo}}
    json.dump(out, open(f'{root}/stations/stations.json', 'w'), separators=(',', ':'))
    print('network base', round(base, 2), [(x['k'], x['mean_leg'], [cands[c].get('name') or (cands[c]['a'], cands[c]['e'], cands[c]['i']) for c in x['stations']]) for x in sets])
    print('fuel saved top', [(cands[c].get('name') or cands[c]['a'], single['hydrogen/oxygen'][c]) for c in np.argsort(-np.array(single['hydrogen/oxygen']))[:3]])
    print('reach best', [(cands[c]['a'], cands[c]['e'], round(float(reach[c]), 2)) for c in np.argsort(reach)[:3]])
    print('control best', [(cands[c]['a'], cands[c]['e'], round(float(ctrl[c]), 3)) for c in np.argsort(-ctrl)[:3]])

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '/work')
