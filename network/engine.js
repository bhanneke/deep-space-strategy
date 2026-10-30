// Temporal transfer network for the solar-system atlas.
// Works in Node (module.exports) and in the browser (window.NetEngine).
// Units: AU, days, km/s for reported delta-v.
(function (root) {
  const K2 = 2.9591220828e-4;          // Sun GM, AU^3/day^2 (Gaussian constant squared)
  const GAUSS_DEG_DAY = 0.9856076686;  // mean motion of a 1 AU orbit, deg/day
  const AUD = 1731.456837;             // km/s per AU/day
  const D2R = Math.PI / 180;

  /* ---------- ephemeris from the atlas dataset ---------- */
  function elementsAt(b, jd) {
    const o = b.orbit || {};
    if (o.j2000_elements) {                    // JPL approximate positions, table 1
      const T = (jd - 2451545.0) / 36525, el = o.j2000_elements, r = o.rates_per_julian_century;
      const a = el.a_au + r.a_au_per_cy * T, e = el.e + r.e_per_cy * T, I = el.I_deg + r.I_deg_per_cy * T;
      const L = el.L_deg + r.L_deg_per_cy * T, wb = el.long_peri_deg + r.long_peri_deg_per_cy * T, Om = el.long_node_deg + r.long_node_deg_per_cy * T;
      const M = ((L - wb) % 360 + 540) % 360 - 180;
      return { a, e, I: I * D2R, node: Om * D2R, w: (wb - Om) * D2R, M: M * D2R };
    }
    if (o.elements && o.elements.a_au && o.elements.epoch_jd_tdb) {   // JPL SBDB osculating elements, two-body propagation
      const el = o.elements, n = GAUSS_DEG_DAY / Math.pow(el.a_au, 1.5);
      const M = ((el.M_deg + n * (jd - el.epoch_jd_tdb)) % 360 + 360) % 360;
      return { a: el.a_au, e: el.e, I: el.i_deg * D2R, node: el.node_deg * D2R, w: el.arg_peri_deg * D2R, M: M * D2R };
    }
    if (o.kepler) return orbitFromKepler(o.kepler, jd);   // synthetic station orbits
    return null;
  }
  function orbitFromKepler(k, jd) {
    const n = GAUSS_DEG_DAY / Math.pow(k.a, 1.5);
    const M = ((k.M0 + n * (jd - k.epoch)) % 360 + 360) % 360;
    return { a: k.a, e: k.e, I: (k.i || 0) * D2R, node: (k.node || 0) * D2R, w: (k.w || 0) * D2R, M: M * D2R };
  }
  function keplerE(M, e) { let E = e < 0.8 ? M : Math.PI; for (let k = 0; k < 40; k++) { const d = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E)); E -= d; if (Math.abs(d) < 1e-13) break; } return E; }
  function rot(k, xp, yp) {
    const cw = Math.cos(k.w), sw = Math.sin(k.w), cO = Math.cos(k.node), sO = Math.sin(k.node), cI = Math.cos(k.I), sI = Math.sin(k.I);
    return [(cw * cO - sw * sO * cI) * xp + (-sw * cO - cw * sO * cI) * yp,
            (cw * sO + sw * cO * cI) * xp + (-sw * sO + cw * cO * cI) * yp,
            (sw * sI) * xp + (cw * sI) * yp];
  }
  function stateFromElements(k) {
    const E = keplerE(k.M, k.e), r = k.a * (1 - k.e * Math.cos(E)), sq = Math.sqrt(1 - k.e * k.e);
    const pos = rot(k, k.a * (Math.cos(E) - k.e), k.a * sq * Math.sin(E));
    const f = Math.sqrt(K2 * k.a) / r;
    const vel = rot(k, -f * Math.sin(E), f * sq * Math.cos(E));
    return { r: pos, v: vel };
  }
  function stateAt(b, jd) { const k = elementsAt(b, jd); return k ? stateFromElements(k) : null; }

  /* ---------- vectors ---------- */
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const len = a => Math.sqrt(dot(a, a));

  /* ---------- Lambert, universal variables, single revolution, prograde ---------- */
  function stumpC(z) { if (z > 1e-6) { const s = Math.sqrt(z); return (1 - Math.cos(s)) / z; } if (z < -1e-6) { const s = Math.sqrt(-z); return (Math.cosh(s) - 1) / (-z); } return 0.5; }
  function stumpS(z) { if (z > 1e-6) { const s = Math.sqrt(z); return (s - Math.sin(s)) / (s * s * s); } if (z < -1e-6) { const s = Math.sqrt(-z); return (Math.sinh(s) - s) / (s * s * s); } return 1 / 6; }
  function lambert(r1, r2, tof, mu) {
    const R1 = len(r1), R2 = len(r2);
    let dnu = Math.acos(Math.max(-1, Math.min(1, dot(r1, r2) / (R1 * R2))));
    if (cross(r1, r2)[2] < 0) dnu = 2 * Math.PI - dnu;            // prograde about the ecliptic pole
    if (Math.abs(Math.sin(dnu)) < 1e-3) return null;               // near 0/180 deg: plane undefined
    const A = Math.sin(dnu) * Math.sqrt(R1 * R2 / (1 - Math.cos(dnu)));
    const yf = z => R1 + R2 + A * (z * stumpS(z) - 1) / Math.sqrt(stumpC(z));
    const sq = Math.sqrt(mu);
    const F = z => { const y = yf(z); return Math.pow(y / stumpC(z), 1.5) * stumpS(z) + A * Math.sqrt(y) - sq * tof; };
    let lo = -60, hi = 4 * Math.PI * Math.PI - 1e-4;
    if (yf(lo) < 0) { let a = lo, b = hi; if (yf(b) < 0) return null; for (let i = 0; i < 60; i++) { const m = (a + b) / 2; if (yf(m) < 0) a = m; else b = m; } lo = b; }
    if (!(F(lo) <= 0 && F(hi) >= 0)) return null;
    for (let i = 0; i < 60; i++) { const m = (lo + hi) / 2; if (F(m) < 0) lo = m; else hi = m; }
    const z = (lo + hi) / 2, y = yf(z), f = 1 - y / R1, g = A * Math.sqrt(y / mu), gd = 1 - y / R2;
    return { v1: mul(sub(r2, mul(r1, f)), 1 / g), v2: mul(sub(mul(r2, gd), r1), 1 / g) };
  }

  /* ---------- node costs: hyperbolic excess to or from a parking orbit ---------- */
  // park: {mu km^3/s^2, r km}. Burn = sqrt(vinf^2 + 2mu/r) - sqrt(mu/r) (Oberth burn at periapsis).
  function wellCost(vinf, park) { if (!park) return vinf; return Math.sqrt(vinf * vinf + 2 * park.mu / park.r) - Math.sqrt(park.mu / park.r); }

  function leg(A, B, jdDep, tof) {
    const sA = stateAt(A.body, jdDep), sB = stateAt(B.body, jdDep + tof);
    if (!sA || !sB) return null;
    const L = lambert(sA.r, sB.r, tof, K2);
    if (!L || !isFinite(L.v1[0])) return null;
    const vD = len(sub(L.v1, sA.v)) * AUD, vA = len(sub(L.v2, sB.v)) * AUD;
    const dep = wellCost(vD, A.park), arr = B.arrive === 'flyby' ? 0 : wellCost(vA, B.park);
    return { dv: dep + arr, dep, arr, vinfDep: vD, vinfArr: vA, tof };
  }
  function hohmannDays(a1, a2) { const at = (a1 + a2) / 2; return Math.PI * Math.sqrt(at * at * at / K2); }
  function tofGrid(A, B, jd, steps) {
    const a1 = elementsAt(A.body, jd).a, a2 = elementsAt(B.body, jd).a, tH = hohmannDays(a1, a2);
    const lo = Math.max(25, 0.3 * tH), hi = Math.min(3200, 1.8 * tH), out = [];
    for (let i = 0; i < steps; i++) out.push(lo * Math.pow(hi / lo, i / (steps - 1)));
    return out;
  }
  function bestLeg(A, B, jdDep, tofs) {
    let best = null;
    for (const t of tofs) { const l = leg(A, B, jdDep, t); if (l && (!best || l.dv < best.dv)) best = l; }
    return best;
  }

  const api = { K2, AUD, elementsAt, stateAt, lambert, wellCost, leg, bestLeg, tofGrid, hohmannDays, stateFromElements, orbitFromKepler };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.NetEngine = api;
})(typeof window !== 'undefined' ? window : globalThis);
