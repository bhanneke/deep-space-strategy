// Search candidate heliocentric depot orbits: which orbits reach many destinations cheaply over 2027-2047?
const fs = require('fs'), E = require('./engine'), { makeNodes } = require('./nodes');
const atlas = require('../atlas/atlas_dataset.json');
const N = Object.fromEntries(makeNodes(atlas).map(n => [n.id, n]));
const DEST = ['earth', 'venus', 'mars', 'ceres', 'vesta', 'psyche', '101955', '433'];
const EPOCH = Date.parse('2027-01-01') / 864e5 + 2440587.5;
const DATES = []; for (let j = EPOCH; j < EPOCH + 20 * 365.25; j += 60) DATES.push(j);
const WAIT = [0, 40, 80, 120];        // a depot can wait up to ~4 months for a better departure
function station(k, label) { return { id: 'stn', name: label, body: { orbit: { kepler: { ...k, epoch: EPOCH } } }, park: null }; }
function tofs(A, B, j) { return E.tofGrid(A, B, j, 12); }
function bestFrom(A, B, j) { let b = Infinity; for (const w of WAIT) { const l = E.bestLeg(A, B, j + w, tofs(A, B, j + w)); if (l && l.dv < b) b = l.dv; } return b; }
function score(k, label) {
  const S = station(k, label), per = {};
  let tot = 0, cnt = 0, earthIn = 0;
  for (const d of DEST) {
    let s = 0;
    for (const j of DATES) s += Math.min(bestFrom(S, N[d], j), 30);
    per[d] = s / DATES.length; tot += per[d]; cnt++;
  }
  for (const j of DATES) earthIn += Math.min(bestFrom(N.earth, S, j), 30);
  return { k, label, reach: tot / cnt, per, earthIn: earthIn / DATES.length };
}
const cands = [];
for (const a of [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])
  for (const e of [0, 0.15, 0.3, 0.45, 0.6])
    for (const w of (e === 0 ? [0] : [0, 90, 180, 270]))
      for (const M0 of [0, 90, 180, 270]) {
        if (a * (1 - e) < 0.45) continue;                       // keep perihelion outside 0.45 AU
        cands.push({ a, e, i: 0, node: 0, w, M0 });
      }
const ref = [
  [{ a: 1.5237, e: 0.0934, i: 1.85, node: 49.56, w: 286.5, M0: 0 }, 'Mars-like orbit (phase 0)'],
  [{ a: 1.6, e: 0.375, i: 0, node: 0, w: 0, M0: 0 }, 'Aldrin-cycler-like (q 1.0, Q 2.2 AU)'],
  [{ a: 1.0, e: 0.0, i: 0, node: 0, w: 0, M0: 0 }, 'Earth-like circular 1 AU']
];
const t0 = Date.now(), res = [];
for (const [k, l] of ref) res.push(score(k, l));
cands.forEach((k, i) => { res.push(score(k, `a ${k.a} e ${k.e} w ${k.w} M0 ${k.M0}`)); if (i % 100 === 0) console.log(i, '/', cands.length, ((Date.now() - t0) / 1000).toFixed(0) + 's'); });
res.sort((x, y) => x.reach - y.reach);
fs.writeFileSync('stations.json', JSON.stringify({ meta: { epoch_jd: EPOCH, dates: DATES.length, step_days: 60, wait_days: WAIT, dest: DEST,
  method: 'Score = average over 20 years (60-day steps) and over 8 destinations of the cheapest impulsive transfer from a depot on the candidate orbit (rendezvous vinf at the depot; parking-orbit burns at planets), allowing up to 120 days of waiting. earthIn = average cost to reach the depot from low Earth orbit. Coplanar candidates (i = 0), phase M0 at 2027-01-01.' }, results: res }));
console.log('done', res.length, 'in', ((Date.now() - t0) / 1000).toFixed(0), 's');
for (const r of res.slice(0, 12)) console.log(r.label.padEnd(40), 'reach', r.reach.toFixed(2), 'earthIn', r.earthIn.toFixed(2));
for (const r of res.filter(r => /like/.test(r.label))) console.log('REF', r.label, 'rank', res.indexOf(r) + 1, 'reach', r.reach.toFixed(2), 'earthIn', r.earthIn.toFixed(2));
