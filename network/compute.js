// Precompute the temporal network: best transfer per route for each departure date.
const fs = require('fs'), E = require('./engine'), { makeNodes } = require('./nodes');
const atlas = require('../atlas/atlas_dataset.json');
const nodes = makeNodes(atlas);
const jd0 = Date.parse('2026-11-01') / 864e5 + 2440587.5, jd1 = Date.parse('2048-12-31') / 864e5 + 2440587.5, STEP = 10, TOFS = 30;
const dates = []; for (let j = jd0; j <= jd1; j += STEP) dates.push(j);
const routes = [];
const t0 = Date.now();
for (const A of nodes) for (const B of nodes) {
  if (A === B) continue;
  const dep = [], arr = [], tof = [];
  for (const j of dates) {
    const l = E.bestLeg(A, B, j, E.tofGrid(A, B, j, TOFS));
    dep.push(l ? +l.dep.toFixed(2) : null); arr.push(l ? +l.arr.toFixed(2) : null); tof.push(l ? Math.round(l.tof) : null);
  }
  routes.push({ from: A.id, to: B.id, dep, arr, tof });
}
const out = {
  meta: { built: new Date().toISOString().slice(0, 10), jd0, step_days: STEP, n: dates.length, tof_steps: TOFS,
    method: 'Single-revolution prograde Lambert transfers between JPL-element positions (atlas dataset); best of a log-spaced flight-time grid (0.3 to 1.8 x Hohmann time, 25 to 3200 days). Departure and arrival burns from the stated parking orbits via the Oberth relation sqrt(vinf^2 + 2mu/r) - sqrt(mu/r); small bodies use rendezvous vinf. Impulsive burns, no gravity assists, no low-thrust.' },
  nodes: nodes.map(n => ({ id: n.id, name: n.name, park: n.parkLabel })),
  routes
};
fs.writeFileSync('network.json', JSON.stringify(out));
console.log('routes', routes.length, 'dates', dates.length, 'seconds', ((Date.now() - t0) / 1000).toFixed(1), 'bytes', fs.statSync('network.json').size);
