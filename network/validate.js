const E = require('./engine'), { makeNodes } = require('./nodes');
const atlas = require('../atlas/atlas_dataset.json');
const N = Object.fromEntries(makeNodes(atlas).map(n => [n.id, n]));
const jd = d => Date.parse(d) / 864e5 + 2440587.5, ds = j => new Date((j - 2440587.5) * 864e5).toISOString().slice(0, 10);
function scan(a, b, from, to, step = 5) {
  let best = null;
  for (let j = jd(from); j <= jd(to); j += step) { const l = E.bestLeg(N[a], N[b], j, E.tofGrid(N[a], N[b], j, 40)); if (l && (!best || l.dv < best.dv)) best = { ...l, dep: l.dep, date: ds(j) }; }
  return best;
}
const t0 = Date.now();
for (const [a, b, f, t] of [['earth', 'mars', '2026-08-01', '2027-03-01'], ['earth', 'mars', '2028-09-01', '2029-04-01'], ['earth', 'venus', '2027-09-01', '2028-08-01'], ['earth', 'ceres', '2026-11-01', '2029-11-01'], ['earth', 'psyche', '2026-11-01', '2029-11-01'], ['earth', '101955', '2026-11-01', '2030-11-01'], ['mars', 'ceres', '2026-11-01', '2030-11-01'], ['earth', 'jupiter', '2026-11-01', '2028-11-01']]) {
  const r = scan(a, b, f, t);
  console.log(`${a}->${b}: best ${r.date} total ${r.dv.toFixed(2)} km/s (dep ${r.dep.toFixed(2)}, arr ${r.arr.toFixed(2)}, vinf dep ${r.vinfDep.toFixed(2)}) tof ${Math.round(r.tof)} d`);
}
console.log('ms', Date.now() - t0);
