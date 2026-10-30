// Network nodes: atlas bodies plus their parking orbit or arrival convention.
// GM values: NASA NSSDC planetary fact sheets (km^3/s^2). Parking radii are modelling choices, stated here.
const PARK = {
  earth:   { mu: 398600.4, r: 6378 + 300, label: 'low Earth orbit, 300 km' },
  mars:    { mu: 42828.4,  r: 3396 + 300, label: 'low Mars orbit, 300 km' },
  venus:   { mu: 324858.6, r: 6052 + 300, label: 'low Venus orbit, 300 km' },
  mercury: { mu: 22031.8,  r: 2440 + 200, label: 'low Mercury orbit, 200 km' },
  jupiter: { mu: 126686534, r: 1882700,   label: 'circular orbit at Callisto distance' },
  saturn:  { mu: 37931187,  r: 1221900,   label: 'circular orbit at Titan distance' },
  ceres:   { mu: 62.63,     r: 470 + 50,  label: 'low Ceres orbit, 50 km' }
};
const NODE_IDS = ['earth', 'venus', 'mercury', 'mars', 'ceres', 'vesta', 'psyche', '433', '101955', '162173', 'jupiter', 'saturn'];
function makeNodes(atlas) {
  const BY = Object.fromEntries(atlas.bodies.map(b => [b.id, b]));
  return NODE_IDS.map(id => ({ id, name: BY[id].name, body: BY[id], park: PARK[id] || null, parkLabel: PARK[id] ? PARK[id].label : 'rendezvous (gravity negligible)' }));
}
module.exports = { PARK, NODE_IDS, makeNodes };
