# Convergence of the Sims-Flanagan leg towards continuous thrust: the independent integration's mismatch should fall
# roughly with the square of the number of segments if the gap is discretisation error.
import sys, time
sys.path.insert(0, '/work')
from engine import lowthrust as LT
from engine.ship import PRESETS
ship = PRESETS['solar-electric cargo']
for nseg in (24, 48):
    t = time.time()
    leg = LT.heliocentric('earth', 'mars', ship, ship.mass(ship.fuel_capacity_kg), (10600, 10900), (300, 1200), nseg=nseg, starts=8)
    if not leg:
        print(nseg, 'no feasible leg', round(time.time() - t), 's'); continue
    c = leg['check']
    print(nseg, 'segments | tof', round(leg['tof_days']), 'd | fuel', round(leg['fuel_kg']), 'kg | dv eq', round(leg['dv_equivalent_kms'], 3),
          '| check: pos', round(c['position_error_km']), 'km, vel', round(c['velocity_error_ms'], 1), 'm/s, mass', round(c['mass_error_kg'], 1), 'kg |',
          round(time.time() - t), 's')
