# Low-thrust check numbers for the documentation: Earth-Mars, departures Jan-Nov 2029, 24 and 48 segments.
import sys, time
sys.path.insert(0, '/work')
from engine import lowthrust as LT
from engine.ship import Ship
ship = Ship('sep test', dry_kg=8000, payload_kg=20000, fuel_capacity_kg=12000, thrust_n=10.0, isp_s=2500, propulsion='electric')
for nseg in (24, 48):
    t = time.time()
    leg = LT.heliocentric('earth', 'mars', ship, ship.mass(ship.fuel_capacity_kg), (10600, 10900), (300, 1200), nseg=nseg, starts=8)
    c = leg['check']
    print(f"{nseg} segments: tof {leg['tof_days']:.0f} d, fuel {leg['fuel_kg']:.0f} kg (best Sims-Flanagan {leg['best_sims_flanagan_fuel_kg']:.0f}), "
          f"dv eq {leg['dv_equivalent_kms']:.3f} km/s, slack {leg['vinf_dep_ms']:.2f}/{leg['vinf_arr_ms']:.2f} m/s | check: "
          f"{c['position_error_km']:.4f} km, {c['velocity_error_ms']:.6f} m/s, {c['mass_error_kg']:.5f} kg | {time.time() - t:.0f} s")
