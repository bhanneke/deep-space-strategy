# Probe (not a test): does the continuous-thrust solver converge for a hypothetical high-thrust, high-Isp ship?
# Earth (300 km) to Mars (300 km), 2028-2029 window.
import sys, time
sys.path.insert(0, '/work'); sys.path.insert(0, '/work/flyby')
from engine.ship import Ship
from engine import lowthrust as LT

def mjd(y):
    return (y - 2000) * 365.25 - 0.5

for thrust in (500.0, 2000.0, 5000.0):
    ship = Ship('high-Isp probe', dry_kg=40000, payload_kg=20000, fuel_capacity_kg=40000, thrust_n=thrust, isp_s=10000)
    t = time.time()
    r = LT.transfer('earth', 'mars', ship, (mjd(2028), mjd(2029.5)), (40, 400))
    leg = r.get('leg') or {}
    print(f'thrust {thrust:.0f} N, a0 {1000 * thrust / ship.mass(ship.fuel_capacity_kg):.1f} mm/s2, low_thrust {ship.low_thrust}, '
          f'feasible {r["feasible"]}, days {r.get("arrive_mjd2000", 0) - r.get("depart_mjd2000", 0):.0f}, fuel {r.get("fuel_kg", 0):.0f} kg, '
          f'check {leg.get("check", {}).get("position_error_km")} km, {time.time() - t:.0f} s')
    for p in r['phases']:
        print('   ', p['phase'], round(p['dv_kms'], 2), 'km/s', round(p['days'], 1), 'days', round(p['fuel_kg']), 'kg')
