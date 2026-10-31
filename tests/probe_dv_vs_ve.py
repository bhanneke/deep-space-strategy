# Probe (not a test): does the heliocentric continuous-thrust leg's equivalent delta-v depend on the exhaust speed at a
# fixed initial acceleration? If not, a grid over acceleration alone serves every engine (fuel from the rocket equation).
# Earth to Mars, departures 2029, flight 300-1200 days, start mass 40 t, no planet-side phases.
import sys, time
sys.path.insert(0, '/work'); sys.path.insert(0, '/work/flyby')
from engine.ship import Ship
from engine import lowthrust as LT
M0 = 40000.0
for a0 in (0.25e-3, 2e-3, 20e-3):                          # m/s^2
    row = []
    for isp in (2000, 4000, 10000):
        ship = Ship('probe', dry_kg=10000, payload_kg=10000, fuel_capacity_kg=20000, thrust_n=a0 * M0, isp_s=isp)
        t = time.time()
        leg = LT.heliocentric('earth', 'mars', ship, M0, (10600, 10900), (300, 1200), nseg=24, starts=8)
        row.append(f"Isp {isp}: {leg['dv_equivalent_kms']:.3f} km/s, {leg['tof_days']:.0f} d" if leg else f'Isp {isp}: none')
    print(f'a0 {1000 * a0:g} mm/s2 | ' + ' | '.join(row), flush=True)
