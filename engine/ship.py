"""Ships and their mass accounting.

A ship carries dry mass (structure, engines, tanks), payload and fuel. Every burn obeys the rocket equation for the
ship's current total mass: m_after = m_before * exp(-dv / ve), ve = Isp * g0. The engine keeps fuel, payload and dry
mass separately, so a journey can unload payload, refuel at a stop and carry the remaining fuel to the next leg.

Thrust enters in two places. For impulsive manoeuvres it gives the burn duration, which the engine compares with the
parking-orbit period: a burn that lasts longer than about a tenth of an orbit loses efficiency (gravity loss), which
the impulsive model does not capture, so the result carries a warning. For low-thrust ships the heliocentric leg uses
the continuous-thrust model in lowthrust.py instead.
"""
from dataclasses import dataclass, asdict, field
import math

G0 = 9.80665                          # m/s^2, standard gravity (defines Isp)
FINITE_BURN_LIMIT = 0.1               # burn duration as share of the parking-orbit period above which a warning shows

@dataclass
class Ship:
    name: str
    dry_kg: float                     # structure, engines, tanks, without payload or fuel
    payload_kg: float
    fuel_capacity_kg: float
    thrust_n: float
    isp_s: float
    propulsion: str = 'chemical'      # 'chemical' or 'nuclear thermal' (impulsive) or 'electric' (low thrust)

    @property
    def ve_kms(self):
        return self.isp_s * G0 / 1000.0

    @property
    def mdot_kgs(self):
        return self.thrust_n / (self.isp_s * G0)

    @property
    def low_thrust(self):
        return self.propulsion == 'electric'

    def mass(self, fuel_kg, payload_kg=None):
        return self.dry_kg + (self.payload_kg if payload_kg is None else payload_kg) + fuel_kg

    def dv_capacity_kms(self, fuel_kg=None, payload_kg=None):
        """largest velocity change with the given fuel (full tanks by default)"""
        f = self.fuel_capacity_kg if fuel_kg is None else fuel_kg
        m0 = self.mass(f, payload_kg)
        return self.ve_kms * math.log(m0 / (m0 - f))

    def to_dict(self):
        return asdict(self)

@dataclass
class BurnResult:
    dv_kms: float
    fuel_before_kg: float
    fuel_used_kg: float
    fuel_after_kg: float
    mass_before_kg: float
    mass_after_kg: float
    duration_s: float
    feasible: bool
    warnings: list = field(default_factory=list)

def burn(ship, fuel_kg, dv_kms, payload_kg=None, parking_period_s=None):
    """one impulsive manoeuvre with the ship's current fuel; infeasible when the fuel does not cover it"""
    m0 = ship.mass(fuel_kg, payload_kg)
    m1 = m0 * math.exp(-dv_kms / ship.ve_kms)
    used = m0 - m1
    ok = used <= fuel_kg + 1e-9
    duration = used / ship.mdot_kgs if ship.mdot_kgs > 0 else float('inf')
    warn = []
    if parking_period_s and duration > FINITE_BURN_LIMIT * parking_period_s:
        warn.append(f'burn lasts {duration / 60:.0f} min, {100 * duration / parking_period_s:.0f}% of the parking-orbit period: '
                    'finite-burn losses are not modelled')
    return BurnResult(dv_kms, fuel_kg, used, fuel_kg - used if ok else 0.0, m0, m1 if ok else m0 - fuel_kg, duration, ok, warn)

# Illustrative presets (rounded generic values, not specific vehicles). Edit them or define your own ship.
PRESETS = {    # full-tank capacity with payload: methalox 7.2, hydrolox 7.9, nuclear thermal 7.5 km/s
    'methalox tug':        Ship('methalox tug', dry_kg=12000, payload_kg=20000, fuel_capacity_kg=200000, thrust_n=2.0e6, isp_s=370),
    'hydrolox tug':        Ship('hydrolox tug', dry_kg=10000, payload_kg=10000, fuel_capacity_kg=100000, thrust_n=4.0e5, isp_s=450),
    'nuclear thermal':     Ship('nuclear thermal', dry_kg=25000, payload_kg=20000, fuel_capacity_kg=60000, thrust_n=3.0e5, isp_s=900,
                                propulsion='nuclear thermal'),
    'solar-electric cargo': Ship('solar-electric cargo', dry_kg=8000, payload_kg=20000, fuel_capacity_kg=12000, thrust_n=2.0, isp_s=2500,
                                 propulsion='electric'),
}
