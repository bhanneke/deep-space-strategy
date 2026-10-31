# solmaps routing engine

The engine finds the itineraries a given ship can fly between two sites within a departure window, with arrival
dates, fuel use and delivered payload. It combines precomputed transfer alternatives with ship
calculations at query time and refines trajectories on demand.

## Architecture

| Layer | Module | When it runs | What it holds |
| --- | --- | --- | --- |
| End states | `engine/conventions.py` | always | how a ship leaves and enters each node: parking orbit, rendezvous or Earth–Moon L1 |
| Transfer alternatives | `engine/alternatives.py` | `make data`, 30 s | per route and departure date up to six Lambert transfers that trade flight time against cost |
| Ships | `engine/ship.py` | per query | dry mass, payload, fuel capacity, thrust, specific impulse; rocket equation per burn |
| Journeys | `engine/journey.py` | per query, under a second | chronological multi-leg journeys with waiting, deadline and refuelling |
| Low thrust | `engine/lowthrust.py` | on request, 10–60 s | continuous-thrust transfers with changing mass |
| Refinement cache | `engine/refine.py` | on request | exact impulsive transfers at a given date; low-thrust transfers; stored under `data/cache/` |
| Service and page | `engine/api.py`, `engine/ui/index.html` | `make serve` | compare ships on one trip |

## Conventions

Every transfer has a heliocentric part and planet-side parts, and the engine keeps them apart. The excess speeds (v∞) at both ends describe the heliocentric leg. The planet-side burn turns a parking orbit into the departure
hyperbola, or the arrival hyperbola into a parking orbit:

Δv = √(v∞² + 2μ/r) − √(μ/r).

| Node | End state |
| --- | --- |
| Earth, Venus, Mars | circular orbit 300 km above the surface |
| Mercury | 200 km |
| Ceres | 50 km |
| Jupiter, Saturn | circular orbit at the distance of Callisto and Titan |
| Vesta, Psyche, Eros, Bennu, Ryugu | rendezvous: the burn equals v∞ |
| Earth–Moon L1 | the cheaper of a burn at L1 and a drop to a 300 km perigee with a perigee burn; 3.77 km/s to or from low Earth orbit |

The network, the flyby search, the station precompute and the engine import these values from one module. A test
compares them with the JavaScript engine's values (`network/nodes.js`).

## Transfer alternatives

For each ordered pair of the 12 sites and each departure date (every 10 days, November 2026 to 2100), the precompute
solves Lambert's problem for 40 flight times: the network's 30 values (0.3–1.8 Hohmann times) plus 10 longer ones up to
2.5 Hohmann times. It keeps the transfers that are best for some trade-off of arrival date against cost (a Pareto front), thinned to six.
The cheapest alternative equals the network value whenever it flies within the network's flight-time grid. In a third
of the sampled route-dates, a flight longer than the network's limit costs less, by a median of 4.6 km/s, mostly on
dates between launch windows.

## Ships and mass accounting

A ship's mass is dry mass plus payload plus fuel. Every burn obeys the rocket equation for the current mass,
m_after = m_before · exp(−Δv / (Isp · g0)), and a burn is infeasible when the fuel on board does not cover it. A leg
has two burns: the departure burn on the ship's mass at departure, then the arrival burn on the lighter ship. Thrust
gives the burn duration; when a burn lasts longer than a tenth of the parking-orbit period, the result carries a
warning, because the impulsive model ignores the losses of long burns.

The four presets use rounded generic values for illustration:

| Preset | Dry t | Payload t | Fuel t | Thrust | Isp s | Full-tank capacity |
| --- | --- | --- | --- | --- | --- | --- |
| methalox tug | 12 | 20 | 200 | 2 MN | 370 | 7.2 km/s |
| hydrolox tug | 10 | 10 | 100 | 400 kN | 450 | 7.9 km/s |
| nuclear thermal | 25 | 20 | 60 | 300 kN | 900 | 7.5 km/s |
| solar-electric cargo | 8 | 20 | 30 | 10 N | 2,500 | 16.7 km/s |

## Journeys

The search runs on a time-expanded graph. A partial journey (label) holds its node, its arrival date, the fuel on
board and the fuel loaded so far. From a label the ship may leave on any network date within a waiting limit (default
1,100 days) on any transfer alternative. A leg leaves no earlier than the ship's arrival. The ship refills to capacity
at nodes that supply propellant and that the query allows: Earth, Earth–Moon L1, Mars, Ceres, Bennu, Ryugu. The ship
stops only there. The deadline applies to each leg's actual arrival time. The search drops a label when another label at
the same node arrived no later, holds at least as much fuel and has loaded no more. The result lists every journey on
the Pareto front of arrival date, fuel loaded and fuel left. For each journey, `max_payload` finds the largest payload
that fits the same legs and refuels.

The departure window covers the network dates inside it. The engine re-derives each returned itinerary burn by burn with the
ship model and checks that fuel passes from leg to leg.

## Low thrust

A low-thrust transfer has three phases, reported separately:

1. **Departure spiral.** From a parking orbit, the ship spirals out to escape. For a slow tangential spiral, the
   velocity change equals the circular speed of the parking orbit (Edelbaum): 7.73 km/s from 300 km above Earth. Time
   and fuel follow from thrust and mass flow. At small bodies the ship starts at rest relative to the body.
2. **Heliocentric leg.** pykep's Sims–Flanagan transcription (`sf_pl2pl`) finds the route from several random starts.
   Its throttle profile seeds pykep's zero-order-hold model (`zoh_pl2pl`): thrust constant in magnitude and direction
   within each of 24 segments, integrated with the heyoka Taylor integrator together with the falling mass,
   dm/dt = −|T|/ve. IPOPT maximises the arrival mass subject to position, velocity and mass continuity. SLSQP stalls on
   this problem.
3. **Capture spiral.** The reverse spiral into the destination's parking orbit; none at small bodies.

The zero-order-hold problem needs a non-zero range for the excess speeds at both ends. The engine allows up to
3 m/s, which the solution uses as free velocity: about 0.1% of a typical transfer, and reported with every result.

### Checks

The independent check integrates the solution's thrust profile with scipy's DOP853 integrator: two-body gravity plus
thrust, mass falling by |T|/ve. This code is separate from pykep's propagators. The check compares the end of the
integration with the arrival planet and with the optimiser's final mass. Earth to Mars, departures January to
November 2029, a 40 t ship with 10 N at Isp 2,500 s:

| Segments | Flight | Fuel | Best Sims–Flanagan fuel | Position error | Velocity error | Mass error |
| --- | --- | --- | --- | --- | --- | --- |
| 24 | 950 d | 8,335 kg | 8,265 kg | 0.41 km | 0.00004 m/s | 0.00002 kg |
| 48 | 979 d | 8,248 kg | 8,313 kg | 0.26 km | 0.00003 m/s | < 0.00001 kg |

The two transcriptions agree within 1% on fuel. With Sims–Flanagan alone (impulses at segment midpoints), the same
check showed errors of 8.6 million km at 24 segments, 2.4 million km at 48 and 0.28 million km at 96. That pattern of
discretisation error is why the engine refines in the continuous model.

## Earlier review findings, checked against the code (2026-10-31)

| Finding | Verdict | Change |
| --- | --- | --- |
| Inconsistent Ceres arrival conventions | Confirmed: the second flyby run priced Ceres as a rendezvous, the network as capture into a 50 km orbit | One conventions module for every component; the flyby reruns (routes3) use it |
| Independent inbound and outbound hub scoring | Confirmed in the growth map and in the station evaluation | Onward legs leave after the arrival. Earth–Moon L1 leverage 52–84% with slingshots (was 57–90%); station fuel saved 70% for L1 (was 82%) |
| Low thrust approximated by a speed multiplier | Not found; the code had no low-thrust model, only exhaust speeds applied to impulsive Δv | Continuous-thrust model above |
| Only the minimum Δv per date | Confirmed | Transfer alternatives |

## Build and test

```sh
make image      # pinned container: base image by digest, conda and pip packages by version
make data       # network to 2100 and transfer alternatives (under a minute)
make test       # 39 tests, about 40 s
make serve      # comparison page at http://localhost:8077/
make pairs site # station edges and the public pages into docs/
```

The tests cover:
- the conventions in every component, including the JavaScript engine
- the rocket equation and fuel limits
- the alternatives against the network
- journey chronology, deadlines, refuelling, departure windows and maximum payload
- the low-thrust check, mass consistency, throttle bound, spiral formula and the agreement of the two transcriptions
- the cache, the service and its input validation

## Limits

- Impulsive legs are single-revolution Lambert arcs on a 10-day grid; the refinement cache refines one transfer at an
  exact date on request.
- The journey search treats each heliocentric leg independently of the next. Slingshot routes are not yet edges of the
  journey graph, and planetary flybys appear only in the flyby search.
- Low-thrust transfers cover one leg between two sites with spirals at both ends. Low-thrust ships are not part of the
  multi-leg journey search. The spiral uses Edelbaum's formula, without eclipses or the third-body effect of the Moon.
- No finite-burn losses for chemical ships, only the warning; no aerocapture in the engine yet (the cargo planner in
  the 3D atlas assumes it up to 6 km/s, an unsourced assumption).
- The result list of a journey query can hold near-identical journeys 10 days apart; a thinning step would make it easier to read.
- Ship presets are illustrative.
