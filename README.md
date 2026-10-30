# solmaps

A queryable timetable of routes through the solar system, and the strategy analysis built on it.

solmaps computes transfer costs between solar-system sites for every departure date, ranks propellant hubs by how much Earth-launched mass they save, and searches for depot orbits around the Sun. The project plan lives in the Claude Doc "Solar System Railway: Project Plan".

## Contents

| Folder | What it holds |
| --- | --- |
| `atlas/` | Resource atlas: verified dataset (`atlas_dataset.json`, every entry with evidence level and source) and the interactive map (`build.py` renders `solar-system-atlas.html`) |
| `network/` | Temporal transfer network: Lambert solver and parking-orbit costs (`engine.js`), route precomputation (`compute.js`), depot-orbit search (`stations.js`), validation (`validate.js`), growth-map page (`build.py`) |
| `env/` | Container environment with ESA pykep 3.0.0 from conda-forge, for flyby sequences |

## Rebuild

```sh
cd atlas && python3 build.py
cd ../network && node validate.js && node compute.js && node stations.js && python3 build.py
docker build -t solmaps-env env/
```

## Status

- Direct transfers validated: departure from low Earth orbit to Mars 3.61 km/s against a sourced 3.6 km/s; pykep Lambert agrees (vinf 3.03 km/s, 2026-10-30, 308 days).
- Planet positions use JPL Table 1 (valid 1800 to 2050). Release 1 covers departures to 2050.
- Not yet included: gravity assists, low thrust, the learned surrogate.

## Licence

Not yet chosen. pykep is MPL-2.0, which permits a permissive licence for solmaps' own code.
