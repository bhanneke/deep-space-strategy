# solmaps: the engine of Deep Space Strategy

Deep Space Strategy studies strategy for a settled solar system in five modules: Geography, Resources, Topology, Infrastructure and Strategy. solmaps is its open map and route engine: a queryable timetable of routes through the solar system, and the strategy analysis built on it.

solmaps computes transfer costs between solar-system sites for every departure date, ranks propellant hubs by how much Earth-launched mass they save, and searches for depot orbits around the Sun. The technical documentation (concept, formulas, references) is at https://bhanneke.github.io/solmaps/technical.html.

## Contents

| Folder | What it holds |
| --- | --- |
| `atlas/` | Resource atlas: verified dataset (`atlas_dataset.json`, every entry with evidence level and source) and the interactive map (`build.py` renders `solar-system-atlas.html`) |
| `network/` | Temporal transfer network: Lambert solver and parking-orbit costs (`engine.js`), route precomputation (`compute.js`), depot-orbit search (`stations.js`), validation (`validate.js`), growth-map page (`build.py`) |
| `env/` | Container environment with ESA pykep 3.0.0 and pygmo from conda-forge, for flyby sequences |
| `flyby/` | DE440 bodies, the Python direct network to 2100, the flyby search (Earth and Mars origins), cycler metric, mission validation |
| `atlas3d/` | The 3D atlas with transfer arcs and launchable ships |
| `surrogate/` | Route dataset generator and the first learned surrogate (PyTorch, separate container) |
| `site/` | Landing page, documentation and `build_site.py`, which writes the public site into `docs/` for GitHub Pages |

The public site and its documentation ("How solmaps works") live at https://bhanneke.github.io/solmaps/. The site is an experimental research prototype.

## Rebuild

```sh
cd atlas && python3 build.py
cd ../network && node validate.js && node compute.js && node stations.js && python3 build.py
docker build -t solmaps-env env/
```

## Status

- Direct transfers validated: departure from a 300 km Earth orbit to Mars 3.61 km/s against 3.59 km/s derived for the Hohmann transfer; pykep Lambert agrees (vinf 3.03 km/s, 2026-10-30, 308 days).
- Planet positions use JPL Table 1 (valid 1800 to 2050). Release 1 covers departures to 2050.
- Not yet included: gravity assists, low thrust, the learned surrogate.

## Licence

Code: MIT (`LICENSE`). Data and documentation: CC BY 4.0 (`LICENSE-DATA.md`). pykep, used in the container, is MPL-2.0.
