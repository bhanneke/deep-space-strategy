#!/bin/sh
# Return-route flyby searches (belt and Jupiter system back to Earth and Mars), 2027-2098.
# Usage: sh flyby/run_returns.sh "ceres vesta"   (one docker run per origin, in sequence)
cd "$(dirname "$0")/.."
for origin in $1; do
  docker run --rm -v "$(pwd)":/work -w /work/flyby solmaps-env python /work/flyby/routes.py "$origin" earth,mars 2027 2097 8 2000 30 6 > "flyby/run3_$origin.log" 2>&1
done
