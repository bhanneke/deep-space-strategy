#!/bin/sh
# Queue of routes3 flyby searches, 2027-2098, full budget (8 islands, 2,000 generations, 30 individuals, 6 rounds,
# 3 starts; 6 for two-flyby sequences). Usage: sh flyby/run_routes3.sh "ceres:earth,mars psyche:earth,mars"
cd "$(dirname "$0")/.."
for job in $1; do
  origin=${job%%:*}; targets=${job#*:}
  docker run --rm -v "$(pwd)":/work -w /work/flyby solmaps-env python /work/flyby/routes3.py "$origin" "$targets" 2027 2097 8 1000 60 6 3 > "flyby/run3_${origin}_$(echo "$targets" | tr , _).log" 2>&1
done
