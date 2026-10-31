#!/bin/zsh -l
# Second queue on the Mac Mini, next to run_queue_mini.sh: everything else not yet computed.
#   1. flyby search for the pairs no run has covered (Mars to Earth, Venus, Mercury; Earth to Mars, Venus and the
#      near-Earth asteroids; Saturn to Earth and Mars)
#   2. the full Earth-Mars cycler scan, 2027-2097
#   3. the low-thrust atlas: the electric preset on every site pair and two-year window (cached)
#   4. once the first queue has finished: sampled trajectories for every flyby result
# Same disk guard as the first queue; containers capped at 2 GB with 4 pool workers.
# Start: nohup zsh ~/solmaps/flyby/run_queue_mini2.sh > /dev/null 2>&1 &     Log: ~/solmaps/queue2.log
cd ~/solmaps
IMAGE=solmaps-env:2
log() { echo "$(date '+%F %T') $*" >> queue2.log; }
room() { [ "$(df -g / | awk 'NR==2{print $4}')" -ge 40 ]; }
run() {   # run NAME CMD...: one container, guarded, up to three attempts
  name=$1; shift
  for attempt in 1 2 3; do
    while ! room; do log "waiting: root volume below 40 GB free"; sleep 1800; done
    log "start $name (attempt $attempt)"
    docker run --rm --name solmaps2-$name --memory=2000m -e SOLMAPS_POOL=4 -v ~/solmaps:/work -w /work $IMAGE "$@" >> flyby/q2_$name.log 2>&1 &
    pid=$!
    while kill -0 $pid 2>/dev/null; do
      if ! room; then log "stopping $name: root volume below 40 GB free"; docker stop solmaps2-$name > /dev/null; fi
      sleep 600
    done
    wait $pid; code=$?
    log "end $name exit $code"
    [ $code -eq 0 ] && return 0
  done
}
run flyby_mars    python flyby/routes3.py mars earth,venus,mercury 2027 2097 8 1000 60 6 3
run flyby_earth   python flyby/routes3.py earth mars,venus,433,101955,162173 2027 2097 8 1000 60 6 3
run flyby_saturn  python flyby/routes3.py saturn earth,mars 2027 2097 8 1000 60 6 3
run cyclers       python flyby/cyclers.py 2027 2097 8 2000 30 6
run lowthrust     python -m engine.lt_atlas
while ! grep -q "queue finished" queue.log 2>/dev/null; do log "waiting for the first queue before sampling trajectories"; sleep 3600; done
run trajectories  python flyby/trajectories.py
log "queue 2 finished"
