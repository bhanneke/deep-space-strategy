#!/bin/zsh -l
# routes3 flyby searches on the Mac Mini, one Docker run at a time (container capped at 2.5 GB, 6 pool workers).
# Guard: the Mini swaps to its root volume and panicked once when that volume filled (2026-10-19), so every job waits
# for >= 40 GB free there, and a running job stops if free space falls below that; routes3 resumes it later.
# Start: nohup zsh ~/solmaps/flyby/run_queue_mini.sh > /dev/null 2>&1 &     Log: ~/solmaps/queue.log
cd ~/solmaps
JOBS=(ceres:earth,mars vesta:earth,mars psyche:earth,mars jupiter:earth,mars
      earth:jupiter,saturn,mercury earth:ceres,vesta,psyche mars:jupiter,saturn mars:ceres,vesta,psyche)
log() { echo "$(date '+%F %T') $*" >> queue.log; }
room() { [ "$(df -g / | awk 'NR==2{print $4}')" -ge 40 ]; }
for job in $JOBS; do
  origin=${job%%:*}; targets=${job#*:}; tag=${origin}_${targets//,/_}; name=solmaps-$tag
  for attempt in 1 2 3; do
    while ! room; do log "waiting: root volume below 40 GB free"; sleep 1800; done
    log "start $job (attempt $attempt)"
    docker run --rm --name $name --memory=2500m -e SOLMAPS_POOL=6 -v ~/solmaps:/work -w /work/flyby solmaps-env \
      python /work/flyby/routes3.py $origin $targets 2027 2097 8 1000 60 6 3 >> flyby/run3_$tag.log 2>&1 &
    pid=$!
    while kill -0 $pid 2>/dev/null; do
      if ! room; then log "stopping $job: root volume below 40 GB free"; docker stop $name > /dev/null; fi
      sleep 600
    done
    wait $pid; code=$?
    log "end $job exit $code"
    [ $code -eq 0 ] && break
  done
done
log "queue finished"
