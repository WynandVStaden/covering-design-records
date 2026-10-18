#!/bin/bash
# Run the shrink search on many entries in parallel.
#
# usage: bash scripts/run_batch.sh <seconds_per_entry> <jobs_file> [parallel_jobs]
#
# Each line of the jobs file is:  v k t b start_file
#   b          = number of blocks in start_file (the current record cover)
#   start_file = a valid covering to start from (one block per line, points 1..v)
# For each entry the solver shrinks the cover for the given time and saves every full cover it reaches
# into results/C_v_k_t/. The log line reports the smallest full cover found.
SECS=${1:?seconds per entry}; JOBS=${2:?jobs file}; P=${3:-$(nproc)}
BIN=$(dirname "$0")/wls128
mkdir -p results
tr -d '\r' < "$JOBS" | xargs -P "$P" -L 1 bash -c '
  out=$("'"$BIN"'" $0 $1 $2 $3 '"$SECS"' 1 $4 7 0.05 1 2>/dev/null)
  best=$(echo "$out" | grep "SHRINK result" | grep -o "[0-9]* blocks" | grep -o "[0-9]*")
  echo "C($0,$1,$2) start $3 -> best full ${best:-none}"
  mkdir -p results/C_${0}_${1}_${2}
  mv -f shrink_cover_${0}_${1}_${2}_*_1.txt results/C_${0}_${1}_${2}/ 2>/dev/null'
