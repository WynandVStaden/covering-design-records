#!/bin/bash
# Runs (or resumes) the whole random-sample experiment: 300 plain entries, then 300 lottery entries.
# Safe to stop at any time (shutdown, or rs/stop.bat); finished runs are kept in results.tsv and skipped on restart,
# an interrupted or garbled run is repeated from the start. Each phase is repeated until it reports ALL DONE.
RS=/mnt/c/Users/wvs00/claude_chat/coverings/rs
REPO=/mnt/c/Users/wvs00/covering-design-records
BIN=/home/ubuntu/bl
cd $REPO
for kind in plain lottery; do
  for attempt in 1 2 3; do
    python3 baselines/random_sample/run.py $BIN $RS/sample_${kind}300.tsv $RS/out_$kind 120 3 > $RS/run_$kind.log 2>&1
    grep -q "ALL DONE" $RS/run_$kind.log && break
  done
done
echo "EXPERIMENT DONE $(date)" >> $RS/run_all.log
