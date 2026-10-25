#!/bin/bash
# Runs (or resumes) the whole random-sample experiment: 300 plain entries, then 300 lottery entries.
# Safe to stop at any time (shutdown, or rs/stop.bat); finished runs are kept in results.tsv and skipped on restart,
# an interrupted run is repeated from the start.
RS=/mnt/c/Users/wvs00/claude_chat/coverings/rs
REPO=/mnt/c/Users/wvs00/covering-design-records
BIN=/home/ubuntu/bl
cd $REPO
# wait for a run started earlier (e.g. the first 100 plain entries) instead of competing with it
while pgrep -f 'random_sample/run.py' >/dev/null; do sleep 20; done
python3 baselines/random_sample/run.py $BIN $RS/sample_plain300.tsv $RS/out_plain 120 3 >> $RS/run_plain.log 2>&1
while [ ! -f $RS/sample_lottery300.tsv ]; do sleep 60; done
python3 baselines/random_sample/run.py $BIN $RS/sample_lottery300.tsv $RS/out_lottery 120 3 >> $RS/run_lottery.log 2>&1
echo "EXPERIMENT DONE $(date)" >> $RS/run_all.log
