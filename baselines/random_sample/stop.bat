@echo off
rem Stop the random-sample experiment. Finished runs are kept; run resume.bat to continue later.
wsl -e bash -c "pkill -f 'r[u]n_all.sh'; pkill -f 'random_sampl[e]/run.py'; pkill -x wls2; pkill -x wlsm2; pkill -x wlsm; pkill -x sa; echo stopped"
pause
