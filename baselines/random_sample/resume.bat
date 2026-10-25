@echo off
rem Start or resume the random-sample experiment in the background (finished runs are skipped).
wsl -e bash -c "if pgrep -f 'r[u]n_all.sh' >/dev/null; then echo already running; else nohup bash /mnt/c/Users/wvs00/claude_chat/coverings/rs/run_all.sh >/dev/null 2>&1 & sleep 1; echo started; fi"
pause
