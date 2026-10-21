"""Shrink many covers with vktm/wlsm2big under a time limit, then verify (and clean) every improvement with zerocov.

usage: python3 scripts/reshrink_batch.py jobs.json outdir
  jobs.json: a list of {"v","k","t","m","b","path"} (path = an existing full cover with b blocks; plain covers use m = t)
  Writes outdir/results.tsv (one line per job, resumable) and outdir/out/v,k,t,m.txt for every verified improvement.

Covers are grouped by C(v,m): up to 1e6 run 11 at a time on one core each, up to 3e8 three at a time on 4 cores, up to
1.6e9 one at a time on all cores. Each run stops after STALL seconds without a new best cover or MAXSECS in total; a
watchdog kills runs that overrun and keeps whatever cover they saved. If wlsm2big refuses a cover as too expensive
(exit 3), the original cover only gets zero-coverage cleanup by zerocov.

Build first:
  g++ -O3 -march=native -fopenmp -o vktm/wlsm2big vktm/wlsm2big.cpp
  g++ -O3 -march=native -fopenmp -o scripts/zerocov scripts/zerocov.cpp
"""
import json, os, subprocess, sys, glob, re, shutil, time
from math import comb
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVER = os.path.join(ROOT, 'vktm', 'wlsm2big'); CHECKER = os.path.join(ROOT, 'scripts', 'zerocov')
MAXSECS, STALL = 600, 60
SETTINGS = dict(POOL='2000', POOL_REFRESH='500', NORESET='1', BUMP_GROWTH='30', LEVEL_BUDGET='3000', GROWTH='50',
                REVERT_SECS='120', SAVE_SECS='30', STALL=str(STALL))


def main(jobs_file, outdir):
    work, out, res = os.path.join(outdir, 'work'), os.path.join(outdir, 'out'), os.path.join(outdir, 'results.tsv')
    os.makedirs(work, exist_ok=True); os.makedirs(out, exist_ok=True)
    if not os.path.exists(res): open(res, 'w').write('v\tk\tt\tm\tb_old\tb_new\tgain\tstatus\tsecs\n')
    done = {tuple(map(int, l.split('\t')[:4])) for l in open(res).read().splitlines()[1:] if l}
    jobs = json.load(open(jobs_file))
    for j in jobs: j['NM'] = comb(j['v'], j['m'])

    def run(j, threads):
        key = (j['v'], j['k'], j['t'], j['m']); name = '%d,%d,%d,%d' % key
        if key in done: return
        wd = os.path.join(work, name.replace(',', '_')); os.makedirs(wd, exist_ok=True)
        for f in glob.glob(os.path.join(wd, 'mcover_*')): os.remove(f)
        shutil.copy(j['path'], os.path.join(wd, 'start.txt'))
        env = dict(os.environ, OMP_NUM_THREADS=str(threads), **SETTINGS); t0 = time.time(); rc = None
        with open(os.path.join(wd, 'log.txt'), 'w') as lg:
            try:
                rc = subprocess.run([SOLVER, *map(str, key), str(j['b']), str(MAXSECS), '1', 'start.txt', '1000', '5', '1'],
                                    cwd=wd, env=env, stdout=lg, stderr=subprocess.STDOUT, timeout=MAXSECS + 600).returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
        files = glob.glob(os.path.join(wd, 'mcover_*'))
        if rc == 3: files = [os.path.join(wd, 'start.txt')]
        size = lambda f: int(re.findall(r'_(\d+)_\d+\.txt$', f)[0]) if 'mcover_' in f else j['b']
        best = min(files, key=size) if files else None
        status, final = 'no_gain', j['b']
        if best and (size(best) < j['b'] or rc == 3):
            clean = os.path.join(out, name + '.txt')
            r = subprocess.run([CHECKER, *map(str, key), best, clean], capture_output=True, text=True,
                               env=dict(os.environ, OMP_NUM_THREADS=str(threads)))
            m1 = re.search(r'-> (\d+) blocks, uncovered after: 0', r.stdout)
            if r.returncode == 0 and m1 and ' 0 of ' in r.stdout and int(m1.group(1)) < j['b']:
                status, final = 'verified', int(m1.group(1))
            else:
                status = 'VERIFY_FAILED' if r.returncode != 0 else 'no_gain'
                if os.path.exists(clean): os.remove(clean)
        with open(res, 'a') as f:
            f.write('\t'.join(map(str, [*key, j['b'], final, j['b'] - final, status, round(time.time() - t0)])) + '\n')
        print(name, j['b'], '->', final, status, flush=True)

    tiers = (([j for j in jobs if j['NM'] <= 1e6], 11, 1), ([j for j in jobs if 1e6 < j['NM'] <= 3e8], 3, 4),
             ([j for j in jobs if 3e8 < j['NM'] <= 1.6e9], 1, 12))
    for tier, workers, threads in tiers:
        with ThreadPoolExecutor(workers) as ex: list(ex.map(lambda j: run(j, threads), tier))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
