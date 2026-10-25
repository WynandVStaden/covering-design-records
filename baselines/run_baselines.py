"""Same-budget comparison of the weighted search with tabu and annealing baselines, from the same start covers.

usage (inside WSL): python3 baselines/run_baselines.py bindir startdir outdir [seconds] [seeds]
  bindir: compiled wls2, wlsm2, wlsm, sa (see README.md); startdir: previous table covers named C_v_k_t.txt /
  C_v_k_t_m.txt; outdir: one line per run in results.tsv (resumable).
Methods for covering designs: weighted (wls2), tabu (wls2 with NOWEIGHT=1), annealing (sa).
Methods for lottery coverings: weighted (wlsm2), tabu (wlsm2 with NOWEIGHT=1), direct (wlsm, m-subset scoring).
Every run is decremental from the same start cover and is single-threaded; runs are shuffled so that machine load
affects all methods alike. The reported size is the smallest full cover the run saved.
"""
import os, sys, random, re, subprocess, time
from concurrent.futures import ThreadPoolExecutor

PLAIN = [(57, 12, 3), (82, 14, 3), (97, 15, 3), (31, 11, 4), (75, 20, 4), (80, 19, 4), (32, 15, 5), (39, 18, 5),
         (33, 18, 6), (39, 20, 6)]
LOTTERY = [(61, 8, 3, 4), (66, 8, 3, 4), (78, 9, 3, 4), (52, 15, 4, 5), (56, 13, 4, 5), (58, 15, 4, 5)]


def nblocks(p): return sum(1 for l in open(p) if l.strip())


def main(bindir, startdir, outdir, secs=60, seeds=5):
    os.makedirs(outdir, exist_ok=True); res = os.path.join(outdir, 'results.tsv')
    done = set()
    if os.path.exists(res):
        for l in open(res).read().splitlines()[1:]: done.add(tuple(l.split('\t')[:3]))
    else:
        open(res, 'w').write('entry\tmethod\tseed\tstart\tbest\tsecs\n')
    jobs = []
    for e in PLAIN:
        for m in ('weighted', 'tabu', 'annealing'):
            for s in range(1, seeds + 1): jobs.append((e, m, s))
    for e in LOTTERY:
        for m in ('weighted', 'tabu', 'direct'):
            for s in range(1, seeds + 1): jobs.append((e, m, s))
    random.Random(0).shuffle(jobs)

    def run(job):
        e, method, seed = job; name = ','.join(map(str, e))
        if (name, method, str(seed)) in done: return
        start = os.path.join(startdir, 'C_' + '_'.join(map(str, e)) + '.txt'); b = nblocks(start)
        wd = os.path.join(outdir, 'work', '%s_%s_%d' % (name.replace(',', '_'), method, seed)); os.makedirs(wd, exist_ok=True)
        env = dict(os.environ, OMP_NUM_THREADS='1')
        if method == 'tabu': env['NOWEIGHT'] = '1'
        if len(e) == 3:
            v, k, t = e
            if method == 'annealing': cmd = [os.path.join(bindir, 'sa'), v, k, t, b, secs, seed, start, 1]
            else: cmd = [os.path.join(bindir, 'wls2'), v, k, t, b, secs, seed, start, 7, 0.05, 1]
        else:
            v, k, t, m = e
            prog = 'wlsm' if method == 'direct' else 'wlsm2'
            cmd = [os.path.join(bindir, prog), v, k, t, m, b, secs, seed, start, 1000, 5, 1]
        t0 = time.time()
        out = subprocess.run(list(map(str, cmd)), cwd=wd, env=env, capture_output=True, text=True,
                             timeout=secs + 300).stdout
        sizes = [int(x) for x in re.findall(r'FULL (\d+) blocks', out)] + \
                [int(x) for x in re.findall(r'full cover with (\d+) blocks', out)]
        best = min(sizes) if sizes else b
        open(os.path.join(wd, 'log.txt'), 'w').write(out)
        with open(res, 'a') as f: f.write('%s\t%s\t%d\t%d\t%d\t%.0f\n' % (name, method, seed, b, best, time.time() - t0))
        print(name, method, seed, b, '->', best, flush=True)

    with ThreadPoolExecutor(11) as ex: list(ex.map(run, jobs))
    print('ALL DONE', flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], a[2], int(a[3]) if len(a) > 3 else 60, int(a[4]) if len(a) > 4 else 5)
