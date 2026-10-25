"""Run the pre-registered random-sample comparison (PROTOCOL.md).

usage (inside WSL): python3 baselines/random_sample/run.py bindir sample.tsv outdir [seconds=120] [seeds=3]
  sample.tsv: stratum, v, k, t, m, table_value, start (path of the start cover)
Writes outdir/results.tsv (entry, stratum, method, seed, table_value, start_size, best), resumable.
"""
import csv, os, random, re, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor


def main(bindir, sample, outdir, secs=120, seeds=3):
    os.makedirs(outdir, exist_ok=True); res = os.path.join(outdir, 'results.tsv')
    done = set()
    if os.path.exists(res):
        for l in open(res).read().splitlines()[1:]:
            x = l.split('\t')
            if len(x) == 7: done.add((x[0], x[2], x[3]))   # a garbled line means that run is repeated
    else:
        open(res, 'w').write('entry\tstratum\tmethod\tseed\ttable_value\tstart_size\tbest\n')
    rows = list(csv.DictReader(open(sample), delimiter='\t'))
    jobs = []
    for r in rows:
        plain = r['t'] == r['m']
        for m in (('weighted', 'tabu', 'annealing') if plain else ('weighted', 'tabu', 'direct')):
            for s in range(1, seeds + 1): jobs.append((r, m, s))
    random.Random(1).shuffle(jobs)
    lock = threading.Lock()   # runs finish concurrently; without it two result lines can interleave

    def run(job):
        r, method, seed = job
        v, k, t, m = (int(r[c]) for c in ('v', 'k', 't', 'm')); name = '%d,%d,%d,%d' % (v, k, t, m)
        if (name, method, str(seed)) in done: return
        start = r['start']; b = sum(1 for l in open(start) if l.strip())
        wd = os.path.join(outdir, 'work', '%s_%s_%d' % (name.replace(',', '_'), method, seed)); os.makedirs(wd, exist_ok=True)
        env = dict(os.environ, OMP_NUM_THREADS='1')
        if method == 'tabu': env['NOWEIGHT'] = '1'
        if t == m:
            cmd = [os.path.join(bindir, 'sa'), v, k, t, b, secs, seed, start, 1] if method == 'annealing' else \
                  [os.path.join(bindir, 'wls2'), v, k, t, b, secs, seed, start, 7, 0.05, 1]
        else:
            cmd = [os.path.join(bindir, 'wlsm' if method == 'direct' else 'wlsm2'), v, k, t, m, b, secs, seed, start,
                   1000, 5, 1]
        try:
            out = subprocess.run(list(map(str, cmd)), cwd=wd, env=env, capture_output=True, text=True,
                                 timeout=secs + 300).stdout
        except subprocess.TimeoutExpired as e:
            out = (e.stdout or b'').decode() if isinstance(e.stdout, bytes) else (e.stdout or '')
        sizes = [int(x) for x in re.findall(r'FULL (\d+) blocks', out)] + \
                [int(x) for x in re.findall(r'full cover with (\d+) blocks', out)]
        best = min(sizes) if sizes else b
        open(os.path.join(wd, 'log.txt'), 'w').write(out)
        with lock, open(res, 'a') as f:
            f.write('%s\t%s\t%s\t%d\t%s\t%d\t%d\n' % (name, r['stratum'], method, seed, r['table_value'], b, best))

    with ThreadPoolExecutor(11) as ex: list(ex.map(run, jobs))
    print('ALL DONE', flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], a[2], int(a[3]) if len(a) > 3 else 120, int(a[4]) if len(a) > 4 else 3)
