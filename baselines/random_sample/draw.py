"""Draw the random sample of table entries for the method comparison (see PROTOCOL.md). Deterministic: seed 20261025.

usage: python baselines/random_sample/draw.py vkt_table vktm_table outdir
  vkt_table / vktm_table: the Covering Repository's table archives (lines 'v k t m size'), downloaded 2026-10-25.
Writes outdir/order_plain.tsv and outdir/order_lottery.tsv: every eligible entry of each stratum in a random order.
The sample is the first 25 entries of each stratum's order for which a start cover can be obtained (PROTOCOL.md).
"""
import csv, hashlib, os, random, sys
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
SEED = 20261025


def table(p):
    r = {}
    for l in open(p):
        x = l.split()
        if len(x) >= 5 and x[4] != '-':
            try: r[tuple(map(int, x[:4]))] = int(x[4])
            except ValueError: pass
    return r


def ours(live):
    """Entries whose live value is one of our coverings (excluded: attacking them tests nothing)."""
    held = set()
    for f, cols in (('search_plain.tsv', 'vktt'), ('search_lottery.tsv', 'vktm'), ('derived.tsv', 'vktm')):
        for r in csv.DictReader(open(os.path.join(REPO, 'data', f)), delimiter='\t'):
            key = tuple(int(r[c]) for c in ('v', 'k', 't', 't' if cols == 'vktt' else 'm'))
            if live.get(key) == int(r['ours']): held.add(key)
    return held


def main(vkt, vktm, outdir):
    live = table(vkt); live.update(table(vktm)); held = ours(live)
    os.makedirs(outdir, exist_ok=True)
    plain_strata = {'t=3': lambda k: k[2] == 3, 't=4': lambda k: k[2] == 4, 't=5': lambda k: k[2] == 5,
                    't>=6': lambda k: k[2] >= 6}
    lot_strata = {'m-t=1': lambda k: k[3] - k[2] == 1, 'm-t=2': lambda k: k[3] - k[2] == 2,
                  'm-t=3': lambda k: k[3] - k[2] == 3, 'm-t>=4': lambda k: k[3] - k[2] >= 4}
    # eligibility: within the memory and time range of all compared programs, not trivial, not held by us
    plain = [k for k, b in live.items() if k[2] == k[3] and k[2] >= 3 and k[1] < k[0] and 10 <= b <= 2000
             and comb(k[0], k[2]) <= 3e7 and k not in held]
    lot = [k for k, b in live.items() if k[2] < k[3] and k[2] >= 2 and k[1] < k[0] and 10 <= b <= 2000
           and comb(k[0], k[3]) <= 3e7 and k not in held]
    rng = random.Random(SEED)
    for name, pool, strata in (('plain', plain, plain_strata), ('lottery', lot, lot_strata)):
        with open(os.path.join(outdir, 'order_%s.tsv' % name), 'w', newline='\n') as f:
            f.write('stratum\trank\tv\tk\tt\tm\ttable_value\n')
            for s, test in strata.items():
                members = sorted(k for k in pool if test(k)); rng.shuffle(members)
                for i, k in enumerate(members):
                    f.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\n' % (s, i + 1, *k, live[k]))
                print(name, s, 'eligible', len(members))
    for p in (vkt, vktm): print(os.path.basename(p), 'sha1', hashlib.sha1(open(p, 'rb').read()).hexdigest())
    print('excluded as ours:', len(held))


if __name__ == '__main__':
    main(*sys.argv[1:4])
