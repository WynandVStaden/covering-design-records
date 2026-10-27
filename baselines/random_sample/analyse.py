"""Analysis of one half of the random sample (written after the runs; the protocol fixed outcomes, not tests). usage: analyse.py outdir"""
import sys, math, statistics as st
from collections import defaultdict
out = sys.argv[1]
R = defaultdict(lambda: defaultdict(list)); TV = {}; SS = {}; STR = {}
for l in open(out + '/results.tsv').read().splitlines()[1:]:
    x = l.split('\t')
    if len(x) != 7: continue   # garbled by the old concurrent-write bug; that run was repeated
    e, s, m, seed, tv, start, b = x
    R[e][m].append(int(b)); TV[e] = int(tv); SS[e] = int(start); STR[e] = s
meths = sorted({m for e in R for m in R[e]}, key=lambda m: ['weighted', 'tabu', 'annealing', 'direct'].index(m))
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z*z/n; c = (p + z*z/(2*n)) / d; h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return '%d/%d = %.1f%% (%.1f-%.1f)' % (k, n, 100*p, 100*(c-h), 100*(c+h))
def binom_two(k, n):
    if n == 0: return 1.0
    x = min(k, n-k); p = sum(math.comb(n, i) for i in range(x+1)) / 2**n
    return min(1, 2*p)
print('entries', len(R), 'start<table', sum(SS[e] > TV[e] for e in R), 'start==table', sum(SS[e] == TV[e] for e in R))
strata = sorted(set(STR.values()))
for m in meths:
    imp = [e for e in R if min(R[e][m]) < TV[e]]
    print('%-9s improved %s' % (m, wilson(len(imp), len(R))), ' | ', '  '.join('%s %d' % (s, sum(STR[e] == s for e in imp)) for s in strata))
    gaps = [(st.median(R[e][m]) - TV[e]) / TV[e] * 100 for e in R if TV[e] > 0]
    imps = [(TV[e] - min(R[e][m])) / TV[e] * 100 for e in imp]
    print('          median rel gap of median %.3f%%; mean improvement among improved %.2f%% (median %.2f%%)' % (st.median(gaps), st.mean(imps) if imps else 0, st.median(imps) if imps else 0))
any_ = [e for e in R if any(min(R[e][m]) < TV[e] for m in meths)]
print('any method improved', wilson(len(any_), len(R)))
for a in meths:
    for b in meths:
        if a >= b and not (a == 'weighted'): continue
        if a == b: continue
        w = sum(st.median(R[e][a]) < st.median(R[e][b]) for e in R); l = sum(st.median(R[e][a]) > st.median(R[e][b]) for e in R)
        w2 = sum((min(R[e][a]) < TV[e]) and not (min(R[e][b]) < TV[e]) for e in R); l2 = sum((min(R[e][b]) < TV[e]) and not (min(R[e][a]) < TV[e]) for e in R)
        print('%s vs %s: median wins %d ties %d losses %d sign p=%.2g | improve-only %d vs %d McNemar p=%.2g' % (a, b, w, len(R)-w-l, l, binom_two(w, w+l), w2, l2, binom_two(w2, w2+l2)))
for s in strata:
    es = [e for e in R if STR[e] == s]
    print(s, '  '.join('%s %s' % (m, wilson(sum(min(R[e][m]) < TV[e] for e in es), len(es))) for m in meths))
# descriptive features, best method = weighted
def feat(name, f, cuts):
    row = []
    for lo, hi in cuts:
        es = [e for e in R if lo <= f(e) < hi]
        row.append('[%s,%s) %d/%d' % (lo, hi, sum(min(R[e]['weighted']) < TV[e] for e in es), len(es)))
    print(name, '; '.join(row))
P = lambda e: tuple(map(int, e.split(',')))
feat('v', lambda e: P(e)[0], [(0, 30), (30, 50), (50, 70), (70, 999)])
feat('k', lambda e: P(e)[1], [(0, 6), (6, 10), (10, 16), (16, 99)])
feat('t', lambda e: P(e)[2], [(0, 3), (3, 4), (4, 6), (6, 99)])
feat('size', lambda e: TV[e], [(0, 50), (50, 200), (200, 800), (800, 10**9)])
feat('start==table', lambda e: int(SS[e] == TV[e]), [(0, 1), (1, 2)])
def schonheim(v, k, t):
    b = 1
    for i in range(t - 1, -1, -1): b = -(-b * (v - i) // (k - i))
    return b
if all(P(e)[2] == P(e)[3] for e in R):   # plain: ratio of table value to the Schonheim bound
    feat('table/Schonheim x10', lambda e: 10 * TV[e] / schonheim(*P(e)[:3]), [(0, 11), (11, 12), (12, 14), (14, 16), (16, 999)])
