"""Rebuild a derived covering from its recipe and the component coverings, checking every construction step.

usage: python derived/build.py v,k,t,m [out.txt]
       python derived/build.py --all          (rebuild and check every recipe, write nothing)

A recipe (derived/recipes/v,k,t,m.json) is a tree. Its leaves are component coverings in derived/components/,
each also published on the Covering Repository. Every inner node is one construction step whose hypothesis is
asserted while rebuilding (M = any m-subset of the points):

  product  blocks B1|B2 over D1 x D2, D1 on points 1..v1, D2 on v1+1..v. If |M & V1| >= m1 and |M & V2| >= m2,
           the D1 and D2 blocks meet M in >= t1 and >= t2 points, so B1|B2 meets M in >= t1+t2 >= t points.
           The terms' intervals [m1, m-m2] must contain every feasible a = |M & V1|.
  dual     complements of a C(v, v-k, v-k-m+t, v-m) covering, since |B^c & M^c| = v-k-m+|B & M|.
  m-1      a covering for (m-1)-sets also covers m-sets.
  k-1      blocks of size k-1, padded to size k.
  t+1      meeting M in t+1 points implies meeting it in t.
  v+1      delete point v+1 from a (v+1)-point covering.
  v-1m-1   an m-set has at least m-1 points among the first v-1.
  trivial  t=0: one empty block; m>v: nothing to cover; v<=k: the whole point set; t=1: blocks covering v-m+1 points.

The result is a valid C(v,k,t,m) covering whenever every component is (check those with scripts/verify.py or
vktm/mverify). Plain coverings are the case m = t.
"""
import sys, os, json, glob
from itertools import product as iprod

HERE = os.path.dirname(os.path.abspath(__file__))


def read_cover(path):
    return [tuple(sorted(int(x) for x in l.split())) for l in open(path) if l.strip()]


def derive(n):
    v, k, t, m, op = n['v'], n['k'], n['t'], n['m'], n['op']
    assert 0 <= t <= m, n
    if op == 'trivial':
        if t == 0: return {()} if m <= v else set()
        if m > v: return set()
        if v <= k: return {tuple(range(1, v + 1))}
        if t == 1:
            need = v - m + 1
            bl = set(tuple(range(i + 1, min(i + k, need) + 1)) for i in range(0, need, k))
            assert len(set().union(*bl)) >= need and all(len(b) <= k for b in bl)
            return bl
        raise AssertionError('bad trivial node %s' % n)
    if op == 'leaf':
        bl = set(read_cover(os.path.join(HERE, n['src'])))
        assert all(len(b) == k and len(set(b)) == k and 1 <= b[0] and b[-1] <= v for b in bl), n['src']
        return bl
    if op in ('m-1', 'k-1', 't+1', 'v+1', 'v-1m-1'):
        c = n['child']
        want = {'m-1': (v, k, t, m - 1), 'k-1': (v, k - 1, t, m), 't+1': (v, k, t + 1, m),
                'v+1': (v + 1, k, t, m), 'v-1m-1': (v - 1, k, t, m - 1)}[op]
        assert (c['v'], c['k'], c['t'], c['m']) == want and want[2] <= want[3], (op, n['v'], c['v'])
        bl = derive(c)
        if op == 'v+1': bl = set(tuple(p for p in b if p <= v) for b in bl)
        return bl
    if op == 'dual':
        ch = n['child']
        assert (ch['v'], ch['k'], ch['t'], ch['m']) == (v, v - k, max(0, v - k - m + t), v - m)
        full = set(range(1, v + 1)); out = set()
        for b in derive(ch):
            s = set(b); p = 1
            while len(s) < v - k:
                s.add(p); p += 1
            out.add(tuple(sorted(full - s)))
        return out
    if op == 'product':
        v1 = n['v1']; v2 = v - v1
        assert 1 <= v1 < v
        cov = set(); out = set()
        for (k1, t1, m1, k2, t2, m2, c1, c2) in n['terms']:
            assert k1 + k2 == k and t1 + t2 >= t and min(k1, k2, t1, t2) >= 0
            assert (c1['v'], c1['k'], c1['t'], c1['m']) == (v1, k1, t1, m1)
            assert (c2['v'], c2['k'], c2['t'], c2['m']) == (v2, k2, t2, m2)
            D1 = derive(c1); D2 = derive(c2)
            for a in range(max(0, m - v2), min(m, v1) + 1):
                if a >= m1 and m - a >= m2 and (D1 or m1 > v1) and (D2 or m2 > v2): cov.add(a)
            for b1, b2 in iprod(D1, D2): out.add(tuple(sorted(b1 + tuple(p + v1 for p in b2))))
        need = set(range(max(0, m - v2), min(m, v1) + 1))
        assert need <= cov, ('interval gap', v, k, t, m, sorted(need - cov))
        return out
    raise AssertionError('unknown step %s' % op)


def pad(bl, k):
    out = set()
    for b in bl:
        s = set(b); p = 1
        while len(s) < k:
            s.add(p); p += 1
        out.add(tuple(sorted(s)))
    return out


def build(name):
    v, k, t, m = map(int, name.split(','))
    tree = json.load(open(os.path.join(HERE, 'recipes', name + '.json')))
    assert (tree['v'], tree['k'], tree['t'], tree['m']) == (v, k, t, m)
    return sorted(pad(derive(tree), k))


if __name__ == '__main__':
    if sys.argv[1] == '--all':
        for p in sorted(glob.glob(os.path.join(HERE, 'recipes', '*.json'))):
            name = os.path.basename(p)[:-5]
            print(name, len(build(name)), 'blocks, all construction steps check out')
    else:
        bl = build(sys.argv[1])
        out = open(sys.argv[2], 'w') if len(sys.argv) > 2 else sys.stdout
        for b in bl: out.write(' '.join(map(str, b)) + '\n')
        print('C(%s): %d blocks, all construction steps check out' % (sys.argv[1], len(bl)), file=sys.stderr)
