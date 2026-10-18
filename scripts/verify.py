"""Independent brute-force check of a covering file.

usage: python verify.py v k t file.txt

A file is one block per line, points 1..v separated by spaces. The check confirms that every block is a set of
k distinct points in 1..v, that no block repeats, and that every t-subset of {1..v} lies in at least one block.
"""
import sys, itertools
from math import comb


def verify(path, v, k, t):
    blocks = [tuple(sorted(map(int, l.split()))) for l in open(path) if l.strip()]
    if len(set(blocks)) != len(blocks):
        return False, len(blocks), "duplicate blocks"
    if not all(len(set(b)) == k and all(1 <= p <= v for p in b) for b in blocks):
        return False, len(blocks), "block with wrong size or point out of range"
    covered = set()
    for b in blocks:
        covered.update(itertools.combinations(b, t))
    missing = comb(v, t) - len(covered)
    return missing == 0, len(blocks), f"{missing} t-subsets uncovered" if missing else "all t-subsets covered"


if __name__ == "__main__":
    v, k, t = map(int, sys.argv[1:4])
    ok, n, msg = verify(sys.argv[4], v, k, t)
    print(f"C({v},{k},{t}) with {n} blocks: {'VALID' if ok else 'INVALID'} ({msg})")
    sys.exit(0 if ok else 1)
