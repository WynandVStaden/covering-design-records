"""Load current (v,k,t) records and compute Schonheim bounds / slack."""
from math import comb, ceil
def schonheim(v, k, t):
    L = 1
    for i in range(t - 1, -1, -1):
        L = ceil((v - i) / (k - i) * L)
    return L
def load(path="data/vkt_archive.bin"):
    rec = {}
    for line in open(path):
        p = line.split()
        if len(p) < 5: continue
        v, k, t, m, b = map(int, p[:5])
        if m == t: rec[(v, k, t)] = b
    return rec
if __name__ == "__main__":
    rec = load()
    print(len(rec), "entries")
    rows = []
    for (v, k, t), b in rec.items():
        L = schonheim(v, k, t)
        rows.append((b / L, v, k, t, b, L, comb(v, k)))
    # tractable for orbit ILP: C(v,k) <= 3e7
    cand = sorted([r for r in rows if r[6] <= 3e7 and r[3] >= 3], reverse=True)
    print("tractable entries:", len(cand))
    print("largest ratio record/Schonheim among tractable:")
    for r in cand[:25]:
        print(f"  C({r[1]},{r[2]},{r[3]}) = {r[4]}  Schonheim {r[5]}  ratio {r[0]:.3f}  C(v,k)={r[6]:.2e}")
