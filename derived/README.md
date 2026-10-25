# Derived coverings

These coverings are not found by search. They are built from smaller coverings that are already published, by
constructions whose correctness follows from the components. [records/derived.md](../records/derived.md) lists every
derived covering we have published, with the exact components it uses; [records/components.md](../records/components.md)
lists each component and where it comes from.

## Constructions

A lottery covering C(v,k,t,m) is a set of k-subsets (blocks) of v points such that every m-subset M meets some block in
at least t points. Plain coverings C(v,k,t) are the case m = t.

**Product over a split of the points.** Split the points into V1 (v1 points) and V2 (v − v1 points). For coverings
D1 = C(v1,k1,t1,m1) and D2 = C(v−v1,k2,t2,m2) with k1+k2 = k and t1+t2 ≥ t, the *term* D1 × D2 is every union B1 ∪ B2.
If |M ∩ V1| ≥ m1 and |M ∩ V2| ≥ m2, the two blocks meet M in at least t1 and t2 points, so the union meets it in at
least t. A term therefore covers every m-set whose split a = |M ∩ V1| lies in [m1, m − m2]. A set of terms whose
intervals cover every feasible a, from max(0, m − v + v1) to min(m, v1), is a covering of size Σ|D1|·|D2|.
Pieces may themselves be products, so the construction is recursive.

**Complement duality.** B meets M in at least t points exactly when the complement of B meets the complement of M in at
least v − k − m + t points. So the complements of a C(v, v−k, v−k−m+t, v−m) covering form a C(v,k,t,m) covering of the
same size.

**One-step rules.** A covering for (m−1)-sets also covers m-sets; blocks of size k−1 can be padded to size k; meeting M
in t+1 points implies meeting it in t; deleting point v+1 from a C(v+1,k,t,m) covering gives a C(v,k,t,m) covering;
and an m-set has at least m−1 points among the first v−1, so a C(v−1,k,t,m−1) covering also works.

Partition ("sum") constructions and the standard point-addition construction are special cases of the product.

## Files

| Path | Content |
|---|---|
| `recipes/v,k,t,m.json` | the construction tree for one derived covering; leaves point to files in `components/` |
| `components/C_v_k_t_m_b.txt` | not published: download each component from the Covering Repository (for personal use), save it here under this name and check its sha1 against `data/components.tsv` |
| `build.py` | rebuilds a derived covering from its recipe and checks the hypothesis of every step |

## Rebuild and check

```bash
python derived/build.py 81,13,3,4 cover.txt     # rebuild one covering
python derived/build.py --all                   # rebuild every recipe and check every step
```

`build.py` proves the result is a covering *given* that its components are. Every component is published on the
Covering Repository and can be checked on its own (`scripts/verify.py` for plain coverings, `vktm/mverify` for lottery
coverings). Components are small enough for exhaustive checks; most derived coverings are not, which is why the
recipes are published.

## Sources and credit

Components come from our own search records, from the Covering Repository, and from the La Jolla Covering Repository
archive (D. M. Gordon, Zenodo 19735294, CC BY 4.0). Each component's origin is listed in
[records/components.md](../records/components.md). Every component is also available on the Covering Repository at the
same size or smaller, as its rules for derived coverings require.
