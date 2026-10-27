# Covering design records by weighted local search

This repository holds the code behind our records on the [Covering Repository](https://www.coveringrepository.com): plain coverings C(v,k,t) and lottery coverings C(v,k,t,m) found by weighted local search, and derived coverings built from published coverings by explicit constructions. [RECORDS.md](RECORDS.md) has the current counts and links to every list; each entry shows its previous value, its method, and whether it is still the record.

The paper describing the methods, with a pre-registered comparison on a random sample of table entries, is in [`paper/`](paper/).

A (v,k,t) covering design is a set of k-element blocks drawn from v points such that every t-element subset of the points lies in at least one block. C(v,k,t) is the smallest number of blocks possible. For most parameters the exact value is unknown, so the public tables record the best covering anyone has found.

## Repository layout

| Folder | What it is |
|---|---|
| [`v1/`](v1/) | `wls128`: the first solver, kept unchanged so its records can be reproduced. |
| [`v2/`](v2/) | `wls2`: the same search, 11–18× faster thanks to incremental move scores. Current solver. |
| [`vktm/`](vktm/) | `wlsm2big` (current lottery solver, for large covers), `wlsm2` and the first version `wlsm`, and `mverify`. |
| [`derived/`](derived/) | Derived coverings: construction recipes, their component coverings, and `build.py` to rebuild and check them. |
| [`records/`](records/) | The record lists: search records (plain and lottery), derived coverings, and their components. |
| [`data/`](data/) | The same lists as tab-separated files, the source the Markdown is generated from. |
| [`scripts/`](scripts/) | `verify.py` (independent checker for (v,k,t)), `table.py` (reads the site's table archive), `update_records.py` (regenerates the lists against the live tables), `zerocov.cpp` (exhaustive check and zero-coverage cleanup), `reshrink_batch.py` (batch shrinking). |
| [`baselines/`](baselines/) | `sa_cover.cpp` (simulated annealing for comparison) and [`random_sample/`](baselines/random_sample/): the pre-registered random-sample comparison (protocol, sample, runner, raw results, analysis). |
| [`paper/`](paper/) | The paper (LaTeX source and PDF) and a map from its sections to the code and data. |

## Method (v1)

The search borrows its main idea from SAT and MaxSAT solvers: **clause weighting**. It treats every t-subset as a clause that must be covered.

1. **State.** The state is b blocks, with a weight on every t-subset (initially 1) and a count of how many blocks cover each t-subset.
2. **Move.** Pick an uncovered t-subset T: usually the heaviest of a small random sample. Among blocks that already contain t−1 points of T, consider swapping one of their other points for the missing point. Choose the swap with the best weighted gain (weight of newly covered t-subsets minus weight of newly uncovered ones). A short tabu list forbids undoing recent swaps.
3. **Weighting.** After every move, add 1 to the weight of every t-subset that is still uncovered. Subsets that stay uncovered become expensive, which pushes the search out of local minima. This is the "breakout" idea from SAT local search and RWLS for unicost set cover (Gao et al., 2015).
4. **Shrink mode.** Start from the current record covering. Whenever the cover is complete, save it, delete the block that covers the fewest t-subsets uniquely, and let the search repair the hole. Every complete cover reached is a candidate record.

On our test entries v1 was 200–700 times faster than parallel tempering and simulated annealing, and those methods found no records. The wins came almost entirely from entries other people had improved recently (2025–26), mostly t = 3 with large v. The oldest table entries (1996–2005) did not improve.

## What changed in v2

v1 rescored every candidate swap from scratch, enumerating about 2(k−t+1)·C(k−1,t−1) t-subsets per block. v2 keeps the scores up to date incrementally, the way fast SAT solvers do:

- **Break scores.** `brk[block][point]` holds the total weight of t-subsets that the block covers alone and that contain the point. A t-subset's unique owner is recovered in O(1) from the sum of the ids of its covering blocks, so the scores stay exact as coverage changes.
- **Make scores.** For a block and an incoming point, one pass over the block's (t−1)-subsets gives the gain of every possible outgoing point at once.

The search itself is unchanged; only the bookkeeping is. A `-DCHECK` build recomputes everything by brute force every 5,000 steps to confirm the scores stay exact.

Two-minute runs from the same starting covers (previous record covers from the table archive):

| Entry | Start | v1 result | v2 result | v1 steps/s | v2 steps/s |
|---|---|---|---|---|---|
| C(48,9,3) | 292 | 275 | **273** | 5,916 | 66,463 |
| C(54,13,4) | 746 | 668 | **666** | 391 | 7,112 |
| C(70,12,3) | 352 | 343 | **340** | 2,027 | 36,166 |
| C(35,18,6) | 223 | 211 | **207** | 12 | 200 |
| C(90,12,3) | 704 | 683 | **679** | 1,472 | 25,745 |

All five v2 covers beat our own v1 records. A full v2 rerun over 3,251 entries then improved 312 table entries: 188 new records and 124 improvements of v1 records. They are marked `v2` in RECORDS.md. A second pass of longer (600 s) runs, started from our own best covers on the 319 entries that had improved, shrank 249 of them further (1,629 blocks in total).

v1's batch runs also stopped each entry after shrinking 40 blocks below the starting value. That cap is gone in v2. It had cut off C(54,13,4) at 706, which v2 takes to 666.

## Lottery coverings (vktm)

A (v,k,t,m) covering asks that every m-subset meet some block in at least t points: the lottery guarantee "if the m drawn numbers are among your v, some ticket matches t".

**wlsm (first version)** puts the breakout weights on m-subsets and scores a swap by enumerating the m-subsets it affects. It is correct, but on record covers it almost never removed a block (1 improvement in 78 records tried).

**wlsm2** rests on one observation: an m-subset is covered exactly when one of its t-subsets lies inside a block. So coverage is tracked at two levels, `ct[T]` (blocks containing each t-subset) and `cnt[M]` (covered t-subsets inside each m-subset), and the m-subset weights are pushed down to t-subsets:

- `A[T]`: total weight of uncovered m-subsets containing T, i.e. the gain from covering T;
- `Bc[T]`: total weight of m-subsets whose only covered t-subset is T, i.e. the loss from uncovering T.

A swap is then scored with about C(k−1,t−1) table lookups, as in v2, and m-subsets are only touched when a t-subset actually flips between covered and uncovered. Bumping every uncovered weight is a single increment of a global counter: an uncovered m-subset weighs `base + E`, and `A[T] = Abase[T] + E·U[T]`, where `U[T]` counts uncovered m-subsets containing T. A `-DCHECK` build verifies `ct`, `cnt`, `A` and `Bc` against a brute-force recount. Every saved cover is checked independently by `mverify`.

What we measured while tuning it on record covers:

| Change | Effect |
|---|---|
| t-level scores + global weight counter | 10–25× more steps per second than wlsm |
| score every candidate block (not a sample of 12) | about twice as many records found |
| `STALL=25` (stop an entry after 25 s without shrinking) | about 3.4× more entries per hour, few hits lost |
| portfolio of 3 runs (seeds, tabu 2/5/10) | finds records single runs miss |
| two-for-one block exchange (`PAIRX`) | correct, but no extra records in our tests |
| restarting a level with another deletion (`LEVEL_STEPS`) | no measurable gain |

**Where the records are.** Hit rate depends far more on *which* entry is attacked than on search settings. On record covers downloaded from the site:

| Group | Hit rate |
|---|---|
| records set or re-polished in 2025–26 | 0% (0 of 102) |
| older (≤ 2023) records by Terenzio Lazzarotto or imported "Known designs", t = 4, m = 5, v ≥ 40 | 47% (15 of 32) |
| same, t = 3, m = 4, v ≥ 40 | 11% |
| same, t = 3, m = 4, v < 40 | 1% |

This produced 75 lottery records, all verified exactly (e.g. C(78,9,3,4) 672 → 647, C(74,8,3,4) 784 → 765); see [records/search_lottery.md](records/search_lottery.md).

**Weighting on lottery coverings.** On a random sample of 300 lottery entries not held by us (2 minutes, 3 seeds each), wlsm2 with the weighting switched off (`NOWEIGHT=1`) improved 13 entries, wlsm2 as is 4, and wlsm 8. At these settings and a two-minute budget, weighting combined with wlsm2's t-level score estimates did worst; for lottery coverings, try `NOWEIGHT=1`. See [baselines/random_sample/](baselines/random_sample/).

**Limits.** wlsm2 stores about 14 bytes per m-subset, so it needs C(v,m) up to roughly 10^9 at most on a 16 GB machine. Most lottery records have far larger m: the median entry held by the table's largest contributor is C(52,5,5,19), with 7.6 × 10^13 m-subsets. Those need a method that never enumerates m-subsets (for example lazy constraint generation with a fast exact checker), which we have not built.

## wlsm2big: the lottery solver for large covers

`vktm/wlsm2big.cpp` is wlsm2 with the same search and the same shrink mode (start from a full cover, delete the block
whose removal costs least, repair with point swaps, repeat), engineered for covers with up to about 1.6·10^9 m-subsets
and hundreds of thousands of blocks. It also runs plain coverings (m = t), though v2 is leaner for small plain covers.

| Change | Why |
|---|---|
| witness byte per m-subset (sum of local indices of its covered t-subsets) | when an m-subset is covered once, the byte names the covering t-subset, so a 1↔2 change updates one score instead of C(m,t) |
| O(m−t) colex ranking of t-subsets | cheaper bookkeeping on every coverage change |
| OpenMP for setup, score rebuilds, each flip's superset enumeration, the candidate scan and candidate scoring | per-m-subset state is touched once per flip, so only shared t-set scores need atomics |
| deletion pool with early exit | re-score ~2,000 cheapest blocks instead of all of them; stop at score 0 |
| geometric weight bump between levels (`BUMP_GROWTH`) | the same as decaying old breakout weights; without it long runs drift |
| revert by undo log (`LEVEL_BUDGET`, `GROWTH`, `REVERT_SECS`) | when a deletion cannot be repaired, return to the last full cover and delete another block |
| cost guard, exit-save on SIGTERM, `STALL` on real progress | refuses covers where one flip costs > 10^8 updates; never loses the best cover |

A `-DCHECK` build recounts all bookkeeping by brute force every 2,000 steps; every change above was run through it.

Measured on C(35,10,9,12) (834,451,800 m-subsets), starting from a 404,385-block constructed cover:

| Version | Steps per second | Result |
|---|---|---|
| wlsm2 | — | out of memory |
| wlsm2big, parallel setup and pool | ~33 | stalled near 366,000 |
| + witness byte, reverts, geometric bump | ~120 | |
| + parallel flips, scoring and scan, pool early exit | ~270–350 | **~363,000 after 10 minutes** |

`scripts/reshrink_batch.py` runs wlsm2big over a list of covers with a time limit (60 s without progress, 10 minutes
at most) and verifies every improvement with `scripts/zerocov.cpp`, an independent exhaustive checker that also removes
blocks with zero exclusive coverage. Re-shrinking our own published covers this way improved 137 of them (marked
`+wlsm2big` in the record lists), e.g. C(60,6,5,7) 223,464 → 201,801 and C(59,6,5,7) 204,392 → 183,232. Constructed
covers gain most: the product construction covers many m-subsets several times. Search records were already tight.

## Derived coverings

Many more records come from combining published coverings: products over a split of the points, complement duality and
one-point extensions ([derived/README.md](derived/README.md)). A scan of the whole table finds entries that some
combination of smaller published coverings beats, and the covering is then built explicitly. For each one we publish the
recipe and identify every component covering by its parameters, size and sha1 hash; the components themselves are
not redistributed (coverings downloaded from the Covering Repository are for personal use), but each can be downloaded
there and the covering rebuilt with `python derived/build.py v,k,t,m` ([records/components.md](records/components.md)). They are
listed separately from the search records because they come from known constructions, not from search.

## Updating the lists

```bash
python scripts/update_records.py --work <working folder>
```

downloads the live tables, adds every newly uploaded derived covering (recipe and components) from the working folder,
and regenerates RECORDS.md and records/. Search records stay listed after they are beaten, with the current value.

## Usage

Compile (needs GCC or Clang for `unsigned __int128`; on Windows use WSL):

```bash
g++ -O3 -march=native -o v2/wls2 v2/wls2.cpp
```

Shrink a C(48,9,3) cover `start.txt` that has 292 blocks (one block per line, points 1..48). The search runs for 600 seconds with seed 1 and goes as low as it can:

```bash
./v2/wls2 48 9 3 292 600 1 start.txt 7 0.05 1
```

Many entries in parallel (each jobs-file line is `v k t b start_file`):

```bash
bash v2/run_batch.sh 60 jobs.txt
```

Verify any covering independently:

```bash
python scripts/verify.py 48 9 3 cover.txt
```

The v1 commands are the same with `v1/wls128`.

Lottery coverings (example: shrink a C(70,8,3,4) record cover of 703 blocks, full candidate scan, tabu 5, stop after 25 s without progress):

```bash
g++ -O3 -march=native -o vktm/wlsm2 vktm/wlsm2.cpp
STALL=25 ./vktm/wlsm2 70 8 3 4 703 120 1 start.txt 1000 5 1
g++ -O3 -o vktm/mverify vktm/mverify.cpp && ./vktm/mverify 70 8 3 4 cover.txt
```

## Credits

Search, verification and uploads: Wynand van Staden. The code was written with AI assistance (Claude). The coverings themselves are hosted by the Covering Repository; each was checked with `scripts/verify.py` and by the site's Validator before upload.

## License

MIT
