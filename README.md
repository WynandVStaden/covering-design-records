# Covering design records by weighted local search

This repository holds the code behind **401 current records for covering numbers C(v,k,t)** (128 from v1, 273 from v2) on the [Covering Repository](https://www.coveringrepository.com) (October 2026). [RECORDS.md](RECORDS.md) lists every entry, its previous value, and which solver found it.

A (v,k,t) covering design is a set of k-element blocks drawn from v points such that every t-element subset of the points lies in at least one block. C(v,k,t) is the smallest number of blocks possible. For most parameters the exact value is unknown, so the public tables record the best covering anyone has found.

## Repository layout

| Folder | What it is | Status |
|---|---|---|
| [`v1/`](v1/) | `wls128`: the first solver. Every record in RECORDS.md so far came from it. | Kept unchanged so the published records can be reproduced. |
| [`v2/`](v2/) | `wls2`: same search, 11–18× faster thanks to incremental move scores. | Current solver. |
| [`vktm/`](vktm/) | `wlsm2` (and the first version `wlsm`): solver for lottery coverings C(v,k,t,m). | Finds records on a targeted subset of the lottery table (see below). |
| [`scripts/`](scripts/) | `verify.py` (independent checker for (v,k,t)), `table.py` (reads the site's table archive). | |

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

All five v2 covers beat our own v1 records. A full v2 rerun over 3,251 entries then improved 312 table entries: 188 new records and 124 improvements of v1 records. They are marked `v2` in RECORDS.md.

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

So far this has produced 54 lottery improvements (e.g. C(78,9,3,4) 672 → 647, C(74,8,3,4) 784 → 765). They are not uploaded yet.

**Limits.** wlsm2 stores about 14 bytes per m-subset, so it needs C(v,m) up to roughly 10^9 at most on a 16 GB machine. Most lottery records have far larger m: the median entry held by the table's largest contributor is C(52,5,5,19), with 7.6 × 10^13 m-subsets. Those need a method that never enumerates m-subsets (for example lazy constraint generation with a fast exact checker), which we have not built.

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
