# Covering design records by weighted local search

This repository holds the code behind **233 current records for covering numbers C(v,k,t)** on the [Covering Repository](https://www.coveringrepository.com) (October 2026). [RECORDS.md](RECORDS.md) lists every entry, its previous value, and which solver found it.

A (v,k,t) covering design is a set of k-element blocks drawn from v points such that every t-element subset of the points lies in at least one block. C(v,k,t) is the smallest number of blocks possible. For most parameters the exact value is unknown, so the public tables record the best covering anyone has found.

## Repository layout

| Folder | What it is | Status |
|---|---|---|
| [`v1/`](v1/) | `wls128`: the first solver. Every record in RECORDS.md so far came from it. | Kept unchanged so the published records can be reproduced. |
| [`v2/`](v2/) | `wls2`: same search, 11–18× faster thanks to incremental move scores. | Current solver. |
| [`vktm/`](vktm/) | `wlsm`: experimental solver for lottery coverings C(v,k,t,m). | Correct but not yet competitive (see below). |
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

All five v2 covers beat our own v1 records. A full v2 rerun over about 3,250 entries is in progress. Records it produces will be marked `v2` in RECORDS.md once they are accepted.

v1's batch runs also stopped each entry after shrinking 40 blocks below the starting value. That cap is gone in v2. It had cut off C(54,13,4) at 706, which v2 takes to 666.

## Lottery coverings (vktm, experimental)

A (v,k,t,m) covering asks that every m-subset meet some block in at least t points: the lottery guarantee "if the m drawn numbers are among your v, some ticket matches t". `wlsm` applies the same weighted search with weights on m-subsets. A swap x_out → x_in in block B only changes m-subsets of the form S ∪ {x} ∪ R, where S is a (t−1)-subset of B and R an (m−t)-subset outside B. The gain and loss of every swap for a block come from one pass over those sets.

Status, stated honestly: every cover `wlsm` produces passes the independent checker `mverify`, but it does not yet beat the lottery records. In a first test it could not remove a single block from any of the record covers we tried. Those records are mostly built with dedicated lotto-design software. We are also finding that from-scratch construction is the weak point of local search in general, so we expect structured constructions to matter more here than raw search speed.

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

## Credits

Search, verification and uploads: Wynand van Staden. The code was written with AI assistance (Claude). The coverings themselves are hosted by the Covering Repository; each was checked with `scripts/verify.py` and by the site's Validator before upload.

## License

MIT
