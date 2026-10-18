# Covering design records by weighted local search

This repository holds the code that found **185 new upper bounds on covering numbers C(v,k,t)**, all accepted by the [Covering Repository](https://www.coveringrepository.com) in October 2026. The full list is in [RECORDS.md](RECORDS.md).

A (v,k,t) covering design is a set of k-element blocks drawn from v points such that every t-element subset of the points lies in at least one block. C(v,k,t) is the smallest number of blocks possible. For most parameters the exact value is unknown, and the public tables record the best covering anyone has found.

| t | records |
|---|---|
| 3 | 141 |
| 4 | 26 |
| 5 | 10 |
| 6 | 6 |
| 7 | 2 |

The largest gains include:
- C(54,13,4): 746 → 706
- C(70,13,4): 2580 → 2540
- C(62,13,4): 1513 → 1474
- C(96,9,3): 2270 → 2248
- C(48,9,3): 292 → 274

## Method

The search borrows its main idea from SAT and MaxSAT solvers: **clause weighting**. It treats every t-subset as a clause that must be satisfied (covered). It is implemented in [`src/wls128.cpp`](src/wls128.cpp).

1. **State.** The state is b blocks, with a weight on every t-subset (initially 1) and a count of how many blocks cover each t-subset.
2. **Move.** Pick an uncovered t-subset T: usually the heaviest of a small random sample, sometimes a random one. Among blocks that already contain t−1 points of T, consider swapping one of their other points for the missing point. Choose the swap with the best weighted gain (weight of newly covered t-subsets minus weight of newly uncovered ones). Ties are broken at random, and a short tabu list forbids undoing recent swaps.
3. **Weighting.** After every move, add 1 to the weight of every t-subset that is still uncovered. Subsets that stay uncovered become expensive to leave uncovered, which pushes the search out of local minima. This is the "breakout" idea from SAT local search and the RWLS algorithm for unicost set cover (Gao et al., 2015).
4. **Shrink mode.** Start from the current record covering. Whenever the cover is complete, save it, delete the block that covers the fewest t-subsets uniquely, and let the weighted search repair the hole. Repeat until time runs out. Every complete cover reached along the way is a candidate record.

Blocks are 128-bit masks (so v ≤ 128). t-subsets are indexed by their combinatorial rank, so each move costs time proportional to the blocks touched, not to C(v,t).

On our test entries this was 200–700 times faster than the parallel tempering and simulated annealing we tried first, and those methods found no records. The wins came almost entirely from entries that other people had improved recently (2025–26), mostly at t = 3 with large v. The oldest table entries (1996–2005) did not improve with any method we tried.

## Usage

Compile the solver (needs GCC or Clang for `unsigned __int128`; on Windows use WSL):

```bash
g++ -O3 -march=native -o wls128 src/wls128.cpp
```

Shrink a cover. `start.txt` is a valid C(48,9,3) covering with 292 blocks, one block per line, points 1..48. The search runs for 600 seconds with seed 1 and goes no lower than 270 blocks:

```bash
./wls128 48 9 3 292 600 1 start.txt 7 0.05 270
```

Run many entries in parallel. Each line of the jobs file is `v k t b start_file`:

```bash
bash scripts/run_batch.sh 150 jobs.txt
```

Verify any covering independently (pure Python, brute force over all t-subsets):

```bash
python scripts/verify.py 48 9 3 cover.txt
```

`scripts/table.py` reads the table archive that the Covering Repository publishes and computes the Schönheim lower bound.

## Credits

Search, verification and uploads: Wynand van Staden. The code was written with AI assistance (Claude). The coverings themselves are hosted by the Covering Repository; every one was checked there and with `scripts/verify.py` before upload.

## License

MIT
