# Random-sample comparison

A pre-registered test of how often each search improves a randomly chosen table entry, reported in the section "A random sample of table entries" of the
paper ([../../paper/](../../paper/)). The protocol, the hypotheses and the sample were committed before any run;
see [PROTOCOL.md](PROTOCOL.md) and the git history.

## Files

| File | What it is |
|---|---|
| `PROTOCOL.md` | Population, strata, methods, budget, outcomes and hypotheses, fixed in advance, with Amendment 1 (75 entries per stratum). |
| `draw.py` | Draws the random order of eligible entries per stratum (seed 20261025) from the table archives. |
| `order_plain.tsv`, `order_lottery.tsv` | Every eligible entry, in random order per stratum. |
| `sample_plain300.tsv`, `sample_lottery300.tsv` | The sample: first 75 entries per stratum (`sample_*.tsv` without 300 are the first 25 per stratum). The `start` column is a local path; start covers are downloaded from the Covering Repository and not redistributed. |
| `run.py`, `run_all.sh` | Runs every (entry, method, seed); resumable. `stop.bat` / `resume.bat` pause and continue it on Windows + WSL. |
| `results_plain.tsv`, `results_lottery.tsv` | Raw results, one line per run: entry, stratum, method, seed, table value, start size, best size. A few lines garbled by an early concurrent-write bug are skipped; those runs were repeated. |
| `analyse.py` | The analysis: improvement rates with Wilson intervals, per stratum, paired comparisons, descriptive features. `python analyse.py <dir>` reads `<dir>/results.tsv`. |
| `verify_plain.txt`, `verify_lottery.txt` | Exhaustive checks of every improved covering (best per entry and method): all valid. |

## Results (300 entries per kind, 3 seeds × 120 s)

| | weighted | tabu (no weighting) | annealing / direct | any |
|---|---|---|---|---|
| covering designs improved | 43 (14.3%) | 25 (8.3%) | 1 (0.3%, annealing) | 54 (18.0%) |
| lottery coverings improved | 4 (1.3%) | 13 (4.3%) | 8 (2.7%, direct) | 16 (5.3%) |

On covering designs the weighting helps (29 entries improved only by weighted, 11 only by tabu). On lottery coverings
it hurts (1 against 10). The improved coverings were measured against the tables of 24 October 2026 and have not been
uploaded.
