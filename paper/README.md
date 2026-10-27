# Paper

**Weighted local search and a construction scan for covering-design tables, with an extension to lottery coverings**

- [`wls_coverings.pdf`](wls_coverings.pdf): the paper.
- [`wls_coverings.tex`](wls_coverings.tex): its LaTeX source (`pdflatex wls_coverings.tex`, run twice).

## Where each part of the paper lives in this repository

| Paper section | Code and data |
|---|---|
| Search for covering designs | [`v1/`](../v1/) (recompute implementation), [`v2/`](../v2/) (`wls2`, incremental scores; `NOWEIGHT=1` switches weighting off) |
| Lottery coverings | [`vktm/wlsm.cpp`](../vktm/wlsm.cpp) (direct m-subset scoring), [`vktm/wlsm2.cpp`](../vktm/wlsm2.cpp) (t-level scoring) |
| Large covers | [`vktm/wlsm2big.cpp`](../vktm/wlsm2big.cpp), [`scripts/reshrink_batch.py`](../scripts/reshrink_batch.py) |
| Construction scan | [`derived/prodscan.cpp`](../derived/prodscan.cpp), [`derived/build.py`](../derived/build.py), recipes in [`derived/recipes/`](../derived/recipes/) |
| Verification | [`scripts/verify.py`](../scripts/verify.py), [`vktm/mverify.cpp`](../vktm/mverify.cpp), [`scripts/zerocov.cpp`](../scripts/zerocov.cpp) |
| Random sample of table entries | [`baselines/random_sample/`](../baselines/random_sample/) (protocol, sample, runner, raw results, analysis); annealing in [`baselines/sa_cover.cpp`](../baselines/sa_cover.cpp) |
| Published coverings | [`RECORDS.md`](../RECORDS.md), [`records/`](../records/), [`data/`](../data/) |

The coverings themselves are on the [Covering Repository](https://www.coveringrepository.com). Component coverings
downloaded from there are identified by parameters, size and sha1 but not redistributed (personal use only).
