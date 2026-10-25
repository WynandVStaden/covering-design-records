# Random-sample comparison: protocol (fixed before any run, 2026-10-25)

## Question

On table entries chosen at random (not ones we picked), how often does each method improve the current record
within a fixed budget, and by how much? Which entry features go with success?

## Population and sample

- Source: the Covering Repository table archives downloaded 2026-10-24 23:06 (sha1 `6bd4f21b…` for plain,
  `538d18d9…` for lottery; full hashes are printed by `draw.py`).
- An entry is one parameter tuple (v,k,t,m); plain covering designs have m = t.
- Eligible: t ≥ 3 (plain) or t ≥ 2 (lottery), k < v, table value between 10 and 2,000 blocks, C(v,t) ≤ 3·10^7
  (plain) or C(v,m) ≤ 3·10^7 (lottery), so that every compared program can run it; and not currently held by one of
  our coverings (998 excluded).
- Strata: plain by t (3, 4, 5, ≥6); lottery by m − t (1, 2, 3, ≥4).
- `draw.py` (seed 20261025) writes, per stratum, every eligible entry in a random order (`order_plain.tsv`,
  `order_lottery.tsv`).
- The sample is the first 25 entries of each stratum's order for which a start cover is available, in that order:
  100 plain and 100 lottery entries. Entries skipped for lack of a start cover are listed in the results.

## Start covers

- Plain: the covering in the La Jolla Covering Repository archive (Gordon, Zenodo 19735294, CC BY 4.0). Where the
  Covering Repository's value is smaller, the run must first reach it; success is still measured against the live
  table value.
- Lottery: the current record covering downloaded from the Covering Repository (for personal use; not redistributed).

## Methods (same decremental framework, single-threaded)

- Plain: `weighted` (wls2), `tabu` (wls2, NOWEIGHT=1), `annealing` (sa_cover, default settings).
- Lottery: `weighted` (wlsm2, 1,000 candidates, tenure 5), `tabu` (same, NOWEIGHT=1), `direct` (wlsm, same settings).

## Budget

120 seconds of search per run, 3 seeds per method and entry, 11 runs at a time on a 6-core / 12-thread machine,
runs in shuffled order.

## Outcomes and analysis (fixed in advance)

1. Per entry and method: best size over seeds, and median size.
2. Primary: share of entries where the method's best size is below the live table value (improvement rate), per
   stratum and overall, with the number of entries.
3. Secondary: median relative gap (size − table value) / table value; paired win/tie/loss counts between methods
   on median size.
4. Features examined for success: stratum, v, k, table value, and for plain entries the ratio of the table value to
   the Schönheim bound. Descriptive only (no model selection).
5. Every covering counted as an improvement is re-checked by `scripts/verify.py` or `vktm/mverify`.

Hypotheses from our earlier records (stated before the run): weighted beats tabu and annealing on plain entries;
on lottery entries the t-level methods (weighted, tabu) beat direct scoring; on lottery entries weighting is not
clearly better than tabu.
