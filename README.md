# huspmRcpp

A compact R + Rcpp implementation of **High-Utility Sequential Pattern Mining (HUSPM)**.

## Semantics

This implementation mines ordered **single-item events**. A pattern such as

`A -> B -> C`

means that A occurs before B, and B occurs before C, within the same sequence.

Utilities must be strictly positive.

For a pattern P in one sequence, its utility is defined as the **maximum utility over valid embeddings** of P in that sequence. The total pattern utility is the sum of these per-sequence maximum utilities.

This definition avoids counting multiple embeddings of the same pattern repeatedly.

The implementation uses a safe extension upper bound. Per sequence, the gain
available after a state ending at position `q` with `b` slots of `max_length`
left is `g[q+1][b]`, the best utility of a subsequence of positions `>= q+1` of
length at most `b`:

```
g[p][j] = max(g[p+1][j], util[p] + g[p+1][j-1])
UB(P)   = sum_s max over states (util(state) + g[end_pos+1][max_length - |P|])
```

If UB is below the active threshold, the subtree is pruned. Respecting
`max_length` in the bound makes it strictly tighter than the plain suffix sum,
increasingly so at greater depths.

Within a sequence, a state `(q1, u1)` dominates `(q2, u2)` when `q1 <= q2` and
`u1 >= u2`; dominated states are dropped. After that pruning, `end_pos` and
`utility` both strictly increase, so the best predecessor of any position is
simply the last surviving state before it.

## Ordering and Top-K

Results are always returned in one deterministic order: **utility descending,
then length ascending, then pattern lexicographic ascending**.

`top_k` keeps the first K patterns *of that order*, so ties are resolved by the
same rule as the display rather than by discovery order. For example, with `C`
and `A -> B` tied at utility 100, `top_k = 1` keeps `C`.

## Scale

Measured on this machine (R 4.6, arm64), single-threaded, peak RSS includes R
itself (~85 MB baseline):

| workload | old core | rewritten core |
|---|---|---|
| 11k events, 200 items, `max_length=4`, exhaustive | 98.1 s | **0.27 s** |
| 30k events, 30 items, `max_length=4`, exhaustive, 74k patterns | 6.8 s | **0.30 s** |
| 30k events, 30 items, `max_length=4`, `top_k=10000` | 12.1 s | **0.24 s** |
| 10^5 events, 500 items, `max_length=5`, `top_k=10000` | — | **1.4 s** (108 MB) |
| 10^6 events, 2000 items, `max_length=5`, `top_k=10000` | — | **12.5 s** (246 MB) |
| 10^6 events, 100 items, `max_length=5`, exhaustive, 1.01M patterns | — | **40.7 s** |

Cost per search node is O(events in the database), so **exhaustive** mining at
10^6 events is bounded by the size of the output rather than by the search: the
same 10^6-event database yields 10k patterns in 3.6 s, 1.0M patterns in 40.7 s,
and 1.4M patterns in 254 s. If you need a large result at that scale, use
**Top-K** or raise `min_utility`.

## Example

A runnable tour lives in `inst/examples/demo.R`:

```sh
Rscript "$(Rscript -e 'cat(system.file("examples/demo.R", package="huspmRcpp"))')"
```

```r
library(huspmRcpp)

x <- data.frame(
  sequence_id = c(1,1,1, 2,2,2, 3,3),
  item_id     = c("A","B","C", "A","B","D", "A","C"),
  utility     = c(5,2,10, 4,3,8, 6,7)
)

huspm(
  x,
  min_utility = 20,
  max_length = 3
)

# Top-K
huspm(
  x,
  min_utility = 0,
  top_k = 100,
  max_length = 5
)
```

## Installation

```r
# download and install
remotes::install_local("~/github/huspmRcpp/",
                       dependencies = TRUE,
                       build_vignettes = FALSE)
```

## Production notes

This is intentionally a **single-item event** engine. It does not yet implement:
- simultaneous itemsets `{A,B} -> C`;
- negative utilities;
- time-gap constraints;
- per-pattern custom utility semantics;
- distributed mining.

Item IDs are integer-coded internally (`prepare_huspm_data()` does the mapping),
and the C++ core avoids R recursion and per-result allocation. The worst-case
search space is still exponential in `max_length`, so the output size — not the
event count — is what bounds an exhaustive run; see **Scale** above.

For large retail data, a production roadmap is:
1. pre-aggregate identical events if business semantics allow it;
2. add max-gap / max-span constraints;
3. add closed/maximal utility-pattern pruning if needed;
4. add parallel root-level mining;
5. persist item and pattern dictionaries;
6. run the randomized brute-force-oracle comparison as part of CI.
