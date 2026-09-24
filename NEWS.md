# huspmRcpp 0.2.0

## Breaking change

- `top_k` now keeps the first K patterns **of the output order** (utility
  descending, then length ascending, then pattern lexicographic ascending),
  instead of whichever K happened to arrive first among tied patterns. With
  `C` and `A -> B` tied at utility 100, `top_k = 1` now returns `C` where it
  previously returned `A -> B`.

## Fixed

- `huspm()` failed on every call with
  `no function 'Rcpp' provided by package 'Rcpp_precious_remove'`. `DESCRIPTION`
  declared `Imports: Rcpp`, but `NAMESPACE` never imported it, so Rcpp's
  namespace and DLL were never loaded and its C-callables stayed unregistered.
  `NAMESPACE` now carries `importFrom(Rcpp, evalCpp)`; `Rcpp::compileAttributes()`
  does not emit that line.
- `huspm()` failed with `Expecting a single value: [extent=0]` whenever `top_k`
  was left at its default: `NULL` was passed straight through to the C++ `int`
  parameter, overriding the generated `0L` default.
- `R CMD check` never ran the test suite: `tests/testthat.R` was missing.

## Changed

- Rewrote the C++ search core. Per node it now scores every candidate in a
  single O(sequence length) walk instead of one `extend_one()` pass per
  candidate (each of which was run twice, once to score and once to expand),
  and it prunes states by Pareto domination. Typical speedups against the
  previous core are 20x-360x.
- The extension upper bound is now `max_length`-aware (`g[p][b]`), strictly
  tighter than the previous suffix sum.
- Top-K uses a bounded heap instead of re-sorting the whole result vector on
  every insertion.
- Results are stored in a flat arena rather than one `vector<int>` per result.
- The R-side pattern-label rebuild is vectorized; it was roughly a third of
  wall-clock time at ~6e5 patterns and is now a few percent.
- Removed the obsolete `CXX_STD = CXX11` from `src/Makevars` and
  `src/Makevars.win`; the code builds under R's default standard.
- Added `.Rbuildignore` so `.claude/` is not packed into the tarball.
- Added `man/` documentation for all exported functions.
- Added `inst/examples/demo.R`, a runnable end-to-end tour of the package.

# huspmRcpp 0.1.1

- Fix missing Rcpp generated binding for `huspm_cpp()`.
- Add explicit `RcppExports.R` and `RcppExports.cpp`.
- Register the native routine with `R_registerRoutines()`.
