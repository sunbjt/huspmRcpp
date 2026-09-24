Production checklist
====================

1. Validate the exact business definition of sequence and utility.
2. Validate the pattern utility semantics on hand-calculated examples.
3. Compare HUSPM against a brute-force enumerator on random tiny datasets.
4. Benchmark memory and runtime on 10^4, 10^5, 10^6+ events.
5. Add max-gap and max-span constraints for retail sessions.
6. Add context partitioning/indexing.
7. Add parallel root-level mining only after deterministic single-threaded results are stable.
8. Consider closed/maximal utility pattern reduction if output explodes.
9. Add serialization/versioning of item dictionaries.
10. Keep a regression suite with fixed expected pattern utilities.
