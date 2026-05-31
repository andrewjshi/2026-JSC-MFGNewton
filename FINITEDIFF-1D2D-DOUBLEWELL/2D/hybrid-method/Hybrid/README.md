# Nh=128 hybrid run — results and timing caveat

`hybrid-2d.py 128` was run to completion: the single-grid Newton-Picard → Global-Newton
hybrid at `Nx=Ny=128`, `Nt=10`, with the **direct** linear solver throughout. Output:
`hybrid2d_per_nu_times_Nh128.csv`.

## Result

Total wall time ≈ **31,433 s ≈ 8.7 h**, reaching ν = 1e-4.

| phase | what | time (s) |
|---|---|---|
| 1 | Newton-Picard descent to ν_limit = 0.01 | 1,407 |
| 2 | Newton-Picard fails at ν = 0.003 | 746 |
| 3 | Global Newton handoff re-solve at ν = 0.01 | 6,867 |
| 4 | Global Newton continuation 0.01 → 1e-4 | 22,414 |
| | **TOTAL** | **31,433** |

Global Newton with the direct solver dominates (~93% of the total); each GN solve took ~1–2 h.

## Timing caveat: these numbers are inflated by background load

Phase 4 (GN continuation 0.01 → 1e-4) came in at **22,414 s** vs the **16,671 s** projected
from the paper's Table 5 cumulative differences — **~34% over**. Per-ν, the GN solves ran
~15–50% slower than the Table 5 per-cell values:

| ν | Table 5 implied (s) | this run (s) | ratio |
|---|---|---|---|
| 0.01 | 5,434 | 6,867 | 1.26× |
| 0.003 | 4,817 | 7,065 | 1.47× |
| 0.001 | 4,469 | 6,335 | 1.42× |
| 3e-4 | 3,883 | 4,501 | 1.16× |
| 1e-4 | 3,502 | 4,512 | 1.29× |

**The most likely cause is contention, not the algorithm.** Table 5 was almost certainly
measured on an idle machine. This run instead shared an 8-core / 16-GB laptop with ~8 hours of
concurrent activity (analysis work, repeated `latexmk` compiles, small test runs). The Nh=128
direct factorization is **memory-bandwidth-bound**, so anything else touching memory slows it.

**The evidence points at contention and rules out thermal throttling.** The per–Newton-iteration
cost (which normalizes out the differing iteration counts and isolates the factorization cost)
*dropped* over the course of the run:

```
1373 → 1413 → 1267 → 1125 → 1128     (s per Newton iteration)
```

The slowest cells (ν = 0.003, 0.001, at ~1.4–1.5×) landed during the busiest interactive
stretch; the fastest ones (3e-4, 1e-4) during a quieter, discussion-heavy tail. Thermal
throttling would make the per-iteration cost *increase* over time — here it decreased.

## Takeaways

1. **These absolute numbers are not clean timing data** — they are inflated by background load.
   Do not quote the 8.7 h as *the* Nh=128 direct cost; the paper's Table 5 is the cleaner
   reference.
2. **For any timing comparison to be fair, run it on an otherwise-idle machine** — don't pile
   other heavy work on top of a benchmark while it times.
