## huspmRcpp demo.
##
## Run with:  Rscript inst/examples/demo.R
## or, once installed:
##   Rscript "$(Rscript -e 'cat(system.file("examples/demo.R", package="huspmRcpp"))')"

library(huspmRcpp)

rule <- function(title) cat("\n", strrep("-", 68), "\n", title, "\n",
                            strrep("-", 68), "\n", sep = "")

## ---------------------------------------------------------------------------
## 1. A tiny example you can check by hand
## ---------------------------------------------------------------------------
##
##   seq 1:  A(5)  B(2)  C(10)
##   seq 2:  A(4)  B(3)  D(8)
##   seq 3:  A(6)  C(7)
##
## u(P) is the SUM OVER SEQUENCES of the MAX-UTILITY EMBEDDING of P in that
## sequence. So for A -> C:  seq 1 takes A(5)+C(10) = 15, seq 3 takes
## A(6)+C(7) = 13, and seq 2 has no C at all. Total 28.

rule("1. hand-checkable example")

x <- data.frame(
  sequence_id = c(1, 1, 1,  2, 2, 2,  3, 3),
  item_id     = c("A", "B", "C",  "A", "B", "D",  "A", "C"),
  utility     = c(5, 2, 10,  4, 3, 8,  6, 7)
)

print(x)

res <- huspm(x, min_utility = 0, max_length = 3)
print(res)

## Check a few by hand:
u <- function(p) res$utility[res$pattern == p]
stopifnot(
  # A appears in all three sequences: 5 + 4 + 6
  u("A") == 15,
  # A -> C: (5 + 10) in seq 1, (6 + 7) in seq 3, nothing in seq 2
  u("A -> C") == 28,
  # A -> B -> C only fits in seq 1: 5 + 2 + 10
  u("A -> B -> C") == 17
)

## ---------------------------------------------------------------------------
## 2. The part that surprises people: u(A) is NOT the sum of A's utilities
## ---------------------------------------------------------------------------
##
## A single sequence with A repeated. u(A) is the best SINGLE embedding, i.e.
## the largest one: 3, not 1+2+3 = 6. A -> A uses two distinct positions, and
## the best pair is (2, 3) = 5.

rule("2. repeated items: u(P) takes the best embedding, not the sum")

y <- data.frame(sequence_id = 1, item_id = "A", utility = c(1, 2, 3))

print(huspm(y, min_utility = 0, max_length = 3))

## ---------------------------------------------------------------------------
## 3. Ordering and top_k
## ---------------------------------------------------------------------------
##
## Results always come back in one deterministic order: utility descending,
## then length ascending, then pattern lexicographic ascending.
##
## Here C and A -> B are exactly tied at 100. The shorter pattern wins the
## tie, so top_k = 1 keeps C.

rule("3. deterministic order and top_k")

z <- data.frame(
  sequence_id = c(1, 1,  2),
  item_id     = c("A", "B",  "C"),
  utility     = c(60, 40,  100)
)

cat("full result:\n")
print(huspm(z, min_utility = 0, max_length = 2))

cat("\ntop_k = 1 (the tied pair is broken by length, so C wins):\n")
print(huspm(z, min_utility = 0, top_k = 1, max_length = 2))

## ---------------------------------------------------------------------------
## 4. Ordering events by a time column
## ---------------------------------------------------------------------------
##
## Rows are used in the order given unless you name a time column. With
## time = "time" the sequence is re-sorted, so B -> A becomes impossible and
## A -> C appears.

rule("4. time column controls within-sequence order")

w <- data.frame(
  sequence_id = c(1, 1, 1),
  item_id     = c("B", "A", "C"),
  utility     = c(1, 5, 10),
  time        = c(3, 1, 2)
)

print(huspm(w, min_utility = 0, max_length = 2, time = "time"))

## ---------------------------------------------------------------------------
## 5. A larger run with a planted pattern
## ---------------------------------------------------------------------------
##
## 2000 sequences over 20 items. About a third of them are "loyal customers"
## who buy coffee, then milk, then sugar, at high utility; the rest of each
## basket is random low-utility noise. A correct miner should recover the
## planted trio as the single top pattern, with support close to the planting
## rate.

rule("5. planted pattern over 2000 sequences")

set.seed(7)
pool <- c("bread", "eggs", "butter", "cheese", "apple", "yogurt", "rice",
          "pasta", "tomato", "onion", "chicken", "beef", "salmon", "juice",
          "tea", "water", "chips")
planted <- c("coffee", "milk", "sugar")

make_basket <- function(s) {
  noise_n <- sample(2:7, 1)
  noise <- data.frame(sequence_id = s,
                      item_id = sample(pool, noise_n, replace = TRUE),
                      utility = sample(5:25, noise_n, replace = TRUE))
  if (runif(1) >= 0.35) return(noise)
  rbind(data.frame(sequence_id = s, item_id = planted,
                   utility = sample(30:60, length(planted))),
        noise)
}

sim <- do.call(rbind, lapply(1:2000, make_basket))

cat("events:", nrow(sim), " sequences:", length(unique(sim$sequence_id)),
    " items:", length(unique(sim$item_id)), "\n\n")

t0 <- proc.time()[["elapsed"]]
top <- huspm(sim, min_utility = 0, top_k = 10, max_length = 4)
cat("top_k = 10, max_length = 4  ->", proc.time()[["elapsed"]] - t0, "s\n")
print(top)

## The planted trio comes back first, and everything below it is a sub-pattern
## or an extension of it. Note how the three 2-item sub-patterns share the same
## support but differ in utility -- each one drops a different item's utility.

## min_support is a cheap filter: support is monotone, so a pattern that fails
## it prunes its entire subtree. At 500 only the planted family survives.
t0 <- proc.time()[["elapsed"]]
freq <- huspm(sim, min_utility = 0, max_length = 4, min_support = 500)
cat("\nmin_support = 500  ->", proc.time()[["elapsed"]] - t0, "s\n")
print(freq)

## include_singletons = FALSE hides the 1-item patterns, which is what you want
## when you care about *sequential* structure rather than item popularity.
t0 <- proc.time()[["elapsed"]]
seqs_only <- huspm(sim, min_utility = 0, top_k = 8, max_length = 3,
                   include_singletons = FALSE)
cat("\ninclude_singletons = FALSE, top_k = 8  ->",
    proc.time()[["elapsed"]] - t0, "s\n")
print(seqs_only)

cat("\ndone.\n")
