test_that("basic HUSPM works", {
  x <- data.frame(
    sequence_id = c(1,1,1, 2,2, 2, 3,3),
    item_id     = c("A","B","C", "A","B","D", "A","C"),
    utility     = c(5, 2, 10, 4, 3, 8, 6, 7)
  )

  ans <- huspm(x, min_utility = 10, max_length = 3)

  expect_s3_class(ans, "huspm")
  expect_true(all(ans$utility >= 10))
  expect_true(any(ans$pattern == "A -> C"))
})

test_that("utility is the sum over sequences of the best embedding", {
  # A -> C: seq1 best embedding A(5)+C(10)=15, seq2 has no C -> 15 total
  x <- data.frame(
    sequence_id = c(1,1,1, 2,2),
    item_id     = c("A","B","C", "A","B"),
    utility     = c(5, 2, 10, 4, 3)
  )
  ans <- huspm(x, min_utility = 0, max_length = 2)
  expect_equal(ans$utility[ans$pattern == "A -> C"], 15)
  expect_equal(ans$support[ans$pattern == "A -> C"], 1L)
  # A alone appears in both sequences
  expect_equal(ans$utility[ans$pattern == "A"], 9)
  expect_equal(ans$support[ans$pattern == "A"], 2L)
})

test_that("time ordering is respected", {
  x <- data.frame(
    sequence_id = c(1,1,1),
    item_id = c("B","A","C"),
    utility = c(1, 5, 10),
    time = c(3, 1, 2)
  )

  ans <- huspm(x, min_utility = 14, max_length = 2, time = "time")
  expect_true(any(ans$pattern == "A -> C"))
  expect_false(any(ans$pattern == "B -> A"))
})

test_that("ties on time keep the input order", {
  x <- data.frame(
    sequence_id = c(1,1,1),
    item_id = c("A","B","C"),
    utility = c(1, 2, 3),
    time = c(1, 1, 1)
  )
  ans <- huspm(x, min_utility = 6, max_length = 3, time = "time")
  expect_true(any(ans$pattern == "A -> B -> C"))
})

test_that("results are ordered by utility desc, length asc, pattern asc", {
  x <- data.frame(
    sequence_id = c(1,1, 2, 3,3),
    item_id     = c("A","B", "C", "A","B"),
    utility     = c(60,40, 100, 10,10)
  )
  ans <- huspm(x, min_utility = 0, max_length = 2)
  ord <- order(-ans$utility, ans$length, ans$pattern)
  expect_identical(seq_len(nrow(ans)), ord)
})

test_that("top_k keeps the first K of the display order, ties included", {
  # A -> B and C are exactly tied at 100. The shorter pattern wins the tie, so
  # K = 1 must keep C, not A -> B.
  x <- data.frame(
    sequence_id = c(1,1, 2),
    item_id     = c("A","B", "C"),
    utility     = c(60,40, 100)
  )
  full <- huspm(x, min_utility = 0, max_length = 2)
  expect_identical(as.character(full$pattern), c("C", "A -> B", "A", "B"))

  k1 <- huspm(x, min_utility = 0, top_k = 1, max_length = 2)
  expect_identical(as.character(k1$pattern), "C")

  k2 <- huspm(x, min_utility = 0, top_k = 2, max_length = 2)
  expect_identical(as.character(k2$pattern), c("C", "A -> B"))
})

test_that("top_k truncation always equals head(untruncated, k)", {
  set.seed(11)
  x <- do.call(rbind, lapply(1:6, function(s) {
    n <- sample(2:6, 1)
    data.frame(sequence_id = s,
               item_id = sample(LETTERS[1:4], n, replace = TRUE),
               utility = sample(1:9, n, replace = TRUE))
  }))
  for (k in c(1L, 2L, 5L, 20L)) {
    full <- huspm(x, min_utility = 0, max_length = 4)
    tk   <- huspm(x, min_utility = 0, top_k = k, max_length = 4)
    want <- head(full, k)
    expect_identical(as.character(tk$pattern), as.character(want$pattern))
    expect_equal(tk$utility, want$utility)
    expect_identical(tk$support, want$support)
  }
})

test_that("top_k returns at most k", {
  x <- data.frame(
    sequence_id = rep(1:3, each = 4),
    item_id = rep(c("A","B","C","D"), 3),
    utility = rep(c(1,2,3,4), 3)
  )
  ans <- huspm(x, min_utility = 0, top_k = 5, max_length = 3)
  expect_lte(nrow(ans), 5)
  expect_identical(nrow(ans), 5L)
})

test_that("top_k combined with min_utility still respects both", {
  x <- data.frame(
    sequence_id = rep(1:3, each = 4),
    item_id = rep(c("A","B","C","D"), 3),
    utility = rep(c(1,2,3,4), 3)
  )
  full <- huspm(x, min_utility = 12, max_length = 3)
  tk   <- huspm(x, min_utility = 12, top_k = 3, max_length = 3)
  expect_true(all(tk$utility >= 12))
  expect_identical(as.character(tk$pattern), as.character(head(full, 3)$pattern))
})

test_that("min_support filters patterns", {
  x <- data.frame(
    sequence_id = c(1,1, 2,2, 3,3),
    item_id     = c("A","B", "A","B", "A","C"),
    utility     = c(1,1, 1,1, 1,1)
  )
  all_p <- huspm(x, min_utility = 0, max_length = 2)
  expect_true(any(all_p$pattern == "A -> C"))

  sup2 <- huspm(x, min_utility = 0, max_length = 2, min_support = 2)
  expect_true(all(sup2$support >= 2))
  expect_false(any(sup2$pattern == "A -> C"))
  expect_true(any(sup2$pattern == "A -> B"))
})

test_that("min_support above the sequence count returns nothing", {
  x <- data.frame(sequence_id = 1, item_id = "A", utility = 5)
  ans <- huspm(x, min_utility = 0, min_support = 2, max_length = 2)
  expect_identical(nrow(ans), 0L)
  expect_identical(names(ans), c("pattern","length","utility","support","utility_ratio"))
})

test_that("include_singletons controls 1-item patterns", {
  x <- data.frame(
    sequence_id = c(1,1,1),
    item_id = c("A","B","C"),
    utility = c(1, 2, 3)
  )
  with_s <- huspm(x, min_utility = 0, max_length = 2)
  without <- huspm(x, min_utility = 0, max_length = 2, include_singletons = FALSE)
  expect_true(any(with_s$length == 1))
  expect_false(any(without$length == 1))
  expect_true(any(without$pattern == "A -> B"))
})

test_that("max_length = 1 with include_singletons = FALSE is empty", {
  x <- data.frame(sequence_id = c(1,1), item_id = c("A","B"), utility = c(1,2))
  ans <- huspm(x, min_utility = 0, max_length = 1, include_singletons = FALSE)
  expect_identical(nrow(ans), 0L)
  expect_s3_class(ans, "huspm")
  expect_output(print(ans), "<huspm> 0 patterns")
})

test_that("empty results print and keep their columns", {
  x <- data.frame(sequence_id = 1, item_id = "A", utility = 5)
  ans <- huspm(x, min_utility = 999, max_length = 3)
  expect_identical(nrow(ans), 0L)
  expect_identical(names(ans), c("pattern","length","utility","support","utility_ratio"))
  expect_output(print(ans), "<huspm> 0 patterns")
})

test_that("a zero-row input returns an empty result", {
  e <- data.frame(sequence_id = integer(0), item_id = character(0),
                  utility = numeric(0))
  ans <- huspm(e, min_utility = 0)
  expect_identical(nrow(ans), 0L)
  expect_s3_class(ans, "huspm")
  expect_identical(names(ans), c("pattern","length","utility","support","utility_ratio"))
})

test_that("max_length beyond the longest sequence is harmless", {
  x <- data.frame(sequence_id = c(1,1,2), item_id = c("A","B","A"),
                  utility = c(1,2,3))
  ans <- huspm(x, min_utility = 0, max_length = 99)
  expect_true(max(ans$length) <= 2)
  expect_true(any(ans$pattern == "A -> B"))
})

test_that("a single event mines to one singleton", {
  x <- data.frame(sequence_id = 1, item_id = "A", utility = 5)
  ans <- huspm(x, min_utility = 0, max_length = 3)
  expect_identical(nrow(ans), 1L)
  expect_identical(ans$pattern, "A")
  expect_equal(ans$utility, 5)
  expect_identical(ans$support, 1L)
})

test_that("repeated items in one sequence use distinct positions", {
  x <- data.frame(
    sequence_id = c(1,1,1,1),
    item_id = c("A","A","A","A"),
    utility = c(1, 2, 3, 4)
  )
  ans <- huspm(x, min_utility = 0, max_length = 3)
  expect_equal(ans$utility[ans$pattern == "A -> A -> A"], 2 + 3 + 4)
  expect_equal(ans$utility[ans$pattern == "A"], 4)  # max single embedding, not the sum
})

test_that("numeric and factor item ids are labelled", {
  n <- data.frame(sequence_id = c(1,1,2), item_id = c(10, 20, 10), utility = c(1,2,3))
  expect_true(any(huspm(n, min_utility = 0, max_length = 2)$pattern == "10 -> 20"))

  f <- data.frame(sequence_id = c(1,1,2),
                  item_id = factor(c("x","y","x")), utility = c(1,2,3))
  expect_true(any(huspm(f, min_utility = 0, max_length = 2)$pattern == "x -> y"))
})

test_that("utility_ratio is the share of total utility", {
  x <- data.frame(
    sequence_id = c(1,1,1),
    item_id = c("A","B","C"),
    utility = c(1, 2, 3)
  )
  ans <- huspm(x, min_utility = 0, max_length = 2)
  expect_equal(ans$utility_ratio, ans$utility / sum(x$utility))
})

test_that("patterns are never duplicated", {
  set.seed(5)
  x <- do.call(rbind, lapply(1:8, function(s) {
    n <- sample(2:7, 1)
    data.frame(sequence_id = s,
               item_id = sample(LETTERS[1:5], n, replace = TRUE),
               utility = sample(1:9, n, replace = TRUE))
  }))
  ans <- huspm(x, min_utility = 0, max_length = 4)
  expect_identical(anyDuplicated(ans$pattern), 0L)
})

test_that("invalid utility is rejected", {
  x <- data.frame(
    sequence_id = 1:2,
    item_id = c("A","B"),
    utility = c(1, -1)
  )
  expect_error(huspm(x, min_utility = 1))
})

test_that("invalid parameters are rejected", {
  x <- data.frame(sequence_id = 1, item_id = "A", utility = 1)
  expect_error(huspm(x, min_utility = Inf), "finite")
  expect_error(huspm(x, min_utility = -1), "non-negative")
  expect_error(huspm(x, min_utility = 0, max_length = 0), "max_length")
  expect_error(huspm(x, min_utility = 0, min_support = 0), "min_support")
  expect_error(huspm(x, min_utility = 0, top_k = 0), "top_k")
  expect_error(huspm(x, min_utility = 0, top_k = 1.5), "top_k")
})

test_that("missing columns are reported", {
  x <- data.frame(sequence_id = 1, item_id = "A")
  expect_error(huspm(x, min_utility = 0), "Missing columns")
})
