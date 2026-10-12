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

test_that("prepare_huspm_data: NA values are rejected", {
  x <- data.frame(
    sequence_id = c(1, 1, NA),
    item_id = c("A", "B", "C"),
    utility = c(1, 2, 3)
  )
  expect_error(prepare_huspm_data(x), "cannot contain NA")
})

test_that("prepare_huspm_data: non-numeric utility is rejected", {
  x <- data.frame(
    sequence_id = c(1, 1),
    item_id = c("A", "B"),
    utility = c("1", "2")
  )
  expect_error(prepare_huspm_data(x), "must be numeric")
})

test_that("prepare_huspm_data: non-finite or non-positive utility is rejected", {
  x1 <- data.frame(
    sequence_id = c(1, 1),
    item_id = c("A", "B"),
    utility = c(1, Inf)
  )
  expect_error(prepare_huspm_data(x1), "strictly positive utilities")

  x2 <- data.frame(
    sequence_id = c(1, 1),
    item_id = c("A", "B"),
    utility = c(1, 0)
  )
  expect_error(prepare_huspm_data(x2), "strictly positive utilities")

  x3 <- data.frame(
    sequence_id = c(1, 1),
    item_id = c("A", "B"),
    utility = c(1, -1)
  )
  expect_error(prepare_huspm_data(x3), "strictly positive utilities")
})

test_that("prepare_huspm_data: output structure is correct", {
  x <- data.frame(
    sequence_id = c(1,1,1, 2,2, 2, 3,3),
    item_id     = c("A","B","C", "A","B","D", "A","C"),
    utility     = c(5, 2, 10, 4, 3, 8, 6, 7)
  )
  pd <- prepare_huspm_data(x)
  expect_type(pd, "list")
  expect_named(pd, c("sequence", "item", "utility", "sequence_labels", "item_labels"))
  expect_equal(length(pd$sequence), nrow(x))
  expect_equal(length(pd$item), nrow(x))
  expect_equal(length(pd$utility), nrow(x))
  expect_equal(length(pd$sequence_labels), 3L)
  expect_equal(length(pd$item_labels), 4L)
})

test_that("prepare_huspm_data: time ordering is respected", {
  x <- data.frame(
    sequence_id = c(1,1,1),
    item_id = c("B","A","C"),
    utility = c(1, 5, 10),
    time = c(3, 1, 2)
  )
  pd <- prepare_huspm_data(x, time = "time")
  expect_equal(pd$item_labels[pd$item + 1L], c("A", "C", "B"))
})

test_that("prepare_huspm_data: ties on time keep input order", {
  x <- data.frame(
    sequence_id = c(1,1,1),
    item_id = c("B","A","C"),
    utility = c(1, 5, 10),
    time = c(1, 1, 1)
  )
  pd <- prepare_huspm_data(x, time = "time")
  expect_equal(pd$item_labels[pd$item + 1L], c("B", "A", "C"))
})

test_that("label rebuilding fallback works for labels containing unit separator", {
  # Test fallback when one label contains \x1f
  x <- data.frame(
    sequence_id = c(1,1, 2,2),
    item_id = c(paste0("A", "\x1f", "B"), "C", paste0("A", "\x1f", "B"), "D"),
    utility = c(5, 10, 4, 8)
  )
  ans <- huspm(x, min_utility = 0, max_length = 2)
  expect_s3_class(ans, "huspm")
  expect_true(any(grepl(paste0("A", "\x1f", "B"), ans$pattern)))
  expect_true(any(grepl("A", ans$pattern)))
  # Check that the combined pattern is correctly formatted
  expect_true(any(ans$pattern == paste0("A", "\x1f", "B", " -> ", "C")))
})

test_that("item labels with -> in them work correctly", {
  # This is fine in the vectorized path, but test it anyway
  x <- data.frame(
    sequence_id = c(1,1, 2,2),
    item_id = c("A -> B", "C", "A -> B", "D"),
    utility = c(5, 10, 4, 8)
  )
  ans <- huspm(x, min_utility = 0, max_length = 2)
  expect_s3_class(ans, "huspm")
  expect_true(any(grepl("A -> B", ans$pattern)))
  expect_true(any(ans$pattern == "A -> B -> C"))
})

# Brute-force oracle for correctness verification
# Verified to match HUSPM semantics: for each pattern, sum the *maximum* utility
# embedding per sequence, and support is the number of sequences that have at least
# one embedding.
bruteforce_huspm <- function(data, min_utility, max_length, min_support, include_singletons) {
  # Prepare data
  x <- prepare_huspm_data(data)
  n_seqs <- length(x$sequence_labels)

  # Extract per-sequence items and utilities (item ids are 0-based matching x$item)
  seqs <- list()
  for (s in seq_len(n_seqs)) {
    idx <- which(x$sequence == (s - 1L))
    seqs[[s]] <- data.frame(
      item = x$item[idx],
      util = x$utility[idx],
      stringsAsFactors = FALSE
    )
  }

  # List to collect total results
  # each entry: list(total_utility, support)
  # name: comma-separated 0-based item ids
  total_results <- list()

  # Recursive search per sequence
  # For each sequence, generate all possible patterns and keep max utility
  # Then accumulate across all sequences
  for (s in seq_len(n_seqs)) {
    seq <- seqs[[s]]
    n <- nrow(seq)
    # Track maximum utility for each pattern in this sequence - use environment for pass-by-reference
    seq_max <- new.env(hash = TRUE)

    recurse <- function(last_idx, pattern, current_util) {
      len <- length(pattern)
      if (len > 0) {
        key <- paste(pattern, collapse = ",")
        if (!exists(key, envir = seq_max) || current_util > get(key, envir = seq_max)) {
          assign(key, current_util, envir = seq_max)
        }
      }

      if (len >= max_length) return()

      # Extend with any item after last_idx: last_idx is 1-based position in this sequence
      start <- if (len == 0) 1 else last_idx + 1L
      if (start <= n) {
        for (i in start:n) {
          # i is position, item value is seq$item[i]
          new_util <- current_util + seq$util[i]
          new_pattern <- c(pattern, seq$item[i])
          recurse(i, new_pattern, new_util)
        }
      }
    }

    recurse(0, integer(0), 0.0)

    # After processing the entire sequence, add the max for each pattern to the total
    for (key in ls(seq_max)) {
      max_util <- get(key, envir = seq_max)
      if (!(key %in% names(total_results))) {
        total_results[[key]] <- list(total_utility = max_util, support = 1L)
      } else {
        total_results[[key]]$total_utility <- total_results[[key]]$total_utility + max_util
        total_results[[key]]$support <- total_results[[key]]$support + 1L
      }
    }
  }

  # Build result data frame
  keys <- names(total_results)
  if (!length(keys)) {
    return(huspmRcpp:::empty_huspm())
  }

  res_df <- do.call(rbind, lapply(keys, function(key) {
    items <- as.integer(strsplit(key, ",")[[1]])
    pattern <- paste(x$item_labels[items + 1L], collapse = " -> ")
    total_util <- total_results[[key]]$total_utility
    support <- total_results[[key]]$support
    data.frame(
      pattern = pattern,
      length = length(items),
      utility = total_util,
      support = support,
      stringsAsFactors = FALSE
    )
  }))

  # Filter based on parameters
  if (!include_singletons) {
    res_df <- res_df[res_df$length > 1, , drop = FALSE]
  }
  res_df <- res_df[res_df$utility >= min_utility & res_df$support >= min_support, , drop = FALSE]

  if (nrow(res_df) == 0) {
    return(empty_huspm())
  }

  # Compute utility ratio and sort
  res_df$utility_ratio <- res_df$utility / sum(x$utility)
  # Sort order matches main algorithm: utility descending, length ascending, pattern ascending
  res_df <- res_df[order(-res_df$utility, res_df$length, res_df$pattern),
                   c("pattern", "length", "utility", "support", "utility_ratio")]
  rownames(res_df) <- NULL
  class(res_df) <- c("huspm", "data.frame")
  return(res_df)
}

test_that("results match brute-force oracle on small random data", {
  set.seed(42)
  for (i in 1:10) {
    # Small dataset so brute force can enumerate everything
    n_seqs <- sample(2:4, 1)
    n_items <- sample(3:5, 1)
    max_len <- sample(2:3, 1)
    x <- do.call(rbind, lapply(1:n_seqs, function(s) {
      n <- sample(2:5, 1)
      data.frame(
        sequence_id = s,
        item_id = sample(LETTERS[1:n_items], n, replace = TRUE),
        utility = sample(1:5, n, replace = TRUE)
      )
    }))
    huspm_res <- huspm(x, min_utility = 0, max_length = max_len, min_support = 1)
    bf_res <- bruteforce_huspm(x, min_utility = 0, max_length = max_len,
                               min_support = 1, include_singletons = TRUE)

    expect_equal(nrow(huspm_res), nrow(bf_res))
    expect_equal(huspm_res$pattern, bf_res$pattern)
    expect_equal(huspm_res$utility, bf_res$utility, tolerance = 1e-8)
    expect_equal(huspm_res$support, bf_res$support)
  }
})

test_that("results match brute-force with min_utility and min_support filtering", {
  set.seed(123)
  x <- data.frame(
    sequence_id = c(1,1,1, 2,2, 2, 3,3),
    item_id     = c("A","B","C", "A","B","D", "A","C"),
    utility     = c(5, 2, 10, 4, 3, 8, 6, 7)
  )
  huspm_res <- huspm(x, min_utility = 15, max_length = 3, min_support = 1)
  bf_res <- bruteforce_huspm(x, min_utility = 15, max_length = 3,
                             min_support = 1, include_singletons = TRUE)
  expect_equal(nrow(huspm_res), nrow(bf_res))
  expect_equal(huspm_res$pattern, bf_res$pattern)
  expect_equal(huspm_res$utility, bf_res$utility)
})
