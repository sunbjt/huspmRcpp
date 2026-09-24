#' Prepare sequence data for HUSPM
#'
#' @param data data.frame with sequence_id, item_id, utility.
#' @param sequence_id column name for sequence id.
#' @param item_id column name for item id.
#' @param utility column name for positive utility.
#' @param time optional column name used to order events within each sequence.
#' @return A list suitable for the C++ backend.
#' @export
prepare_huspm_data <- function(data,
                                sequence_id = "sequence_id",
                                item_id = "item_id",
                                utility = "utility",
                                time = NULL) {
  stopifnot(is.data.frame(data))
  needed <- c(sequence_id, item_id, utility)
  miss <- setdiff(needed, names(data))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  d <- data[, needed, drop = FALSE]
  names(d) <- c("sequence_id", "item_id", "utility")

  if (!is.null(time)) {
    if (!time %in% names(data)) stop("time column not found: ", time)
    d$..time <- data[[time]]
    d <- d[order(d$sequence_id, d$..time, seq_len(nrow(d))), , drop = FALSE]
    d$..time <- NULL
  } else {
    d <- d[order(d$sequence_id, seq_len(nrow(d))), , drop = FALSE]
  }

  if (anyNA(d$sequence_id) || anyNA(d$item_id) || anyNA(d$utility))
    stop("sequence_id, item_id and utility cannot contain NA.")
  if (!is.numeric(d$utility))
    stop("utility must be numeric.")
  if (any(!is.finite(d$utility)) || any(d$utility <= 0))
    stop("This implementation requires finite, strictly positive utilities.")

  # Keep IDs as character labels; C++ receives integer-coded IDs.
  seq_key <- as.character(d$sequence_id)
  item_key <- as.character(d$item_id)

  seq_levels <- unique(seq_key)
  item_levels <- unique(item_key)

  seq_int <- match(seq_key, seq_levels) - 1L
  item_int <- match(item_key, item_levels) - 1L

  list(
    sequence = as.integer(seq_int),
    item = as.integer(item_int),
    utility = as.numeric(d$utility),
    sequence_labels = seq_levels,
    item_labels = item_levels
  )
}

#' High-Utility Sequential Pattern Mining
#'
#' Pattern utility is defined as the sum over sequences of the maximum-utility
#' embedding of the pattern in that sequence. Utilities must be positive.
#'
#' Patterns are returned in a deterministic order: utility descending, then
#' length ascending, then pattern lexicographic ascending.
#'
#' @param data data.frame containing sequence_id, item_id, utility.
#' @param min_utility minimum total pattern utility; one finite non-negative number.
#' @param top_k if supplied, keep only the top K patterns. "Top K" means the first
#'   K patterns of the output order above, so ties are resolved by the same rule
#'   rather than by discovery order.
#' @param max_length maximum number of items in a pattern.
#' @param sequence_id sequence identifier column.
#' @param item_id item identifier column.
#' @param utility utility column.
#' @param time optional ordering column.
#' @param min_support optional minimum number of sequences containing a pattern.
#' @param include_singletons include 1-item patterns.
#' @return data.frame with pattern, length, utility, support, utility_ratio.
#' @export
huspm <- function(data,
                  min_utility,
                  top_k = NULL,
                  max_length = 5L,
                  sequence_id = "sequence_id",
                  item_id = "item_id",
                  utility = "utility",
                  time = NULL,
                  min_support = 1L,
                  include_singletons = TRUE) {
  if (!is.numeric(min_utility) || length(min_utility) != 1L ||
      !is.finite(min_utility) || min_utility < 0)
    stop("min_utility must be one finite non-negative number.")
  if (!is.null(top_k)) {
    if (length(top_k) != 1L || !is.numeric(top_k) || top_k < 1 ||
        top_k != as.integer(top_k))
      stop("top_k must be a positive integer or NULL.")
    top_k <- as.integer(top_k)
  }
  max_length <- as.integer(max_length)
  if (length(max_length) != 1L || max_length < 1L)
    stop("max_length must be >= 1.")
  min_support <- as.integer(min_support)
  if (length(min_support) != 1L || min_support < 1L)
    stop("min_support must be >= 1.")

  x <- prepare_huspm_data(data, sequence_id, item_id, utility, time)

  # An input with no events has no sequences and no items; the C++ entry point
  # rejects those counts as invalid, so short-circuit here.
  if (!length(x$sequence)) return(empty_huspm())

  out <- huspm_cpp(
    x$sequence, x$item, x$utility,
    n_sequences = length(x$sequence_labels),
    n_items = length(x$item_labels),
    min_utility = min_utility,
    # NULL would override the 0L default on the C++ int parameter and blow up with
    # "Expecting a single value: [extent=0]".
    top_k = if (is.null(top_k)) 0L else top_k,
    max_length = max_length,
    min_support = min_support,
    include_singletons = isTRUE(include_singletons)
  )

  if (!nrow(out)) return(empty_huspm())

  # Label rebuild, vectorized. At ~6e5 patterns the per-element vapply this
  # replaces cost more than the whole C++ search did.
  parts <- strsplit(out$pattern_id, ",", fixed = TRUE)
  lens  <- lengths(parts)
  labs  <- x$item_labels[as.integer(unlist(parts, use.names = FALSE)) + 1L]

  if (any(grepl("\x1f", labs, fixed = TRUE))) {
    # A label contains the join sentinel; fall back to a per-pattern join.
    grp <- rep.int(seq_along(parts), lens)
    out$pattern <- vapply(split(labs, grp), paste, character(1), collapse = " -> ")
  } else {
    # Append the separator to every item, swap it for a sentinel at each group
    # end, concatenate once, then split the result back into patterns.
    sep <- rep.int(" -> ", length(labs))
    sep[cumsum(lens)] <- "\x1f"
    out$pattern <- strsplit(paste0(labs, sep, collapse = ""), "\x1f", fixed = TRUE)[[1]]
  }

  out$pattern_id <- NULL
  out$utility_ratio <- out$utility / sum(x$utility)

  out <- out[, c("pattern", "length", "utility", "support", "utility_ratio")]
  rownames(out) <- NULL
  class(out) <- c("huspm", "data.frame")
  out
}

#' @noRd
empty_huspm <- function() {
  structure(
    data.frame(pattern = character(), length = integer(),
               utility = numeric(), support = integer(),
               utility_ratio = numeric(),
               stringsAsFactors = FALSE),
    class = c("huspm", "data.frame")
  )
}

#' Print a huspm result
#'
#' @param x a `huspm` object.
#' @param ... passed on to [print.data.frame()].
#' @return `x`, invisibly.
#' @export
print.huspm <- function(x, ...) {
  cat("<huspm> ", nrow(x), " patterns\n", sep = "")
  if (nrow(x)) print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}
