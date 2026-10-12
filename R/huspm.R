#' 预处理高效用频繁序列挖掘的输入数据
#'
#' 将输入数据框转换为算法需要的内部格式。输入数据框必须包含序列ID、项目ID和项目效用。
#'
#' @param data 输入数据框，包含 `sequence_id`, `item_id`, `utility` 列。
#' @param sequence_id 表示序列ID的列名。
#' @param item_id 表示项目ID的列名。
#' @param utility 表示项目效用（权重）的列名。
#' @param time 可选列名，用于排序同一个序列中的事件顺序。如果不提供，则按输入数据的顺序排序。
#' @return 包含以下组件的列表：
#'   \item{sequence}{0-based整数序列ID}
#'   \item{item}{0-based整数项目ID}
#'   \item{utility}{数值型项目效用}
#'   \item{sequence_labels}{原始序列标签字符向量}
#'   \item{item_labels}{原始项目标签字符向量}
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

#' 高效用频繁序列模式挖掘 (High-Utility Sequential Pattern Mining)
#'
#' 挖掘数据集中所有满足最小效用阈值的序列模式。模式效用定义为每个序列中该模式的最大嵌入效用之和。
#' 所有输入项目的效用必须为正。
#'
#' 挖掘结果按以下规则确定顺序：先按模式效用降序排列，再按模式长度升序排列，最后按模式字典序升序排列。
#'
#' 当有多个根项目时，算法会自动利用多线程并行处理每个根项目，自动利用多个CPU核心加速计算。
#'
#' @param data 输入数据框，包含 `sequence_id`, `item_id`, `utility` 列。
#' @param min_utility 最小总模式效用，必须是一个有限非负数。所有满足条件的模式都会被输出。
#' @param top_k 如果提供，则只保留前 K 个满足条件的模式。这里的"前 K"指的是按输出顺序保留前 K 个，
#'   因此平局会按照相同规则解决，而不是按照发现顺序。如果不提供，会输出所有满足条件的模式。
#' @param max_length 模式中允许的最大项目个数，必须是正整数。
#' @param sequence_id 表示序列ID的列名。
#' @param item_id 表示项目ID的列名。
#' @param utility 表示项目效用（权重）的列名。
#' @param time 可选列名，用于排序同一个序列中的事件顺序。如果不提供，则按输入数据的顺序排序。
#' @param min_support 可选参数，表示一个模式至少需要出现在多少个不同序列中才能被保留。必须是正整数。
#' @param include_singletons 是否输出长度为1的单项目模式。默认为 TRUE。
#' @return 数据框，包含以下列：
#'   \item{pattern}{字符型模式标签，项目之间用 ` -> ` 分隔}
#'   \item{length}{整型，模式中包含的项目个数}
#'   \item{utility}{数值型，模式的总效用（每个序列中最大嵌入效用之和）}
#'   \item{support}{整型，包含这个模式的序列个数}
#'   \item{utility_ratio}{数值型，模式总效用除以所有输入事件效用之和。注意：由于模式效用是每个序列中的最大嵌入，
#'     分母是所有输入事件的总和，因此这个比例不会加起来等于 1.0，并且可能会对出现在多个模式中的事件重复计数。}
#' @export
huspm <- function(data,
                  min_utility = 0,
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
  if (!is.numeric(max_length) || length(max_length) != 1L || max_length != as.integer(max_length) || max_length < 1L) {
    stop("max_length must be a positive integer.")
  }
  max_length <- as.integer(max_length)
  if (!is.numeric(min_support) || length(min_support) != 1L || min_support != as.integer(min_support) || min_support < 1L) {
    stop("min_support must be a positive integer.")
  }
  min_support <- as.integer(min_support)
  if (!is.logical(include_singletons) || length(include_singletons) != 1L || is.na(include_singletons)) {
    stop("include_singletons must be TRUE or FALSE.")
  }

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
    include_singletons = include_singletons
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
  # Re-sort to match documented ordering: utility descending, length ascending, pattern ascending
  out <- out[order(-out$utility, out$length, out$pattern), ]
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

#' 打印 huspm 挖掘结果
#'
#' @param x `huspm` 对象。
#' @param max_rows 最多打印多少个模式。默认为 20，使用 `max_rows = Inf` 打印所有结果。
#' @param ... 传递给 `print.data.frame()` 的额外参数。
#' @return 原对象 `x`，不可见返回。
#' @importFrom utils head
#' @export
print.huspm <- function(x, max_rows = 20, ...) {
  cat("<huspm> ", nrow(x), " patterns\n", sep = "")
  if (nrow(x) > 0) {
    if (!is.infinite(max_rows) && nrow(x) > max_rows) {
      x_head <- head(x, max_rows)
      print.data.frame(x_head, row.names = FALSE, ...)
      cat("... (", nrow(x) - max_rows, " more rows omitted, use max_rows = Inf to see all)\n", sep = "")
    } else {
      print.data.frame(x, row.names = FALSE, ...)
    }
  }
  invisible(x)
}
