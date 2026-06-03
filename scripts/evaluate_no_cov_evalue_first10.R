base_dir <- "/Users/eilpis/Desktop/projects/AMRfinder-working"
input_root <- file.path(base_dir, "simulation_amrfinder_input")

runs_env <- Sys.getenv("RUNS", unset = "1,2,3,4,5,6,7,8,9,10")
runs <- as.integer(strsplit(runs_env, ",", fixed = TRUE)[[1]])
evalue_cutoff <- as.numeric(Sys.getenv("EVALUE_CUTOFF", unset = "20"))
score_column <- Sys.getenv("EVALUE_SCORE_COLUMN", unset = "e_adjust")
e_bh_alpha <- as.numeric(Sys.getenv("E_BH_ALPHA", unset = "0.0005"))
cutoff_label <- gsub("[^0-9A-Za-z]+", "_", format(evalue_cutoff, trim = TRUE, scientific = FALSE))
score_label <- gsub("[^0-9A-Za-z]+", "_", score_column)
alpha_label <- gsub("[^0-9A-Za-z]+", "_", format(e_bh_alpha, trim = TRUE, scientific = FALSE))

calc_evalue_for_region <- function(region_dat, y_group) {
  sample_mean <- colMeans(region_dat, na.rm = TRUE)
  control_x <- sample_mean[y_group == 0]
  test_x <- sample_mean[y_group == 1]
  control_x <- control_x[is.finite(control_x)]
  test_x <- test_x[is.finite(test_x)]

  if (length(control_x) < 2 || length(test_x) < 2) {
    return(1)
  }

  group_var <- function(x) mean((x - mean(x))^2)
  control_var <- group_var(control_x)
  test_var <- group_var(test_x)
  pooled_x <- c(control_x, test_x)
  pooled_var <- group_var(pooled_x)

  if (!is.finite(pooled_var) || pooled_var <= 0) {
    return(1)
  }
  if (!is.finite(control_var) || control_var <= 0 || !is.finite(test_var) || test_var <= 0) {
    return(.Machine$double.xmax)
  }

  n_control <- length(control_x)
  n_test <- length(test_x)
  log_e_value <- ((n_control + n_test) / 2) * log(pooled_var) -
    (n_control / 2) * log(control_var) -
    (n_test / 2) * log(test_var)

  if (!is.finite(log_e_value) || log_e_value <= 0) {
    return(1)
  }
  if (log_e_value >= log(.Machine$double.xmax)) {
    return(.Machine$double.xmax)
  }

  exp(log_e_value)
}

adjust_evalue_bh <- function(e_value) {
  p_value <- rep(1, length(e_value))
  finite_id <- is.finite(e_value) & e_value > 0
  p_value[finite_id] <- pmin(1, 1 / e_value[finite_id])
  p_value[is.infinite(e_value) & e_value > 0] <- 0

  adjusted_p <- p.adjust(p_value, method = "BH")
  e_adjust <- rep(1, length(adjusted_p))
  zero_id <- adjusted_p == 0
  positive_id <- adjusted_p > 0
  e_adjust[zero_id] <- Inf
  e_adjust[positive_id] <- 1 / adjusted_p[positive_id]
  e_adjust[!is.finite(e_adjust) & !zero_id] <- 1
  e_adjust
}

e_bh_significant <- function(e_value, alpha = 0.05) {
  significant <- integer(length(e_value))
  valid_id <- which(is.finite(e_value) & e_value > 0)
  if (length(valid_id) == 0) {
    return(significant)
  }

  ordered_id <- valid_id[order(e_value[valid_id], decreasing = TRUE)]
  e_sorted <- e_value[ordered_id]
  k_sequence <- seq_along(e_sorted)
  k_total <- length(e_value)
  valid_k <- which((k_sequence * e_sorted / k_total) >= (1 / alpha))
  if (length(valid_k) == 0) {
    return(significant)
  }

  k_star <- max(valid_k)
  significant[ordered_id[seq_len(k_star)]] <- 1L
  significant
}

add_evalue <- function(pred, dat, y_group) {
  if (!"e_value" %in% names(pred)) {
    pred$e_value <- vapply(seq_len(nrow(pred)), function(i) {
      idx <- which(dat$chr == pred$chr[i] & dat$pos > pred$start[i] & dat$pos <= pred$end[i])
      if (length(idx) == 0) {
        return(1)
      }
      calc_evalue_for_region(dat[idx, -c(1, 2), drop = FALSE], y_group)
    }, numeric(1))
  }

  if (!"e_adjust" %in% names(pred)) {
    pred$e_adjust <- adjust_evalue_bh(pred$e_value)
  }
  if (!"e_bh_significant" %in% names(pred)) {
    pred$e_bh_significant <- e_bh_significant(pred$e_value, alpha = e_bh_alpha)
  }

  pred
}

calc_auc <- function(labels, scores) {
  labels <- as.integer(labels)
  ok <- is.finite(scores) & labels %in% c(0L, 1L)
  labels <- labels[ok]
  scores <- scores[ok]
  n_pos <- sum(labels == 1L)
  n_neg <- sum(labels == 0L)
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

evaluate_evalue <- function(run) {
  run_dir <- file.path(input_root, as.character(run))
  results_dir <- file.path(run_dir, "amrfinder_results")
  input_file <- file.path(run_dir, "amrfinder_common_site_input.tsv")
  truth_file <- file.path(run_dir, "truth", "DMRs_unDMRs_signal.bed")
  result_file <- file.path(results_dir, sprintf("run%d_dmr.no.cov.tsv", run))
  evalue_file <- file.path(results_dir, sprintf("run%d_dmr.no.cov.evalue.tsv", run))
  prefix <- if (score_column == "e_bh_significant") {
    sprintf("%s_alpha_%s", score_label, alpha_label)
  } else {
    score_label
  }
  filtered_file <- file.path(results_dir, sprintf("run%d_dmr.no.cov.%s_gt_%s.tsv", run, prefix, cutoff_label))
  coverage_file <- file.path(results_dir, sprintf("run%d_dmr.no.cov.%s_gt_%s.testR_truth_coverage.tsv", run, prefix, cutoff_label))
  metrics_file <- file.path(results_dir, sprintf("run%d_dmr.no.cov.testR_metrics_%s_gt_%s.tsv", run, prefix, cutoff_label))

  if (!file.exists(result_file)) {
    stop("Missing no.cov result file for run ", run, ": ", result_file)
  }

  pred_all <- read.delim(result_file, check.names = FALSE)
  dat <- read.delim(input_file, check.names = FALSE)
  y_group <- c(rep(0, 8), rep(1, 8))

  pred_all <- add_evalue(pred_all, dat, y_group)
  if (!score_column %in% names(pred_all)) {
    stop("Missing score column: ", score_column)
  }
  write.table(pred_all, evalue_file, sep = "\t", quote = FALSE, row.names = FALSE)

  score_value <- pred_all[[score_column]]
  pred <- pred_all[!is.na(score_value) & score_value > evalue_cutoff, , drop = FALSE]
  write.table(pred, filtered_file, sep = "\t", quote = FALSE, row.names = FALSE)

  truth <- read.table(truth_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(truth) <- c("chr", "start", "end", "isDMR")
  truth_eval <- truth
  truth_eval$isFind <- 0L
  truth_eval$covrate <- 0

  for (i in seq_len(nrow(truth_eval))) {
    b_chr <- truth_eval$chr[i]
    b_start <- truth_eval$start[i]
    b_end <- truth_eval$end[i]
    hits <- pred[pred$chr == b_chr & pred$start < b_end & pred$end > b_start, , drop = FALSE]
    if (nrow(hits) > 0) {
      covlen <- pmax(0, pmin(hits$end, b_end) - pmax(hits$start, b_start))
      q_temp <- sum(covlen)
      b_len <- b_end - b_start
      truth_eval$covrate[i] <- q_temp / b_len
      if (q_temp > b_len / 2) {
        truth_eval$isFind[i] <- 1L
      }
    }
  }

  tp <- sum(truth_eval$isDMR == 1 & truth_eval$isFind == 1)
  fn <- sum(truth_eval$isDMR == 1 & truth_eval$isFind == 0)
  fp <- sum(truth_eval$isDMR == 0 & truth_eval$isFind == 1)
  tn <- sum(truth_eval$isDMR == 0 & truth_eval$isFind == 0)

  metrics <- data.frame(
    run = run,
    method = paste0("dmr.no.cov.", score_column),
    cutoff = evalue_cutoff,
    score_column = score_column,
    e_bh_alpha = e_bh_alpha,
    ACC = (tp + tn) / (tp + fn + fp + tn),
    FDR = ifelse((fp + tp) == 0, NA_real_, fp / (fp + tp)),
    Type_I_error = ifelse((fp + tn) == 0, NA_real_, fp / (fp + tn)),
    Power = ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn)),
    AUC = calc_auc(truth_eval$isDMR, truth_eval$covrate),
    TP = tp,
    FP = fp,
    TN = tn,
    FN = fn,
    total_regions = nrow(pred_all),
    significant_regions = nrow(pred),
    result_file = evalue_file,
    filtered_file = filtered_file,
    coverage_file = coverage_file,
    stringsAsFactors = FALSE
  )

  write.table(truth_eval, coverage_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(metrics, metrics_file, sep = "\t", quote = FALSE, row.names = FALSE)
  metrics
}

summary_dir <- file.path(input_root, "batch_results")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
summary_prefix <- if (score_column == "e_bh_significant") {
  sprintf("%s_alpha_%s", score_label, alpha_label)
} else {
  score_label
}
summary_file <- file.path(summary_dir, sprintf("run10_no_cov_%s_gt_%s_metrics_summary.tsv", summary_prefix, cutoff_label))

all_metrics <- do.call(rbind, lapply(runs, evaluate_evalue))
write.table(all_metrics, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
print(all_metrics[, c("run", "method", "ACC", "FDR", "Type_I_error", "Power", "AUC", "significant_regions")])
message("Wrote summary: ", summary_file)
