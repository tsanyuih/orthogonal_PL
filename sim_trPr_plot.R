## Summarize saved treatment-allocation sensitivity results from trPr.rds.
## Exports trPr_n<n>.pdf with panels for J=3,5,8 and returns summaries in memory.
## Monte Carlo standard errors quantify uncertainty in simulation means.
## Run after sim_trPr.R, using identical profile and sensitivity options.

suppressPackageStartupMessages(library(ggplot2))
trpr_plot_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(files)) tail(files, 1)[[1]] else
    if (length(arg)) sub("^--file=", "", arg[[1]]) else
      stop("Run with source('/path/to/sim_trPr_plot.R') or Rscript.")
})
local({
  saved_options <- options(opl.trPr.functions_only = TRUE)
  on.exit(options(saved_options))
  source(file.path(dirname(normalizePath(trpr_plot_file, winslash = "/")),
                   "sim_trPr.R"), local = .GlobalEnv)
})
estimator_labels <- c(plugin = "Plug-in", direct_if = "Direct IF",
  smoothed_plugin = "Smoothed plug-in", orthogonal_smoothed = "Orthogonal smoothed")


validate_trPr <- function(results) {
  params <- expand.grid(Cu = Cu_grid, rate = sensitivity_rate, n = n, J = J)
  if (!identical(attr(results, "run_config"), run_config) ||
      !valid_result(results, nrow(params) * R))
    stop("Results are incomplete or have a different configuration.")
  for (t in seq_len(nrow(params))) {
    part <- results[results$n == params$n[t] & results$J == params$J[t] &
                    results$Cu == params$Cu[t] & results$rate == params$rate[t], ]
    if (!valid_scenario(part, R, params$n[t], params$J[t],
                        params$rate[t], params$Cu[t]))
      stop("Missing/duplicate repetitions in scenario ", t)
  }
  probability_columns <- grep("\\.(treat_rate|misclass_rate)$", names(results), value = TRUE)
  if (any(vapply(results[probability_columns],
                 function(x) any(x < 0 | x > 1), logical(1))))
    stop("Treatment/misclassification rates must lie in [0,1].")
  invisible(TRUE)
}

build_trPr_summary <- function(results) {
  validate_trPr(results)
  groups <- split(seq_len(nrow(results)),
    interaction(results$J, results$n, results$rate, results$Cu, drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(index) {
    part <- results[index, , drop = FALSE]
    do.call(rbind, lapply(estimator_names, function(estimator) {
      value <- function(metric) part[[paste(estimator, metric, sep = ".")]]
      mcse <- function(x) sd(x) / sqrt(length(x))
      error <- value("Rsup_estimate") - value("Rsup_truth")
      data.frame(J = part$J[1], n = part$n[1], rate = part$rate[1],
        Cu = part$Cu[1], estimator = estimator, R = nrow(part),
        treat_rate = mean(value("treat_rate")),
        treat_rate_mcse = mcse(value("treat_rate")),
        Rsup_estimate = mean(value("Rsup_estimate")),
        Rsup_truth = mean(value("Rsup_truth")),
        signed_bias = mean(error), rmse = sqrt(mean(error^2)),
        excess_regret = mean(value("excess")), excess_mcse = mcse(value("excess")),
        misclass_rate = mean(value("misclass_rate")),
        misclass_rate_mcse = mcse(value("misclass_rate")))
    }))
  }))
  rownames(out) <- NULL
  out[order(out$J, out$n, out$rate, out$Cu,
            match(out$estimator, estimator_names)), ]
}

plot_trPr <- function(summary, n_current) {
  dat <- summary[summary$n == n_current, ]
  dat$J_panel <- factor(dat$J, levels = c(3, 5, 8),
                        labels = c("J = 3", "J = 5", "J = 8"))
  dat$estimator <- factor(dat$estimator, estimator_names)
  line_dat <- dat[ave(dat$Cu, dat$J_panel, dat$estimator, FUN = length) > 1L, ]
  plot <- ggplot(dat, aes(Cu, treat_rate, color = estimator, linetype = estimator,
                         shape = estimator)) +
    geom_line(data = line_dat) + geom_point(size = 2.15, stroke = 0.35) +
    facet_grid(. ~ J_panel, drop = FALSE) +
    scale_color_manual(values = c(plugin = "#0072B2", direct_if = "#8E44AD",
                                  smoothed_plugin = "#009E73", orthogonal_smoothed = "#D62728"),
                       labels = estimator_labels, drop = FALSE) +
    scale_linetype_manual(values = c(plugin = "dashed", direct_if = "dotdash",
                                     smoothed_plugin = "dotted", orthogonal_smoothed = "solid"),
                          labels = estimator_labels, drop = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    scale_shape_manual(values = c(plugin = 16, direct_if = 18,
                                  smoothed_plugin = 15, orthogonal_smoothed = 17),
                       labels = estimator_labels, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
    labs(x = expression(C[u]), y = "Treated proportion", color = NULL,
         linetype = NULL, shape = NULL) +
    guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1),
           shape = guide_legend(ncol = 1)) +
    theme_bw(base_size = 13) + theme(
      text = element_text(size = 13),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 13),
      strip.text = element_text(size = 13),
      plot.title = element_text(size = 13),
      plot.subtitle = element_text(size = 13),
      plot.caption = element_text(size = 13),
      plot.tag = element_text(size = 13),
      legend.position = "inside",
      legend.position.inside = c(1, 1),
      legend.justification.inside = c(1, 1),
      legend.direction = "vertical",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = grid::unit(10, "pt"),
      legend.key.spacing.y = grid::unit(1, "pt"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(4, 4, 4, 4),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10))
  ## Anchor to the actual first panel cell, not a fraction of the total width:
  ## facet spacing and axis widths would make that fractional position inexact.
  grob <- ggplotGrob(plot)
  legend_panel <- which(grob$layout$name == "panel-1-1")
  legend <- which(grob$layout$name == "guide-box-inside")
  stopifnot(length(legend_panel) == 1L, length(legend) == 1L)
  grob$layout[legend, c("t", "l", "b", "r")] <-
    grob$layout[legend_panel, c("t", "l", "b", "r")]
  grob
}

run_trPr_plot <- function() {
  ## USER PATH: replace the empty string with your absolute figure directory.
  ## Input results use dir_out in sim_trPr.R. Use forward slashes on Windows.
  figure_dir <- require_user_path(getOption("opl.trPr.figure_dir", ""),
    "figure_dir in sim_trPr_plot.R")
  file <- result_path("trPr.rds")
  recover_rds_backup(file)
  if (!file.exists(file)) stop("Run sim_trPr.R first: ", file)
  summary <- build_trPr_summary(readRDS(file))
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  figure_path <- function(name) file.path(figure_dir, paste0(output_prefix, name))
  for (n_current in sort(unique(summary$n))) {
    plot <- plot_trPr(summary, n_current)
    output <- figure_path(sprintf("trPr_n%d.pdf", n_current))

    ggsave(output, plot, width = 10, height = 3.6, units = "in")
    message("Saved: ", output)
  }
  invisible(summary)
}
if (!isTRUE(getOption("opl.trPr.plot_functions_only", FALSE))) run_trPr_plot()


