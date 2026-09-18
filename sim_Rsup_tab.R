## Publication tables and figures from the existing Rsup.rds.
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

## Locate this script to load shared helpers with source() or Rscript.
script_file <- local({
  source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(source_files)) tail(source_files, 1)[[1]] else
    if (length(arg)) sub("^--file=", "", arg[[1]]) else
      stop("Run with source('/path/to/this_script.R') or Rscript.")
})
previous_source_options <- options(opl.Rsup.functions_only = TRUE)
tryCatch(
  source(file.path(dirname(normalizePath(script_file, winslash = "/")), "data_gen.R")),
  finally = options(previous_source_options))

estimator_labels <- c(
  plugin = "Plug-in",
  direct_if = "Direct IF",
  smoothed_plugin = "Smoothed plug-in",
  orthogonal_smoothed = "Orthogonal smoothed"
)
display_scale <- 100
## USER PATHS: input inherits dir_out from data_gen.R.
## Replace the empty string below with your absolute table/figure directory.
## Use forward slashes on Windows. The R option may override this path.
rsup_data_dir <- dir_out
publication_dir <- path.expand(getOption("opl.Rsup.figure_dir", ""))
publication_path <- function(name) {
  require_user_path(publication_dir, "publication_dir in sim_Rsup_tab.R")
  dir.create(publication_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(publication_dir, paste0(output_prefix, name))
}

build_Rsup_summary <- function(results) {
  required <- c("n", "J", "rate", "rep",
    unlist(lapply(names(estimator_labels), function(estimator) {
      paste0(estimator, c(".Rsup_estimate", ".Rsup_truth", ".excess", ".misclass_rate"))
    }), use.names = FALSE))
  stopifnot(is.data.frame(results), all(required %in% names(results)))

  long <- do.call(rbind, lapply(names(estimator_labels), function(estimator) {
    estimate <- results[[paste0(estimator, ".Rsup_estimate")]]
    truth <- results[[paste0(estimator, ".Rsup_truth")]]
    data.frame(
      n = results$n,
      J = results$J,
      r = results$rate,
      rep = results$rep,
      estimator = unname(estimator_labels[[estimator]]),
      error = estimate - truth,
      excess = results[[paste0(estimator, ".excess")]],
      misclass_rate = results[[paste0(estimator, ".misclass_rate")]]
    )
  }))
  numeric_columns <- long[vapply(long, is.numeric, logical(1))]
  if (anyNA(long) || !all(vapply(numeric_columns, function(x) all(is.finite(x)), logical(1))))
    stop("Missing or nonfinite results; failed replications are not silently dropped.")

  long %>%
    group_by(J, n, r, estimator) %>%
    summarise(
      R = n(),
      signed_bias = mean(error),
      absolute_bias = abs(mean(error)),
      rmse = sqrt(mean(error^2)),
      excess_regret = mean(excess),
      excess_mcse = sd(excess) / sqrt(n()),
      misclassification_rate = mean(misclass_rate),
      .groups = "drop"
    ) %>%
    mutate(estimator = factor(estimator, levels = unname(estimator_labels))) %>%
    arrange(J, n, r, estimator) %>%
    mutate(estimator = as.character(estimator))
}

write_Rsup_excel <- function(summary, file) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Excel export requires the R package 'openxlsx'. Install it with ",
         "install.packages('openxlsx').")
  }

  estimators <- unname(estimator_labels)
  required <- c("J", "n", "r", "estimator", "absolute_bias", "rmse",
                "excess_regret")
  if (!is.data.frame(summary) || !all(required %in% names(summary)))
    stop("The summary does not contain all columns required for Excel export.")

  export <- summary[, required]
  J.values <- sort(unique(export$J))
  n.values <- sort(unique(export$n))
  r.values <- sort(unique(export$r))
  if (length(J.values) != 3L ||
      any(!estimators %in% unique(export$estimator)))
    stop("Expected three J values and all four requested estimators.")

  clean_rounded <- function(value, digits = 3L) {
    rounded <- round(value * display_scale, digits)
    rounded[rounded == 0] <- 0
    rounded
  }

  workbook <- openxlsx::createWorkbook(creator = "R")
  font.family <- "Arial"
  title.style <- openxlsx::createStyle(
    fontName = font.family, fontSize = 15, textDecoration = "bold",
    fontColour = "#111827"
  )
  subtitle.style <- openxlsx::createStyle(
    fontName = font.family, fontSize = 10, textDecoration = "italic",
    fontColour = "#4B5563"
  )
  header.style <- openxlsx::createStyle(
    fontName = font.family, fontSize = 10, textDecoration = "bold",
    fontColour = "#FFFFFF", fgFill = "#1F4E78",
    halign = "center", valign = "center", wrapText = TRUE,
    border = c("top", "bottom", "left", "right"),
    borderStyle = "thin", borderColour = "#FFFFFF"
  )
  body.style <- openxlsx::createStyle(
    fontName = font.family, fontSize = 10, fontColour = "#111827",
    halign = "center", valign = "center"
  )
  note.style <- openxlsx::createStyle(
    fontName = font.family, fontSize = 9, fontColour = "#4B5563"
  )
  integer.style <- openxlsx::createStyle(numFmt = "0")
  rate.style <- openxlsx::createStyle(numFmt = "0.00")
  result.style <- openxlsx::createStyle(numFmt = "0.000")
  separator.style <- openxlsx::createStyle(
    border = "top", borderStyle = "medium", borderColour = "#7F8C8D"
  )
  bottom.style <- openxlsx::createStyle(
    border = "bottom", borderStyle = "medium", borderColour = "#7F8C8D"
  )
  alternate.style <- openxlsx::createStyle(fgFill = "#F6F8FB")

  for (J.current in J.values) {
    sheet <- paste0("J", J.current)
    openxlsx::addWorksheet(workbook, sheet, gridLines = FALSE)
    openxlsx::mergeCells(workbook, sheet, cols = 1:14, rows = 1)
    openxlsx::writeData(
      workbook, sheet, paste0("Rsup simulation results: J = ", J.current),
      startCol = 1, startRow = 1, colNames = FALSE
    )
    openxlsx::mergeCells(workbook, sheet, cols = 1:14, rows = 2)
    openxlsx::writeData(
      workbook, sheet,
      sprintf("%d repetitions; all reported values are multiplied by 100.", R),
      startCol = 1, startRow = 2, colNames = FALSE
    )

    openxlsx::mergeCells(workbook, sheet, cols = 1, rows = 4:5)
    openxlsx::mergeCells(workbook, sheet, cols = 2, rows = 4:5)
    openxlsx::writeData(workbook, sheet, "n", 1, 4, colNames = FALSE)
    openxlsx::writeData(workbook, sheet, "r", 2, 4, colNames = FALSE)
    for (estimator_index in seq_along(estimators)) {
      start.col <- 3L + (estimator_index - 1L) * 3L
      openxlsx::mergeCells(workbook, sheet, cols = start.col:(start.col + 2L),
                           rows = 4)
      openxlsx::writeData(workbook, sheet, estimators[estimator_index],
                          start.col, 4, colNames = FALSE)
      openxlsx::writeData(
        workbook, sheet, matrix(c("Absolute bias", "RMSE", "Excess regret"), nrow = 1L),
        start.col, 5, colNames = FALSE, rowNames = FALSE
      )
    }

    body <- vector("list", length(n.values) * length(r.values))
    body.index <- 0L
    for (n.current in n.values) {
      for (r.current in r.values) {
        body.index <- body.index + 1L
        output.row <- c(n.current, r.current)
        for (estimator_current in estimators) {
          selected <- export$J == J.current & export$n == n.current &
            export$r == r.current & export$estimator == estimator_current
          if (sum(selected) != 1L) {
            stop("Missing/duplicate row: J=", J.current, ", n=", n.current,
                 ", r=", r.current, ", estimator=", estimator_current)
          }
          output.row <- c(
            output.row,
            clean_rounded(export$absolute_bias[selected]),
            clean_rounded(export$rmse[selected]),
            clean_rounded(export$excess_regret[selected])
          )
        }
        body[[body.index]] <- output.row
      }
    }
    body <- do.call(rbind, body)
    first.data.row <- 6L
    last.data.row <- first.data.row + nrow(body) - 1L
    openxlsx::writeData(workbook, sheet, body, 1, first.data.row,
                        colNames = FALSE, rowNames = FALSE)

    openxlsx::addStyle(workbook, sheet, title.style, rows = 1, cols = 1:14,
                       gridExpand = TRUE)
    openxlsx::addStyle(workbook, sheet, subtitle.style, rows = 2, cols = 1:14,
                       gridExpand = TRUE)
    openxlsx::addStyle(workbook, sheet, header.style, rows = 4:5, cols = 1:14,
                       gridExpand = TRUE)
    openxlsx::addStyle(workbook, sheet, body.style,
                       rows = first.data.row:last.data.row, cols = 1:14,
                       gridExpand = TRUE)
    openxlsx::addStyle(workbook, sheet, integer.style,
                       rows = first.data.row:last.data.row, cols = 1,
                       gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(workbook, sheet, rate.style,
                       rows = first.data.row:last.data.row, cols = 2,
                       gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(workbook, sheet, result.style,
                       rows = first.data.row:last.data.row, cols = 3:14,
                       gridExpand = TRUE, stack = TRUE)

    for (n.index in seq_along(n.values)) {
      group.row <- first.data.row + (n.index - 1L) * length(r.values)
      openxlsx::addStyle(workbook, sheet, separator.style,
                         rows = group.row, cols = 1:14,
                         gridExpand = TRUE, stack = TRUE)
      if (n.index %% 2L == 0L) {
        openxlsx::addStyle(
          workbook, sheet, alternate.style,
          rows = group.row:(group.row + length(r.values) - 1L), cols = 1:14,
          gridExpand = TRUE, stack = TRUE
        )
      }
    }
    openxlsx::addStyle(workbook, sheet, bottom.style,
                       rows = last.data.row, cols = 1:14,
                       gridExpand = TRUE, stack = TRUE)

    note.row <- last.data.row + 2L
    notes <- c(
      "Absolute bias = abs(mean(estimate - truth)); RMSE = sqrt(mean((estimate - truth)^2)).",
      "All error rates shown in the figures are included; divide reported metrics by 100 to match the figure axes.",
      "Excess regret is relative to the same-replication oracle depth-2 tree."
    )
    for (note.index in seq_along(notes)) {
      current.row <- note.row + note.index - 1L
      openxlsx::mergeCells(workbook, sheet, cols = 1:14, rows = current.row)
      openxlsx::writeData(workbook, sheet, notes[note.index],
                          1, current.row, colNames = FALSE)
      openxlsx::addStyle(workbook, sheet, note.style,
                         rows = current.row, cols = 1:14, gridExpand = TRUE)
      openxlsx::setRowHeights(workbook, sheet, rows = current.row, heights = 19)
    }

    openxlsx::setRowHeights(workbook, sheet, rows = 1, heights = 26)
    openxlsx::setRowHeights(workbook, sheet, rows = 2, heights = 22)
    openxlsx::setRowHeights(workbook, sheet, rows = 4, heights = 24)
    openxlsx::setRowHeights(workbook, sheet, rows = 5, heights = 34)
    openxlsx::setRowHeights(workbook, sheet,
                            rows = first.data.row:last.data.row, heights = 22)
    openxlsx::setColWidths(workbook, sheet, cols = 1, widths = 11)
    openxlsx::setColWidths(workbook, sheet, cols = 2, widths = 7)
    openxlsx::setColWidths(workbook, sheet,
                           cols = c(3, 4, 6, 7, 9, 10, 12, 13), widths = 12)
    openxlsx::setColWidths(workbook, sheet,
                           cols = c(5, 8, 11, 14), widths = 15)
    openxlsx::freezePane(workbook, sheet, firstActiveRow = 6,
                         firstActiveCol = 3)
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(workbook, file, overwrite = TRUE)
  if (!file.exists(file)) stop("Excel export returned without creating: ", file)
  invisible(file)
}

anchor_legend_to_first_panel <- function(plot, J.current) {
  legend.table <- ggplot2::ggplotGrob(
    plot + ggplot2::theme(legend.position = "right")
  )
  legend.index <- which(legend.table$layout$name == "guide-box-right")
  if (length(legend.index) != 1L ||
      inherits(legend.table$grobs[[legend.index]], "zeroGrob"))
    stop("Could not extract the plot legend.")
  legend.grob <- legend.table$grobs[[legend.index]]

  panel.table <- ggplot2::ggplotGrob(
    plot + ggplot2::theme(legend.position = "none")
  )
  panel.index <- which(panel.table$layout$name == "panel-1-1")
  if (length(panel.index) != 1L)
    stop("Could not identify panel (1,1).")
  panel.cell <- panel.table$layout[panel.index, ]

  if (J.current == 8) {
    legend.y <- grid::unit(1, "npc")
    legend.justification <- c("left", "top")
  } else {
    ## Lift the J3/J5 legend slightly so its last row clears the panel border.
    legend.y <- grid::unit(0, "npc") + grid::unit(1.5, "mm")
    legend.justification <- c("left", "bottom")
  }

  ## The viewport uses the legend grob's full bounding box. J8 is anchored by
  ## its upper-left corner; J3/J5 use the lower-left corner with a small lift.
  anchored.legend <- grid::grobTree(
    ## Use a transparent rectangle for the legend's background area.
    grid::rectGrob(gp = grid::gpar(
      fill = grDevices::adjustcolor("white", alpha.f = 0), col = NA
    )),
    legend.grob,
    vp = grid::viewport(
      x = grid::unit(0, "npc"), y = legend.y,
      width = grid::grobWidth(legend.grob),
      height = grid::grobHeight(legend.grob),
      just = legend.justification
    )
  )
  gtable::gtable_add_grob(
    panel.table, anchored.legend,
    t = panel.cell$t, l = panel.cell$l,
    b = panel.cell$b, r = panel.cell$r,
    z = Inf, clip = "off", name = "anchored-panel-1-1-legend"
  )
}

## One combined 3-row ggplot2 facet grid for each J. Columns are n; rows are
## excess regret, absolute bias and RMSE. All rows share the same error-rate
## positions and each metric shares its y scale across sample-size columns.
write_Rsup_panel_pdf <- function(summary, J.current, file) {
  z <- summary[summary$J == J.current, , drop = FALSE]
  estimators <- unname(estimator_labels)
  legend.labels <- estimators
  n.values <- sort(unique(z$n))
  r.values <- sort(unique(z$r))
  x.breaks <- (1:5) / 10
  expected_rows <- length(n.values) * length(r.values) * length(estimators)
  if (!nrow(z) || nrow(z) != expected_rows ||
      anyDuplicated(z[c("n", "r", "estimator")]) ||
      !setequal(unique(z$estimator), estimators))
    stop("Incomplete estimator-n-r grid for J=", J.current)

  n.labels <- sprintf("n = %d", as.integer(n.values))
  plot.data <- dplyr::bind_rows(
    z %>% transmute(
      n = factor(sprintf("n = %d", as.integer(n)), levels = n.labels),
      r, estimator = factor(estimator, levels = estimators),
      metric = "Excess regret", value = excess_regret
    ),
    z %>% transmute(
      n = factor(sprintf("n = %d", as.integer(n)), levels = n.labels),
      r, estimator = factor(estimator, levels = estimators),
      metric = "Absolute bias", value = absolute_bias
    ),
    z %>% transmute(
      n = factor(sprintf("n = %d", as.integer(n)), levels = n.labels),
      r, estimator = factor(estimator, levels = estimators),
      metric = "RMSE", value = rmse
    )
  ) %>%
    mutate(metric = factor(metric,
                           levels = c("Excess regret", "Absolute bias", "RMSE")))

  y_breaks <- function(limits) {
    ticks <- pretty(limits, n = 5)
    ticks[ticks >= 0 & ticks <= max(limits)]
  }
  format_ticks <- function(values) {
    nonzero <- abs(values[is.finite(values) & values != 0])
    digits <- if (!length(nonzero)) 2L else
      max(2L, min(3L, ceiling(-log10(min(nonzero)))))
    formatC(values, format = "f", digits = digits)
  }
  colors <- stats::setNames(c("#0072B2", "#8E44AD", "#009E73", "#D62728"),
                            estimators)
  shapes <- stats::setNames(c(16, 18, 15, 17), estimators)
  line.types <- stats::setNames(c("dashed", "dotdash", "dotted", "solid"),
                                estimators)

  p <- ggplot2::ggplot(
    plot.data,
    ggplot2::aes(x = r, y = value, color = estimator,
                 linetype = estimator, shape = estimator, group = estimator)
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 2.15, stroke = 0.35) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(metric), cols = ggplot2::vars(n),
      scales = "free_y", labeller = ggplot2::label_value
    ) +
    ggplot2::scale_x_continuous(
      breaks = x.breaks, labels = sprintf("%.1f", x.breaks),
      limits = range(c(x.breaks, r.values)),
      expand = ggplot2::expansion(mult = c(0.06, 0.06))
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks, labels = format_ticks,
      expand = ggplot2::expansion(mult = c(0.06, 0.08))
    ) +
    ggplot2::scale_color_manual(
      values = colors, breaks = estimators, labels = legend.labels
    ) +
    ggplot2::scale_linetype_manual(
      values = line.types, breaks = estimators, labels = legend.labels
    ) +
    ggplot2::scale_shape_manual(
      values = shapes, breaks = estimators, labels = legend.labels
    ) +
    ggplot2::labs(x = expression(paste("Error rate ", italic(r))), y = NULL,
                  color = NULL, linetype = NULL, shape = NULL) +
    ggplot2::theme_bw(base_size = 13, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(size = 13, color = "#111827"),
      axis.text = ggplot2::element_text(size = 13, color = "#111827"),
      axis.text.x = ggplot2::element_text(size = 13, angle = 0, hjust = 0.5),
      axis.title = ggplot2::element_text(size = 13, color = "#111827"),
      strip.text = ggplot2::element_text(size = 13, face = "bold", color = "#111827"),
      strip.background = ggplot2::element_rect(
        fill = "#D9D9D9", color = "#4B5563", linewidth = 0.6
      ),
      panel.border = ggplot2::element_rect(
        fill = NA, color = "#4B5563", linewidth = 0.6
      ),
      panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10, color = "#111827"),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      legend.key.height = grid::unit(12, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
      plot.margin = ggplot2::margin(6, 6, 6, 6, unit = "pt")
    )

  grDevices::cairo_pdf(file, width = 9.3, height = 7.8,
                       family = "Helvetica", pointsize = 13)
  device <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(device), add = TRUE)
  panel.figure <- anchor_legend_to_first_panel(p, J.current)
  grid::grid.draw(panel.figure)
  invisible(file)
}

run_Rsup_tables <- function(export_xlsx = TRUE, export_pdf = TRUE) {
  require_user_path(rsup_data_dir, "dir_out in data_gen.R")
  if (export_xlsx || export_pdf)
    require_user_path(publication_dir, "publication_dir in sim_Rsup_tab.R")
  result_file <- file.path(rsup_data_dir, paste0(output_prefix, "Rsup.rds"))
  if (!file.exists(result_file)) stop("Run sim_Rsup.R first: ", result_file)
  results <- readRDS(result_file)
  if (!identical(attr(results, "run_config"), run_config))
    stop("Result configuration differs from the current profile.")

  expected <- expand.grid(r = rate, n = n, J = J)
  if (nrow(results) != nrow(expected) * R)
    stop("Simulation is incomplete; finish/resume sim_Rsup.R before exporting.")
  for (i in seq_len(nrow(expected))) {
    reps <- results$rep[results$n == expected$n[i] &
                        results$J == expected$J[i] &
                        results$rate == expected$r[i]]
    if (!identical(sort(reps), as.numeric(seq_len(R))))
      stop("Missing or duplicate repetitions in scenario ", i)
  }

  summary <- build_Rsup_summary(results)
  if (export_xlsx) {
    output_file <- publication_path("Rsup_tab.xlsx")
    write_Rsup_excel(summary, output_file)
    message("Saved: ", output_file)
  }
  if (export_pdf) {
    for (J.current in sort(unique(summary$J))) {
      output_file <- publication_path(sprintf("Rsup_plot_J%d.pdf", J.current))
      write_Rsup_panel_pdf(summary, J.current, output_file)
      message("Saved: ", output_file)
    }
  }
  invisible(summary)
}

if (!isTRUE(getOption("opl.Rsup.functions_only", FALSE))) run_Rsup_tables()
