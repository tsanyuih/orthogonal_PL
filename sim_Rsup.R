
## version 4.6.1 (2026-06-24)
## Operating system: Windows 11 (64-bit)
## Intel(R) Xeon(R) Gold 6226R CPU @ 2.90GHz
## Approximate runtime on the above system: 1 hour

suppressPackageStartupMessages({
  library(foreach)
  library(doFuture)
  library(future)
  library(progressr)
  library(doRNG)
  library(policytree)
})

## Locate this script to load shared helpers with source() or Rscript.
script_file <- local({
  source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(source_files)) tail(source_files, 1)[[1]] else
    if (length(arg)) sub("^--file=", "", arg[[1]]) else
      stop("Run with source('/path/to/this_script.R') or Rscript.")
})
## USER PATHS: set dir_in and dir_out in data_gen.R before running this script.
## Keep shared helper scripts beside this file; they are located automatically.
## Load shared settings and utilities without running data generation.
previous_source_options <- options(opl.Rsup.functions_only = TRUE)
tryCatch(
  source(file.path(dirname(normalizePath(script_file, winslash = "/")), "data_gen.R")),
  finally = options(previous_source_options))

## G^beta is the softmax-weighted average of the input values.
## compute_scores() combines the original and zero-augmented versions.
smooth_max <- function(values, beta) {
  z <- beta * values
  w <- exp(z - max(z))
  w <- w / sum(w)
  sum(w * values)
}

## Gradient of pure G^beta with respect to its input values.
smooth_max_grad <- function(values, beta) {
  z <- beta * values
  w <- exp(z - max(z))
  w <- w / sum(w)
  g <- sum(w * values)
  w * (1 + beta * (values - g))
}

## Calculate the four estimator scores.
compute_scores <- function(J, Cu, rate, dat) {
  A <- dat$A
  Y <- dat$Y
  n <- length(Y)
  m1_hat <- plogis(qlogis(dat$prob_Y1) + h * rnorm(n * J, 1 / n^rate, 1 / n^rate))
  m0_hat <- plogis(qlogis(dat$prob_Y0) + h * rnorm(n * J, -1 / n^rate, 1 / n^rate))
  if (is.null(dat$e.true) || length(dat$e.true) != n ||
      any(!is.finite(dat$e.true)) || any(dat$e.true <= 0 | dat$e.true >= 1))
    stop("Missing or invalid true propensity scores; regenerate data with data_gen.R.")
  ## Projection preserves the controlled n^(-r) rate and enforces the fixed
  ## estimated-positivity condition used by the theoretical results.
  e_hat <- plogis(qlogis(dat$e.true) + h * rnorm(n, 1 / n^rate, 1 / n^rate))
  e_hat <- pmin(pmax(e_hat, propensity_clip), 1 - propensity_clip)
  
  delta_L <- compute_delta_lower(m1_hat, m0_hat, J, Cu)
  delta_U <- compute_delta_upper(m1_hat, m0_hat, J, Cu)
  active_L <- apply(delta_L, 1, which.max)
  active_U <- apply(delta_U, 1, which.min)
  row_index <- seq_len(n)
  max_L <- delta_L[cbind(row_index, active_L)]
  min_U <- delta_U[cbind(row_index, active_U)]
  
  Y_ind <- outer(Y, seq_len(J), FUN = "==")
  if_a1 <- A / e_hat * (Y_ind - m1_hat)
  if_a0 <- (1 - A) / (1 - e_hat) * (Y_ind - m0_hat)

  ## Influence-function corrections for each lower- and upper-bound branch.
  if_L <- sapply(1:J, function(j) {
    rowSums(if_a1[, j:J, drop = FALSE]) -
      rowSums(if_a0[, j:J, drop = FALSE])
  })
  if_U <- sapply(1:J, function(j) {
    treated <- if (j < J) {
      rowSums(if_a1[, (j + 1):J, drop = FALSE])
    } else {
      rep(0, n)
    }
    treated - rowSums(if_a0[, j:J, drop = FALSE])
  })
  
  ## Plug-in estimator.
  score_plugin <- pmax(min_U, 0) + pmin(max_L, 0)
  
  ## Smoothed estimators.
  beta <- 2 * h * n^(max(0.25, 2 * rate - 0.5))
  delta_L_aug <- cbind(delta_L, 0)
  delta_U_aug <- cbind(delta_U, 0)
  ## Step 1: softmax approximation of score function
  psi_L_base <- apply(delta_L, 1, smooth_max, beta)
  psi_L_aug <- apply(delta_L_aug, 1, smooth_max, beta)
  psi_L <- psi_L_base - psi_L_aug
  grad_L_base <- t(apply(delta_L, 1, smooth_max_grad, beta))
  grad_L_aug <- t(apply(delta_L_aug, 1, smooth_max_grad, beta))
  psi_U_base <- -apply(-delta_U, 1, smooth_max, beta)
  psi_U_aug <- -apply(-delta_U_aug, 1, smooth_max, beta)
  psi_U <- psi_U_base - psi_U_aug
  grad_U_base <- t(apply(delta_U, 1, smooth_max_grad, -beta))
  grad_U_aug <- t(apply(delta_U_aug, 1, smooth_max_grad, -beta))
  ## Step 2: adjustment term
  grad_psi_L <- grad_L_base - grad_L_aug[, 1:J]
  grad_psi_U <- grad_U_base - grad_U_aug[, 1:J]
  outcome_residual <- Y_ind - (A * m1_hat + (1 - A) * m0_hat)
  tail_residual <- sapply(1:J, function(j) {
    rowSums(outcome_residual[, j:J, drop = FALSE])
  })
  contrast_weight <- (A - e_hat) / (e_hat * (1 - e_hat))
  treated_weight <- A / e_hat
  weighted_tail <- rowSums((grad_psi_U + grad_psi_L) * tail_residual)
  weighted_category <- rowSums(grad_psi_U * outcome_residual)
  if_correction <- contrast_weight * weighted_tail - treated_weight * weighted_category
  score_smoothed_plugin <- psi_L + psi_U  # ablation: same smoothing and same nuisances
  score_orthogonal_smoothed <- score_smoothed_plugin + if_correction
  
  ## Levis-style direct IF-based estimator for the nonsmoothed score.
  ## The IF correction is applied only to the active hard max/min branch;
  ## the outer truncation contributes only when that branch is below/above zero.
  active_if_L <- if_L[cbind(row_index, active_L)]
  active_if_U <- if_U[cbind(row_index, active_U)]
  direct_if_correction <- (max_L < 0) * active_if_L + (min_U > 0) * active_if_U
  score_direct_if <- score_plugin + direct_if_correction
  
  list(plugin = score_plugin,
       direct_if = score_direct_if,
       smoothed_plugin = score_smoothed_plugin,
       orthogonal_smoothed = score_orthogonal_smoothed)
}

evaluate_tree_estimator <- function(train, test, score_train, score_test) {
  rewards_train <- cbind(0, score_train)
  tree <- policy_tree(train$X, rewards_train, depth = tree_depth)
  policy_hat <- predict(tree, test$X) - 1
  Rsup_estimate <- -mean(policy_hat * score_test)
  Rsup_truth <- -mean(policy_hat * test$psi.tilde)
  
  treat_rate <- mean(policy_hat == 1)
  oracle_policy <- (test$psi.tilde > 0)
  misclass_rate <- mean(policy_hat != oracle_policy)
  
  c(Rsup_estimate = Rsup_estimate, Rsup_truth = Rsup_truth,
    treat_rate = treat_rate, misclass_rate = misclass_rate)
}

oracle_tree_risk <- function(train, test) {
  ## Internal benchmark for excess regret; not a fifth estimator output.
  rewards_train <- cbind(0, train$psi.tilde)
  tree <- policy_tree(train$X, rewards_train, depth = tree_depth)
  policy_hat <- predict(tree, test$X) - 1
  -mean(policy_hat * test$psi.tilde)
}

estimator_names <- c("plugin", "direct_if", "smoothed_plugin", "orthogonal_smoothed")
result_columns <- c("n", "J", "rate", "rep",
  unlist(lapply(estimator_names, function(estimator) {
    paste(estimator,
          c("Rsup_estimate", "Rsup_truth", "treat_rate", "misclass_rate", "excess"),
          sep = ".")
  }), use.names = FALSE))

valid_result <- function(x, expected_rows) {
  is.data.frame(x) && nrow(x) == expected_rows &&
    identical(names(x), result_columns) &&
    all(vapply(x, function(v) is.numeric(v) && all(is.finite(v)), logical(1)))
}

valid_scenario <- function(x, rows, n_current, J_current, rate_current, first = 1L) {
  valid_result(x, rows) && all(x$n == n_current) && all(x$J == J_current) &&
    all(x$rate == rate_current) &&
    identical(x$rep, as.numeric(seq.int(from = first, length.out = rows)))
}

plain_result <- function(x) {
  attributes(x) <- attributes(x)[c("names", "row.names", "class")]
  rownames(x) <- NULL
  x
}

empty_result <- function() {
  as.data.frame(setNames(rep(list(numeric()), length(result_columns)), result_columns),
                check.names = FALSE)
}


run_Rsup <- function() {
  require_user_path(dir_in, "dir_in in data_gen.R")
  require_user_path(dir_out, "dir_out in data_gen.R")
  params <- expand.grid(rate = rate, n = n, J = J)
  scenario_count <- nrow(params)
  input_paths <- unique(mapply(path_in, n = params$n, J = params$J))
  for (path in input_paths) recover_rds_backup(path)
  if (!all(file.exists(input_paths)))
    stop("Generate the observational data first by running data_gen.R in this profile.")

  result_file <- result_path("Rsup.rds")
  recover_rds_backup(result_file)
  parts <- rep(list(NULL), scenario_count)
  completed <- integer(scenario_count)
  rng_states <- rep(list(NULL), scenario_count)
  if (file.exists(result_file)) {
    saved <- readRDS(result_file)
    if (!identical(attr(saved, "run_config"), run_config) ||
        !valid_result(saved, nrow(saved)))
      stop("Incompatible combined result; archive it or choose a new profile: ", result_file)
    for (t in seq_len(scenario_count)) {
      selected <- saved$n == params$n[t] & saved$J == params$J[t] &
        saved$rate == params$rate[t]
      part <- plain_result(saved[selected, , drop = FALSE])
      completed[t] <- nrow(part)
      if (completed[t] > R ||
          !valid_scenario(part, completed[t], params$n[t], params$J[t], params$rate[t]))
        stop("Invalid or noncontiguous repetitions in combined result, scenario ", t)
      if (completed[t] > 0L) parts[[t]] <- part
    }
    if (sum(completed) != nrow(saved)) stop("Unknown scenarios in combined result.")
    state <- attr(saved, "resume_state")
    if (is.null(state)) {
      if (!all(completed == R))
        stop("Incomplete result has no embedded resume state: ", result_file)
    } else {
      if (!identical(state$layout_version, 1L) ||
          !identical(state$scenarios, params) ||
          !identical(state$completed, completed) ||
          !is.list(state$rng_states) || length(state$rng_states) != scenario_count)
        stop("Incompatible embedded resume state: ", result_file)
      rng_states <- state$rng_states
      for (t in which(completed > 0L & completed < R))
        if (!valid_rng_state(rng_states[[t]]))
          stop("Missing RNG state for unfinished scenario ", t)
    }
  }

  workers <- if (all(completed == R)) 0L else simulation_workers(scenario_count)
  info <- execution_info("simulation", workers = workers)
  make_snapshot <- function() {
    present <- which(completed > 0L)
    out <- if (length(present)) do.call(rbind, parts[present]) else empty_result()
    rownames(out) <- NULL
    attr(out, "run_config") <- run_config
    attr(out, "resume_state") <- list(layout_version = 1L, scenarios = params,
      completed = completed, rng_states = rng_states, complete = all(completed == R))
    attr(out, "execution_info") <- info
    out
  }
  persist <- function() atomic_rds_save(make_snapshot(), result_file)

  accept_packet <- function(packet) {
    t <- packet$scenario
    if (length(t) != 1L || !is.numeric(t) || !is.finite(t) ||
        t != as.integer(t) || t < 1L || t > scenario_count)
      stop("Invalid incoming scenario index.")
    first <- packet$first
    last <- packet$last
    if (length(first) != 1L || length(last) != 1L ||
        !is.numeric(first) || !is.numeric(last) ||
        !is.finite(first) || !is.finite(last) ||
        first != as.integer(first) || last != as.integer(last) ||
        first < 1L || last < first || last > R ||
        !valid_scenario(packet$results, last - first + 1L,
                        params$n[t], params$J[t], params$rate[t], first) ||
        !valid_rng_state(packet$rng_state))
      stop("Invalid incoming result batch for scenario ", t)
    if (last <= completed[t]) {
      ## Future completion can deliver a second copy of an already saved packet.
      existing <- plain_result(parts[[t]][seq.int(first, last), , drop = FALSE])
      if (!identical(existing, packet$results))
        stop("Conflicting duplicate result batch for scenario ", t)
      return(invisible(NULL))
    }
    if (first != completed[t] + 1L)
      stop("Out-of-order result batch for scenario ", t)
    parts[[t]] <<- if (completed[t] == 0L) packet$results else
      plain_result(rbind(parts[[t]], packet$results))
    completed[t] <<- as.integer(last)
    rng_states[[t]] <<- packet$rng_state
    persist()
    invisible(NULL)
  }

  persist()
  if (all(completed == R)) {
    message("REUSED complete combined result: ", result_file)
    return(invisible(make_snapshot()))
  }
  ## Fixed starting snapshots: updates received by the writer do not alter
  ## what was sent to workers. doRNG still allocates every scenario's stream.
  resume_counts <- completed
  resume_states <- rng_states
  registerDoFuture()
  if (workers == 1L) plan(sequential) else plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)
  start_time <- Sys.time()
  message(sprintf("Rsup: profile=%s, R=%d, workers=%d, scenarios=%d; update one result every %d repetitions",
                  run_profile, R, workers, scenario_count, save_every))
  finished <- progressr::with_progress({
    progress <- progressr::progressor(steps = scenario_count)
    foreach(t = seq_len(scenario_count), .inorder = TRUE, .packages = "policytree",
            .options.future = list(chunk.size = 1L),
            .options.RNG = simulation_seed) %dorng% {
      n_current <- params$n[t]
      J_current <- params$J[t]
      rate_current <- params$rate[t]
      done <- resume_counts[t]
      packets <- list()
      if (done == R) {
        progress(message = sprintf("REUSED Rsup J=%d n=%d r=%.2f: %d/%d repetitions",
                                   J_current, n_current, rate_current, R, R))
      } else {
        input <- readRDS(path_in(n_current, J_current))
        if (!identical(input$metadata, data_metadata(n_current, J_current)) ||
            length(input$dats) != total)
          stop("Input data do not match the observational 0903 configuration.")
        dats <- input$dats
        if (done > 0L) {
          restore_rng_state(resume_states[[t]])
        }
        out <- vector("list", R)
        first <- done + 1L
        for (tt in seq.int(done + 1L, R)) {
          train <- dats[[tt]]
          test <- dats[[R + tt]]
          scores_train <- compute_scores(J_current, Cu, rate_current, train)
          scores_test <- compute_scores(J_current, Cu, rate_current, test)
          estimator_results <- setNames(lapply(estimator_names, function(estimator) {
            evaluate_tree_estimator(train, test,
                                    scores_train[[estimator]], scores_test[[estimator]])
          }), estimator_names)
          oracle_risk <- oracle_tree_risk(train, test)
          estimator_results <- lapply(estimator_results, function(result) {
            c(result, excess = result[["Rsup_truth"]] - oracle_risk)
          })
          out[[tt]] <- c(n = n_current, J = J_current, rate = rate_current, rep = tt,
                         unlist(estimator_results, use.names = TRUE))
          if (tt %% save_every == 0L || tt == R) {
            block <- as.data.frame(do.call(rbind, out[seq.int(first, tt)]), check.names = FALSE)
            rownames(block) <- NULL
            packet <- list(scenario = t, first = first, last = tt, results = block,
                           rng_state = get(".Random.seed", .GlobalEnv))
            packets[[length(packets) + 1L]] <- packet
            ## Silent batches still reach the main writer; print only when
            ## this scenario is complete and its final batch has been saved.
            progress(amount = as.integer(tt == R), rsup_packet = packet,
              message = if (tt == R)
                sprintf("DONE Rsup J=%d n=%d r=%.2f: %d/%d repetitions in %s",
                  J_current, n_current, rate_current, tt, R, basename(result_file)) else "")
            first <- tt + 1L
          }
        }
      }
      ## Retain a completion copy as an integrity check / delivery fallback.
      list(scenario = t, packets = packets)
    }
  }, handlers = scenario_progress_handler(on_packet = accept_packet),
     enable = TRUE, delay_stdout = FALSE, delay_conditions = character())

  for (item in finished) for (packet in item$packets) accept_packet(packet)
  if (!all(completed == R)) stop("Some repetitions were not received; saved results can be resumed.")
  results <- make_snapshot()
  stopifnot(valid_result(results, scenario_count * R),
            identical(readRDS(result_file), results))
  message("Verified combined result: ", result_file)
  print(Sys.time() - start_time)
  invisible(results)
}

if (!isTRUE(getOption("opl.Rsup.functions_only", FALSE))) run_Rsup()
