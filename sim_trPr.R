
## version 4.6.1 (2026-06-24)
## Operating system: Windows 11 (64-bit)
## Intel(R) Xeon(R) Gold 6226R CPU @ 2.90GHz
## Approximate runtime on the above system: 1.5 hours

trpr_script_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(files)) tail(files, 1)[[1]] else
    if (length(arg)) sub("^--file=", "", arg[[1]]) else
      stop("Run with source('/path/to/sim_trPr.R') or Rscript.")
})
## Keep option restoration local: nested sources also use previous_source_options.
local({
  saved_options <- options(opl.Rsup.functions_only = TRUE)
  on.exit(options(saved_options))
  source(file.path(dirname(normalizePath(trpr_script_file, winslash = "/")),
                   "sim_Rsup.R"), local = .GlobalEnv)
})
## USER PATHS: replace the empty strings with absolute directory paths.
## dir_in must match the dataset directory used in data_gen.R.
## dir_out stores trPr.rds and is inherited by sim_trPr_plot.R.
## Use forward slashes on Windows. R options may override these paths.
dir_in <- path.expand(getOption("opl.trPr.data_dir", ""))
dir_out <- path.expand(getOption("opl.trPr.result_dir", ""))
rsup_simulation_workers <- simulation_workers
simulation_workers <- function(scenarios) {
  configured <- Sys.getenv("OPL_TRPR_WORKERS", "")
  if (!nzchar(configured)) return(rsup_simulation_workers(scenarios))
  requested <- suppressWarnings(as.numeric(configured))
  if (length(requested) != 1L || !is.finite(requested) ||
      requested != floor(requested) || requested < 1 || requested > 64)
    stop("OPL_TRPR_WORKERS must be an integer from 1 to 64.")
  min(as.integer(scenarios), as.integer(requested))
}
Cu_grid <- as.numeric(getOption("opl.trPr.Cu",
  c(0, 0.1, seq(0.2, 0.7, by = 0.05), 0.8, 0.9, 1)))
sensitivity_rate <- as.numeric(getOption("opl.trPr.rate", 0.4))
stopifnot(length(Cu_grid) > 0L, all(is.finite(Cu_grid)),
          all(Cu_grid >= 0 & Cu_grid <= 1), !anyDuplicated(Cu_grid),
          length(sensitivity_rate) == 1L, is.finite(sensitivity_rate),
          sensitivity_rate > 0, sensitivity_rate <= .6)
Cu_grid <- sort(Cu_grid)
## Leave global Cu=.35 intact: data_metadata() must match the source datasets.
run_config <- c(run_config, list(sensitivity_version = "0903-trPr-v1",
  Cu_grid = Cu_grid, sensitivity_rate = sensitivity_rate,
  coupling = "shared observational data; independent nuisance streams by cost"))
source_names <- c("data_gen.R", "sim_Rsup.R", "sim_trPr.R")
run_config$sensitivity_source_md5 <- setNames(
  unname(tools::md5sum(file.path(root_dir, source_names))), source_names)
result_columns <- append(result_columns, "Cu", after = 3L)

with_cost <- function(dat, J, cost) {
  dat$psi.L <- apply(compute_delta_lower(dat$prob_Y1, dat$prob_Y0, J, cost), 1, max)
  dat$psi.U <- apply(compute_delta_upper(dat$prob_Y1, dat$prob_Y0, J, cost), 1, min)
  dat$psi.tilde <- pmax(dat$psi.U, 0) + pmin(dat$psi.L, 0)
  dat
}

valid_scenario <- function(x, rows, n_current, J_current, rate_current,
                           Cu_current, first = 1L) {
  valid_result(x, rows) && all(x$n == n_current) && all(x$J == J_current) &&
    all(x$rate == rate_current) && all(x$Cu == Cu_current) &&
    identical(x$rep, as.numeric(seq.int(from = first, length.out = rows)))
}

run_trPr <- function() {
  require_user_path(dir_in, "dir_in in sim_trPr.R")
  require_user_path(dir_out, "dir_out in sim_trPr.R")
  params <- expand.grid(Cu = Cu_grid, rate = sensitivity_rate, n = n, J = J)
  scenario_count <- nrow(params)
  input_paths <- unique(mapply(path_in, n = params$n, J = params$J))
  for (path in input_paths) recover_rds_backup(path)
  if (!all(file.exists(input_paths)))
    stop("Generate the observational data first by running data_gen.R in this profile.")

  result_file <- result_path("trPr.rds")
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
        saved$rate == params$rate[t] & saved$Cu == params$Cu[t]
      part <- plain_result(saved[selected, , drop = FALSE])
      completed[t] <- nrow(part)
      if (completed[t] > R ||
          !valid_scenario(part, completed[t], params$n[t], params$J[t], params$rate[t], params$Cu[t]))
        stop("Invalid or noncontiguous repetitions in combined result, scenario ", t)
      if (completed[t] > 0L) parts[[t]] <- part
    }
    if (sum(completed) != nrow(saved)) stop("Unknown scenarios in combined result.")
    state <- attr(saved, "resume_state")
    if (is.null(state)) {
      ## Complete results can be read without a saved resume state.
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

  ## Validate and save worker packets in the main process.
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
                        params$n[t], params$J[t], params$rate[t], params$Cu[t], first) ||
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
  message(sprintf("trPr: profile=%s, R=%d, workers=%d, scenarios=%d; update one result every %d repetitions",
                  run_profile, R, workers, scenario_count, save_every))
  finished <- progressr::with_progress({
    progress <- progressr::progressor(steps = scenario_count)
    foreach(t = seq_len(scenario_count), .inorder = TRUE, .packages = "policytree",
            .options.future = list(chunk.size = 1L),
            .options.RNG = simulation_seed) %dorng% {
      n_current <- params$n[t]
      J_current <- params$J[t]
      rate_current <- params$rate[t]
      Cu_current <- params$Cu[t]
      done <- resume_counts[t]
      packets <- list()
      if (done == R) {
        progress(message = sprintf("REUSED trPr J=%d n=%d r=%.2f Cu=%.2f: %d/%d repetitions",
                                   J_current, n_current, rate_current, Cu_current, R, R))
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
          train <- with_cost(dats[[tt]], J_current, Cu_current)
          test <- with_cost(dats[[R + tt]], J_current, Cu_current)
          scores_train <- compute_scores(J_current, Cu_current, rate_current, train)
          scores_test <- compute_scores(J_current, Cu_current, rate_current, test)
          estimator_results <- setNames(lapply(estimator_names, function(estimator) {
            evaluate_tree_estimator(train, test,
                                    scores_train[[estimator]], scores_test[[estimator]])
          }), estimator_names)
          oracle_risk <- oracle_tree_risk(train, test)
          estimator_results <- lapply(estimator_results, function(result) {
            c(result, excess = result[["Rsup_truth"]] - oracle_risk)
          })
          out[[tt]] <- c(n = n_current, J = J_current, rate = rate_current, Cu = Cu_current, rep = tt,
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
                sprintf("DONE trPr J=%d n=%d r=%.2f Cu=%.2f: %d/%d repetitions in %s",
                  J_current, n_current, rate_current, Cu_current, tt, R, basename(result_file)) else "")
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

if (!isTRUE(getOption("opl.trPr.functions_only", FALSE))) run_trPr()


