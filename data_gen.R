suppressPackageStartupMessages({
  library(foreach)
  library(doFuture)
  library(future)
  library(doRNG)
  library(progressr)
})

## Locate this script for source provenance with source() or Rscript.
script_file <- local({
  source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(source_files)) tail(source_files, 1)[[1]] else
    if (length(arg)) sub("^--file=", "", arg[[1]]) else
      stop("Run with source('/path/to/this_script.R') or Rscript.")
})


root_dir <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)
run_profile <- Sys.getenv("OPL_RSUP_PROFILE", "full")
if (!grepl("^[A-Za-z0-9_-]+$", run_profile)) stop("Invalid run profile name.")
profile_defaults <- switch(run_profile,
  full = list(R = 500L, n = c(500L, 1000L, 5000L), J = c(3L, 5L, 8L),
              rate = 0.10 + 0.05 * (0:8)),
  pilot = list(R = 20L, n = 2000L, J = 5L, rate = (1:6) / 10),
  smoke = list(R = 2L, n = 200L, J = 5L, rate = c(.3, .5)),
  list(R = 500L, n = c(500L, 1000L, 5000L), J = c(3L, 5L, 8L),
       rate = 0.10 + 0.05 * (0:8)))
R <- as.integer(getOption("opl.Rsup.R", profile_defaults$R))
n <- as.integer(getOption("opl.Rsup.n", profile_defaults$n))
J <- as.integer(getOption("opl.Rsup.J", profile_defaults$J))
rate <- as.numeric(getOption("opl.Rsup.rate", profile_defaults$rate))
p <- 2L
Cu <- 0.35
h <- 2
tree_depth <- 2L
data_seed <- 123L
## Separate stage seeds: restarting doRNG with the data seed would reuse the
## data-generating random streams as nuisance noise. Keep these seeds distinct.
simulation_seed <- 67891L
propensity_clip <- 0.10
worker_cap <- 56L
total <- 2L * R
save_every <- as.integer(getOption("opl.Rsup.save_every", 10L))
stopifnot(length(R) == 1L, R >= 2L, all(n >= 20L), all(J %in% c(3L, 5L, 8L)),
          all(rate > 0 & rate <= .6), data_seed != simulation_seed,
          !anyDuplicated(n), !anyDuplicated(J), !anyDuplicated(rate),
          length(save_every) == 1L, !is.na(save_every), save_every >= 1L)
## USER PATHS: replace each empty string with your absolute directory path.
## Use forward slashes, e.g. "D:/project/results" or "/path/to/results".
## dir_in stores generated datasets; dir_out stores the regret results.
## sim_Rsup.R and sim_Rsup_tab.R inherit these paths. R options may override them.
dir_out <- path.expand(getOption("opl.Rsup.result_dir", ""))
dir_in <- path.expand(getOption("opl.Rsup.data_dir", ""))
require_user_path <- function(path, setting) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("Set ", setting, " to your absolute path before running this step.")
  path <- gsub("\\\\", "/", path.expand(path))
  if (!grepl("^(/|[A-Za-z]:/)", path))
    stop(setting, " must be an absolute path; use forward slashes on Windows.")
  path
}
## Keep pilot/smoke files separate without introducing result subdirectories.
output_prefix <- if (run_profile == "full") "" else paste0(run_profile, "_")
result_path <- function(name) file.path(require_user_path(dir_out,
  "the result directory in data_gen.R or sim_trPr.R"), paste0(output_prefix, name))
run_config <- list(version = "0903-rsup-v4-clear-names", profile = run_profile,
                   n = n, J = J, R = R, rate = rate, p = p, Cu = Cu, h = h,
                   assignment_mechanism = "pmin(pmax(0.1, X1^2), 0.9)",
                   tree_depth = tree_depth,
                   nuisance_means = c(treatment = 1, control = -1),
                   propensity_nuisance_mean = 1,
                   propensity_clip = propensity_clip,
                   direct_if = "hard active branches with outer zero-truncation indicators",
                   normalize_probabilities = FALSE,
                   data_seed = data_seed, simulation_seed = simulation_seed)

path_in <- function(n, J) {
  file.path(require_user_path(dir_in, "the data directory in data_gen.R or sim_trPr.R"),
    paste0(output_prefix, sprintf("dat_J%d_n%d.rds", J, n)))
}
data_metadata <- function(n_current, J_current) {
  list(version = "0903-observational-v1", n = n_current, J = J_current, p = p,
       R = R, Cu = Cu,
       assignment_mechanism = "pmin(pmax(0.1, X1^2), 0.9)", seed = data_seed,
       grid_n = n, grid_J = J)
}

## Base-R storage avoids an additional serialization dependency.
## Write/read verification and a backup protect already generated results.
recover_rds_backup <- function(file) {
  backup <- paste0(file, ".bak")
  if (!file.exists(backup)) return(invisible(file))
  if (file.exists(file)) {
    ## A readable new file means the commit finished before backup cleanup.
    tryCatch(readRDS(file), error = function(e)
      stop("Output and backup coexist, but output is unreadable; inspect: ", file))
    if (!file.remove(backup)) stop("Cannot remove committed output backup: ", backup)
  } else if (!file.rename(backup, file)) {
    stop("Cannot recover output backup: ", backup)
  }
  invisible(file)
}

atomic_rds_save <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  recover_rds_backup(file)
  backup <- paste0(file, ".bak")
  tmp <- tempfile(pattern = paste0(basename(file), "_"), tmpdir = dirname(file))
  on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
  saveRDS(x, tmp, compress = FALSE)
  if (!identical(readRDS(tmp), x)) stop("Temporary output failed verification: ", tmp)
  if (file.exists(backup)) stop("Unresolved backup exists; inspect first: ", backup)
  if (file.exists(file) && !file.rename(file, backup)) stop("Cannot back up: ", file)
  if (!file.rename(tmp, file)) {
    if (file.exists(backup)) file.rename(backup, file)
    stop("Cannot finalize: ", file)
  }
  if (file.exists(backup) && !file.remove(backup)) warning("Retained backup: ", backup)
  invisible(file)
}


## Near-live progress reaches the main process. Only this process writes the
## combined Rsup result; workers send small batches through progress conditions.
scenario_progress_handler <- function(on_packet = NULL) {
  reporter <- list(update = function(config, state, progression, ...) {
    if (!is.null(progression$rsup_packet) && !is.null(on_packet))
      on_packet(progression$rsup_packet)
    detail <- paste(progression$message, collapse = " ")
    if (!nzchar(detail)) return(invisible(NULL))
    line <- sprintf("[%s] [%d/%d scenarios] %s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    as.integer(state$step), as.integer(config$max_steps), detail)
    cat(line, "\n")
    flush.console()
  })
  progressr::make_progression_handler("scenario_lines", reporter,
    interval = 0, times = Inf, intrusiveness = 0, clear = FALSE, enable = TRUE)
}

valid_rng_state <- function(x) {
  is.integer(x) && length(x) == 7L && !anyNA(x)  # L'Ecuyer-CMRG from doRNG
}
restore_rng_state <- function(x) {
  if (!valid_rng_state(x)) stop("Invalid saved parallel RNG state.")
  assign(".Random.seed", x, envir = .GlobalEnv)
}
execution_info <- function(stage, workers) {
  files <- file.path(root_dir, c("data_gen.R", "sim_Rsup.R", "sim_Rsup_tab.R"))
  list(stage = stage, run_config = run_config,
                      R = R.version.string, session = sessionInfo(),
                      source_md5 = tools::md5sum(files),
                      workers_used = workers,
                      save_every = save_every, recorded_at = Sys.time())
}

simulation_workers <- function(scenarios) {
  configured <- Sys.getenv("OPL_RSUP_WORKERS", "")
  requested <- if (nzchar(configured)) suppressWarnings(as.integer(configured)) else
    as.integer(future::availableCores())
  if (length(requested) != 1L || is.na(requested) || requested < 1L)
    stop("OPL_RSUP_WORKERS must be a positive integer.")
  min(as.integer(scenarios), requested, worker_cap)
}

## Lower-bound branches for eta, shifted by the utility threshold Cu.
compute_delta_lower <- function(m1, m0, J, Cu) {
  eta_lower <- sapply(1:J, function(j) {
    rowSums(m1[, j:J, drop = FALSE]) - rowSums(m0[, j:J, drop = FALSE])
  })
  eta_lower - Cu
}

## Upper-bound branches for eta, shifted by the utility threshold Cu.
compute_delta_upper <- function(m1, m0, J, Cu) {
  eta_upper <- sapply(1:J, function(j) {
    if (j < J) {
      rowSums(m1[, (j + 1):J, drop = FALSE]) - rowSums(m0[, j:J, drop = FALSE]) + 1
    } else {
      0 - rowSums(m0[, j:J, drop = FALSE]) + 1
    }
  })
  eta_upper - Cu
}

getP <- function(a, X, J) {
  coef.rand <- c(-0.3, -0.1, 0.9, -0.6, -0.7, 0.2, 0.7, 0.5, 0.3, 0.7, -0.4, -0.1, 0, 0.2)
  coef.sele <- coef.rand[1:(p * (J - 1))]
  alpha <- matrix(coef.sele, nrow = p, ncol = J - 1)
  alpha[p, ] <- (2 * a - 1) * alpha[p, ]
  interc <- seq(-0.4, -0.1, by = 0.05)  # intercept
  lin <- cbind(1, X) %*% rbind(interc[1:(J - 1)], alpha)
  exp_lin <- exp(lin)
  den <- 1 + rowSums(exp_lin)
  P <- cbind(1 / den, exp_lin / den)
  P
}

## Vectorized inverse-CDF draw from one categorical distribution per row.
draw_ordinal <- function(prob) {
  n_obs <- nrow(prob)
  J_obs <- ncol(prob)
  u <- runif(n_obs)
  y <- rep.int(J_obs, n_obs)
  cumulative <- numeric(n_obs)
  undecided <- rep.int(TRUE, n_obs)

  for (j in seq_len(J_obs - 1L)) {
    cumulative <- cumulative + prob[, j]
    selected <- undecided & (u <= cumulative)
    y[selected] <- j
    undecided[selected] <- FALSE
  }
  y
}

data_gen <- function(n, J) {
  X <- matrix(runif(n * p, -1, 1), n, p)
  e.true <- pmin(pmax(0.1, X[, 1]^2), 0.9)
  A <- rbinom(n, 1, e.true)
  prob_Y1 <- getP(1, X, J)
  prob_Y0 <- getP(0, X, J)
  Y0 <- draw_ordinal(prob_Y0)
  Y1 <- draw_ordinal(prob_Y1)
  Y <- A * Y1 + (1 - A) * Y0
  
  delta_L <- compute_delta_lower(prob_Y1, prob_Y0, J, Cu)
  delta_U <- compute_delta_upper(prob_Y1, prob_Y0, J, Cu)
  psi.L <- apply(delta_L, 1, max)
  psi.U <- apply(delta_U, 1, min)
  psi.tilde <- pmax(psi.U, 0) + pmin(psi.L, 0)
  
  list(X = X, A = A, Y = Y, prob_Y1 = prob_Y1, prob_Y0 = prob_Y0, 
       psi.tilde = psi.tilde, psi.L = psi.L, psi.U = psi.U,
       e.true = e.true)
}

run_data_generation <- function() {
  require_user_path(dir_in, "dir_in in data_gen.R")
  require_user_path(dir_out, "dir_out in data_gen.R")
  params <- expand.grid(n = n, J = J)
  workers <- simulation_workers(nrow(params))
  registerDoFuture()
  if (workers == 1L) plan(sequential) else plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  info <- execution_info("data", workers = workers)
  start_time <- Sys.time()
  message(sprintf("Data: profile=%s, workers=%d, scenarios=%d; final datasets only",
                  run_profile, workers, nrow(params)))
  files <- progressr::with_progress({
    progress <- progressr::progressor(steps = nrow(params))
    foreach(t = seq_len(nrow(params)), .combine = "c", .inorder = TRUE,
            .options.future = list(chunk.size = 1L),
            .options.RNG = data_seed) %dorng% {
      n_current <- params$n[t]
      J_current <- params$J[t]
      output_file <- path_in(n_current, J_current)
      expected <- data_metadata(n_current, J_current)
      recover_rds_backup(output_file)
      if (file.exists(output_file)) {
        saved <- readRDS(output_file)
        if (!identical(saved$metadata, expected) || length(saved$dats) != total)
          stop("Existing data have different settings; archive them or choose a new run profile: ", output_file)
        progress(message = sprintf("REUSED data J=%d n=%d: %s",
                                   J_current, n_current, basename(output_file)))
      } else {
        dats <- lapply(seq_len(total), function(tt) data_gen(n_current, J_current))
        atomic_rds_save(list(metadata = expected, dats = dats,
                             execution_info = info), output_file)
        progress(message = sprintf("DONE data J=%d n=%d: saved %s",
                                   J_current, n_current, basename(output_file)))
      }
      output_file
    }
  }, handlers = scenario_progress_handler(),
     enable = TRUE, delay_stdout = FALSE, delay_conditions = character())
  stopifnot(length(files) == nrow(params), all(file.exists(files)))
  print(Sys.time() - start_time)
  invisible(files)
}

if (!isTRUE(getOption("opl.Rsup.functions_only", FALSE))) run_data_generation()
