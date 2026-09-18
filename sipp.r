
## version 4.6.1 (2026-06-24)
## Operating system: Windows 11 (64-bit)
## Intel(R) Xeon(R) Gold 6226R CPU @ 2.90GHz
## Approximate runtime on the above system: 6 mins

sipp_script_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (length(files)) tail(files,1)[[1]] else if (length(arg))
    sub("^--file=","",arg[[1]]) else stop("Use source() or Rscript.")
})
sipp_root <- dirname(normalizePath(sipp_script_file,winslash="/"))
sipp_helper_root <- sipp_root

local({
  saved <- options(opl.Rsup.functions_only=TRUE)
  on.exit(options(saved))
  source(file.path(sipp_helper_root,"sim_Rsup.R"),local=.GlobalEnv)
})
suppressPackageStartupMessages(library(nnet))
sipp_covars <- c("age","female","race_black","race_asian","hispanic",
  "educ_cat","povratio","married","insured_any","work_limited")
sipp_labels <- c(plugin="Plug-in",direct_if="Direct IF",
  smoothed_plugin="Smoothed plug-in",orthogonal_smoothed="Orthogonal smoothed")
sipp_diagnostic_columns <- c("ps_min","ps_max","ps_clipped_low","ps_clipped_high")
## Statistical and execution settings live here.
## R options do not override this block. Keep statistical settings unchanged
## when resuming a run.
sipp_parameters <- list(
  R=100L,
  Cu=sort(unique(c(seq(10L,60L,by=5L),24L:34L)))/100,
  workers=56L,               # USER SETTING: reduce to suit your machine (1-56).
  save_every=5L,               # Save every five completed repetition-Cu pairs.
  progress_every=50L,          # Total progress independent of cost completion.
  K=3L, J=5L, depth=2L,
  beta_scale=2, beta_power=.25, clip=.05, seed=123L
)


## USER PATH: replace the empty string below with your absolute SIPP results directory.
sipp_output_dir <- function() require_user_path(getOption("opl.sipp.result_dir", ""),
  "the result directory in sipp_output_dir() in sipp.r")
sipp_settings <- function() {
  cfg <- c(list(version="0903-sipp-december-v3-all-trees"),
    sipp_parameters[!names(sipp_parameters) %in% c("workers","save_every","progress_every")],
    list(estimators=names(sipp_labels),allocation="cyclic-1-nuisance-1-train-1-test",
      fold_stratification="treatment-by-outcome-before-bootstrap",
      bootstrap="person-row-within-fold"))
  cfg$Cu <- sort(cfg$Cu)
  stopifnot(length(cfg$R)==1L,is.finite(cfg$R),cfg$R>=2,cfg$R==floor(cfg$R),
    length(cfg$Cu)>0,all(is.finite(cfg$Cu)),all(cfg$Cu>=0 & cfg$Cu<=1),
    !anyDuplicated(cfg$Cu),cfg$K==3L,cfg$J==5L,
    length(cfg$beta_scale)==1L,is.finite(cfg$beta_scale),cfg$beta_scale>0,
    length(cfg$beta_power)==1L,is.finite(cfg$beta_power),cfg$beta_power>0,
    length(cfg$clip)==1L,is.finite(cfg$clip),cfg$clip>0,cfg$clip<.5)
  cfg
}
prepare_sipp <- function(cfg) {
  ## USER PATH: replace the empty string with the absolute filename of
  ## sipp_data_dec.rds, including its extension. Use forward slashes on Windows.
  input <- require_user_path(getOption("opl.sipp.data_file", ""),
    "input in prepare_sipp() in sipp.r")
  if (!file.exists(input)) stop("Missing December SIPP input: ",input)
  raw <- readRDS(input)
  model_columns <- c("A","Y",sipp_covars)
  required <- c(model_columns,"SSUID","PNUM","person_id","MONTHCODE")
  if (!is.data.frame(raw) || !all(required %in% names(raw)))
    stop("Expected the identifier-preserving sipp_data_dec.rds.")
  if (nrow(raw)<50L || anyNA(raw[,required]) ||
      !all(raw$MONTHCODE==12L) || anyDuplicated(raw$person_id) ||
      !all(as.character(raw$person_id)==paste(raw$SSUID,raw$PNUM,sep=":")))
    stop("Input must contain one December record per person with valid identifiers.")
  # Keep all rows in their saved order. No resampling, recoding, or additional
  # standardization here. Identifiers and original-unit columns are not predictors.
  dat <- raw[,model_columns,drop=FALSE]
  if (!all(vapply(dat,is.numeric,logical(1))) || any(!is.finite(as.matrix(dat))) ||
      !setequal(dat$A,0:1) || !setequal(dat$Y,1:5)) stop("Invalid SIPP analysis variables.")
  rownames(dat) <- NULL
  if (nrow(dat) %% cfg$K!=0L)
    stop("Equal thirds require a sample size divisible by three; no rows are discarded.")
  attr(dat,"preprocessing") <- list(input_md5=unname(tools::md5sum(input)),
    N=nrow(dat),preprocessing="december-all-distinct-persons-v1",
    source=attr(raw,"preprocessing"))
  dat
}

# Randomize within each treatment-by-outcome cell and allocate cyclically.
# Carrying the offset across cells keeps total fold sizes within one row,
# while each cell's counts also differ by at most one across folds.
sipp_folds <- function(dat,K=3L) {
  stopifnot(nrow(dat)>=K,K>=2L,K==as.integer(K))
  cells <- split(seq_len(nrow(dat)),interaction(dat$A,dat$Y,drop=TRUE))
  fold <- integer(nrow(dat))
  fold_order <- sample.int(K)
  offset <- 0L
  for (idx in cells) {
    idx <- idx[sample.int(length(idx))]
    fold[idx] <- fold_order[((offset+seq_along(idx)-1L) %% K)+1L]
    offset <- (offset+length(idx)) %% K
  }
  stopifnot(all(fold %in% seq_len(K)),diff(range(tabulate(fold,nbins=K)))<=1L)
  fold
}

# Split distinct input people first, then draw the same number of rows from
# each fold with replacement. Returned indices retain the original person
# mapping, so copies can be checked for separation across all three roles.
sipp_bootstrap_split <- function(dat,K=3L) {
  stopifnot(nrow(dat) %% K==0L)
  original_fold <- sipp_folds(dat,K)
  index <- unlist(lapply(seq_len(K),function(k) {
    idx <- which(original_fold==k)
    idx[sample.int(length(idx),length(idx),replace=TRUE)]
  }),use.names=FALSE)
  fold <- original_fold[index]
  stopifnot(length(index)==nrow(dat),all(tabulate(fold,nbins=K)==nrow(dat)/K))
  list(index=index,fold=fold)
}
fit_sipp_nuisance <- function(dat) {
  if (!setequal(dat$A,0:1) || !setequal(dat$Y,1:5))
    stop("Nuisance fold lacks an arm/category; increase sample size.")
  fml <- ~ age + I(age^2) + female + race_black + race_asian + hispanic +
    educ_cat + povratio + married + insured_any + work_limited
  ps <- glm(update(fml,A~.),data=dat,family=binomial())
  dat$Y <- factor(dat$Y,levels=1:5)
  mu <- nnet::multinom(update(fml,Y~A+.),data=dat,trace=FALSE,maxit=1000)
  if (!isTRUE(ps$converged) || mu$convergence!=0 ||
      any(!is.finite(coef(ps))) || any(!is.finite(coef(mu))))
    stop("Nuisance fit did not converge or has nonfinite coefficients.")
  list(ps=ps,mu=mu)
}
predict_sipp_nuisance <- function(fit,dat,clip=.05) {
  d1 <- d0 <- dat; d1$A <- 1; d0$A <- 0
  m1 <- as.matrix(predict(fit$mu,d1,type="probs"))[,as.character(1:5),drop=FALSE]
  m0 <- as.matrix(predict(fit$mu,d0,type="probs"))[,as.character(1:5),drop=FALSE]
  e_raw <- predict(fit$ps,dat,type="response")
  e <- pmin(pmax(e_raw,clip),1-clip)
  if (any(!is.finite(c(m1,m0,e)))) stop("Nonfinite nuisance predictions.")
  list(m1=m1,m0=m0,e=e,e_raw=e_raw)
}
## Uses the score formulas in compute_scores() with fitted nuisance predictions.
sipp_scores <- function(dat,nuisance,Cu,beta_scale=2,beta_power=.25) {
  J <- 5L; A <- dat$A; Y <- dat$Y; n <- nrow(dat)
  m1_hat <- nuisance$m1; m0_hat <- nuisance$m0; e_hat <- nuisance$e
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
  n_score <- nrow(dat)
  beta <- beta_scale * n_score^beta_power
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

# A=0 has zero reduced objective. A=1 is evaluated with exactly the same
# held-out score as the corresponding learned policy; lower is better.
sipp_metrics <- function(cfg) {
  old <- c("Rsup_estimate","treat_rate")
  if (identical(cfg$version,"0903-sipp-december-v2-thirds")) old else
    c(old,"constant0_Rsup_estimate","constant1_Rsup_estimate")
}
sipp_result_columns <- function(cfg) c("N","J","beta_scale","beta_power","rep","Cu",
  unlist(lapply(cfg$estimators,function(e) paste(e,sipp_metrics(cfg),sep=".")),
    use.names=FALSE),sipp_diagnostic_columns)

# Called only after complete numeric results and the tree archive are verified.
# Resolve and check the exact child directory before any recursive removal.
remove_sipp_tree_checkpoints <- function(out) {
  parent <- normalizePath(out,winslash="/",mustWork=TRUE)
  expected <- paste0(sub("/+$","",parent),"/sipp_tree_checkpoints")
  if (!dir.exists(expected)) return(invisible(TRUE))
  target <- normalizePath(expected,winslash="/",mustWork=TRUE)
  same_path <- if (.Platform$OS.type=="windows")
    identical(tolower(target),tolower(expected)) else identical(target,expected)
  if (!same_path) stop("Refusing checkpoint cleanup outside the expected directory: ",target)
  status <- unlink(target,recursive=TRUE)
  if (status!=0L || dir.exists(target)) {
    warning("Final outputs are verified, but the checkpoint folder could not be fully removed: ",target)
    return(invisible(FALSE))
  }
  message("REMOVED completed-run tree checkpoints: ",target)
  invisible(TRUE)
}

# Preserve each fitted policytree object, but do not duplicate all-person
# predictions or nuisance fits. These leaf counts distinguish bootstrap rows
# from distinct people. The saved split maps each bootstrap row to the input.
sipp_tree_record <- function(tree,fold_data,split,k,score,pi) {
  train_node <- as.integer(predict(tree,fold_data$xt,type="node.id"))
  test_node <- as.integer(predict(tree,fold_data$xe,type="node.id"))
  train_person <- split$index[split$fold==fold_data$training_fold]
  test_person <- split$index[split$fold==fold_data$evaluation_fold]
  ids <- which(vapply(tree$nodes,function(node) isTRUE(node$is_leaf),logical(1)))
  leaves <- do.call(rbind,lapply(ids,function(id) data.frame(node_id=id,
    action=tree$nodes[[id]]$action-1L,
    training_n=sum(train_node==id),evaluation_n=sum(test_node==id),
    training_people=length(unique(train_person[train_node==id])),
    evaluation_people=length(unique(test_person[test_node==id])))))
  list(tree=tree,leaf_summary=leaves,evaluation=list(N=length(pi),
    treat_rate=mean(pi==1),reduced_Rsup_estimate=-mean(pi*score),
    constant0_Rsup_estimate=0,constant1_Rsup_estimate=-mean(score)))
}

validate_sipp_tree_cost <- function(cost,cfg) {
  if (!is.list(cost) || length(cost$Cu)!=1L || !cost$Cu %in% cfg$Cu ||
      length(cost$rotations)!=cfg$K) stop("Invalid saved tree cost/rotations.")
  for (k in seq_len(cfg$K)) {
    rotation <- cost$rotations[[k]]
    left <- ((k-1L+seq_len(cfg$K-1L)) %% cfg$K)+1L
    if (!identical(rotation$folds,c(nuisance=k,training=left[1],evaluation=left[2])) ||
        !identical(names(rotation$estimators),cfg$estimators))
      stop("Invalid saved tree roles/estimators.")
    for (item in rotation$estimators) {
      if (!inherits(item$tree,"policy_tree") ||
          !identical(item$tree$columns,sipp_covars) ||
          sum(item$leaf_summary$training_n)!=cfg$N/cfg$K ||
          sum(item$leaf_summary$evaluation_n)!=cfg$N/cfg$K ||
          any(item$leaf_summary$training_people>item$leaf_summary$training_n) ||
          any(item$leaf_summary$evaluation_people>item$leaf_summary$evaluation_n) ||
          item$evaluation$N!=cfg$N/cfg$K ||
          item$evaluation$constant0_Rsup_estimate!=0 ||
          any(!is.finite(unlist(item$evaluation)))) stop("Invalid saved policy tree.")
    }
  }
  invisible(TRUE)
}

validate_sipp_tree_replication <- function(part,cfg,id,cost_indices) {
  if (!identical(part$run_config,cfg) || part$rep!=id ||
      length(part$costs)!=length(cfg$Cu) || length(part$split$index)!=cfg$N ||
      length(part$split$fold)!=cfg$N || anyNA(part$split$index) ||
      any(part$split$index<1L | part$split$index>cfg$N) ||
      any(!part$split$fold %in% seq_len(cfg$K)) ||
      any(tabulate(part$split$fold,nbins=cfg$K)!=cfg$N/cfg$K) ||
      any(vapply(split(part$split$fold,part$split$index),function(x)
        length(unique(x))!=1L,logical(1)))) stop("Invalid tree replication/split.")
  for (j in cost_indices) {
    validate_sipp_tree_cost(part$costs[[j]],cfg)
    if (part$costs[[j]]$Cu!=cfg$Cu[j]) stop("Misindexed saved tree cost.")
  }
  invisible(TRUE)
}

sipp_replication <- function(rep_id,dat,cfg,cost_indices=seq_along(cfg$Cu),
                             on_cost=NULL) {
  ## Recreate the same split, within-fold bootstrap, and three nuisance fits
  ## when resuming. Partitions are redrawn reproducibly for each replication.
  ## Saved costs are skipped only after this deterministic preparation.
  RNGkind("L'Ecuyer-CMRG"); set.seed(cfg$seed + rep_id)
  split <- sipp_bootstrap_split(dat,cfg$K)
  boot <- dat[split$index,,drop=FALSE]
  fold <- split$fold
  cols <- unlist(lapply(cfg$estimators,function(e)
    paste(e,sipp_metrics(cfg),sep=".")),use.names=FALSE)
  prepared <- vector("list",cfg$K)
  weights <- numeric(cfg$K)
  diagnostics <- matrix(NA_real_,cfg$K,length(sipp_diagnostic_columns),
    dimnames=list(NULL,sipp_diagnostic_columns))
  for (k in seq_len(cfg$K)) {
    left <- ((k-1L + seq_len(cfg$K-1L)) %% cfg$K) + 1L
    nui <- boot[fold==k,,drop=FALSE]
    train <- boot[fold==left[1],,drop=FALSE]
    test <- boot[fold==left[2],,drop=FALSE]
    stopifnot(nrow(nui)==nrow(train),nrow(train)==nrow(test))
    # Equal training/evaluation sizes also give the same beta in both scores.
    weights[k] <- nrow(test)
    fit <- fit_sipp_nuisance(nui)
    prepared[[k]] <- list(train=train,test=test,
      training_fold=left[1],evaluation_fold=left[2],
      nt=predict_sipp_nuisance(fit,train,cfg$clip),
      ne=predict_sipp_nuisance(fit,test,cfg$clip),
      xt=as.matrix(train[,sipp_covars]),xe=as.matrix(test[,sipp_covars]))
    # Out-of-fold test predictions before clipping, pooled across rotations.
    e_raw <- prepared[[k]]$ne$e_raw
    diagnostics[k,] <- c(min(e_raw),max(e_raw),mean(e_raw<cfg$clip),
      mean(e_raw>1-cfg$clip))
  }
  diagnostic_summary <- c(ps_min=min(diagnostics[,"ps_min"]),
    ps_max=max(diagnostics[,"ps_max"]),
    ps_clipped_low=weighted.mean(diagnostics[,"ps_clipped_low"],weights),
    ps_clipped_high=weighted.mean(diagnostics[,"ps_clipped_high"],weights))
  packets <- vector("list",length(cost_indices))
  for (position in seq_along(cost_indices)) {
    j <- cost_indices[position]
    ## Cost-specific stream also protects resume if a tree backend uses RNG.
    set.seed(cfg$seed + 100000L + (rep_id-1L)*length(cfg$Cu) + j)
    values <- matrix(NA_real_,cfg$K,length(cols))
    rotations <- vector("list",cfg$K)
    for (k in seq_len(cfg$K)) {
      fold_data <- prepared[[k]]
      st <- sipp_scores(fold_data$train,fold_data$nt,cfg$Cu[j],cfg$beta_scale,cfg$beta_power)
      se <- sipp_scores(fold_data$test,fold_data$ne,cfg$Cu[j],cfg$beta_scale,cfg$beta_power)
      records <- setNames(lapply(cfg$estimators,function(e) {
        tree <- policytree::policy_tree(fold_data$xt,cbind(0,st[[e]]),depth=cfg$depth)
        pi <- predict(tree,fold_data$xe)-1
        sipp_tree_record(tree,fold_data,split,k,se[[e]],pi)
      }),cfg$estimators)
      values[k,] <- unlist(lapply(records,function(item) with(item$evaluation,
        c(reduced_Rsup_estimate,treat_rate,constant0_Rsup_estimate,
          constant1_Rsup_estimate))),use.names=FALSE)
      rotations[[k]] <- list(folds=c(nuisance=k,training=fold_data$training_fold,
        evaluation=fold_data$evaluation_fold),estimators=records)
    }
    row <- data.frame(N=nrow(dat),J=cfg$J,
      beta_scale=cfg$beta_scale,beta_power=cfg$beta_power,rep=rep_id,Cu=cfg$Cu[j])
    for (i in seq_along(cols)) row[[cols[i]]] <- weighted.mean(values[,i],weights)
    # Repeated across costs within a replication because nuisances are shared.
    for (nm in sipp_diagnostic_columns) row[[nm]] <- unname(diagnostic_summary[nm])
    packet <- list(row=row,cost=list(Cu=cfg$Cu[j],rotations=rotations))
    packets[[position]] <- packet
    # Live delivery includes the common split for checkpointing. The normal
    # completion return stores it once per replication, not once per tree.
    if (!is.null(on_cost)) on_cost(c(packet,list(split=split)))
  }
  list(rep=rep_id,split=split,packets=packets)
}
validate_sipp <- function(results,cfg,complete=TRUE) {
  stopifnot(cfg$version %in% c("0903-sipp-december-v2-thirds","0903-sipp-december-v3-all-trees"),
    length(cfg$N)==1L,is.finite(cfg$N),cfg$N>=50,cfg$N==floor(cfg$N))
  beta_columns <- c("beta_scale","beta_power")
  columns <- sipp_result_columns(cfg)
  if (!is.data.frame(results) || !identical(names(results),columns) ||
      !all(vapply(results,function(x) is.numeric(x) && all(is.finite(x)),logical(1))))
    stop("Invalid SIPP result columns or values.")
  if (any(results$rep!=floor(results$rep) | results$rep<1 | results$rep>cfg$R) ||
      any(results$N!=cfg$N | results$J!=cfg$J))
    stop("Invalid SIPP repetitions or settings.")
  for (field in beta_columns) if (any(results[[field]]!=cfg[[field]]))
    stop("Invalid smoothing configuration in result rows.")
  cost_index <- match(results$Cu,cfg$Cu)
  if (anyNA(cost_index) || anyDuplicated(data.frame(rep=results$rep,cost=cost_index)))
    stop("Unknown or duplicate repetition-Cu pairs.")
  if (complete && nrow(results)!=cfg$R*length(cfg$Cu))
    stop("Incomplete SIPP results: wait for every repetition-Cu pair.")
  if (any(vapply(results[grep("[.]treat_rate$",names(results))],
    function(x) any(x<0 | x>1),logical(1)))) stop("Invalid treatment proportion.")
  if (any(vapply(results[grep("[.]constant0_Rsup_estimate$",names(results))],
    function(x) any(x!=0),logical(1)))) stop("Invalid constant-zero benchmark.")
  if (any(vapply(results[sipp_diagnostic_columns],
      function(x) any(x<0 | x>1),logical(1))) || any(results$ps_min>results$ps_max) ||
      any(results$ps_clipped_low+results$ps_clipped_high>1+1e-12))
    stop("Invalid propensity diagnostics.")
  invisible(TRUE)
}
## Near-live packets are handled only by the main process; workers never
## write the shared RDS. Completion copies provide delivery verification.
sipp_progress_handler <- function(on_packet) {
  reporter <- list(update=function(config,state,progression,...) {
    if (!is.null(progression$sipp_packet)) on_packet(progression$sipp_packet)
  })
  progressr::make_progression_handler("sipp_costs",reporter,
    interval=0,times=Inf,intrusiveness=0,clear=FALSE,enable=TRUE)
}
run_sipp <- function(resume=FALSE) {
  stopifnot(is.logical(resume),length(resume)==1L,!is.na(resume))
  cfg <- sipp_settings()
  workers <- sipp_parameters$workers
  save_every <- sipp_parameters$save_every
  progress_every <- sipp_parameters$progress_every
  stopifnot(length(workers)==1L,is.finite(workers),workers>=1,workers<=56,workers==floor(workers),
    length(save_every)==1L,is.finite(save_every),save_every>=1,save_every==floor(save_every),
    length(progress_every)==1L,is.finite(progress_every),progress_every>=1,
    progress_every==floor(progress_every))
  out <- sipp_output_dir(); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  dat <- prepare_sipp(cfg)
  cfg$N <- nrow(dat)
  cfg$data <- attr(dat,"preprocessing")
  cfg$policytree_version <- as.character(utils::packageVersion("policytree"))
  sources <- c("sipp.r","sim_Rsup.R","data_gen.R")
  source_paths <- c(file.path(sipp_root,"sipp.r"),
    file.path(sipp_helper_root,c("sim_Rsup.R","data_gen.R")))
  cfg$source_md5 <- setNames(unname(tools::md5sum(source_paths)),sources)
  if (anyNA(cfg$source_md5)) stop("Missing SIPP source/helper file.")
  file <- file.path(out,"sipp_res.rds"); recover_rds_backup(file)
  columns <- sipp_result_columns(cfg)
  results <- as.data.frame(setNames(rep(list(numeric()),length(columns)),columns))
  if (resume && file.exists(file)) {
    saved <- readRDS(file)
    if (!identical(attr(saved,"run_config"),cfg))
      stop("Cannot resume: saved settings/source/data differ. Run without --resume to start fresh.")
    validate_sipp(saved,cfg,complete=FALSE)
    results <- plain_result(saved)
  }
  total <- cfg$R*length(cfg$Cu)
  completed <- matrix(FALSE,cfg$R,length(cfg$Cu))
  if (nrow(results)) completed[cbind(results$rep,match(results$Cu,cfg$Cu))] <- TRUE
  tree_file <- file.path(out,"sipp_bootstrap_trees.rds")
  checkpoint_file <- function(id) file.path(out,"sipp_tree_checkpoints",
    sprintf("rep_%03d.rds",id))
  tree_parts <- vector("list",cfg$R)
  dirty_trees <- rep(FALSE,cfg$R)
  # A complete archive is sufficient for complete-run resume; checkpoint
  # files are needed only for an interrupted run. Extra costs written before
  # a result-checkpoint commit are ignored and deterministically recomputed.
  archive <- NULL
  if (resume && all(completed)) {
    recover_rds_backup(tree_file)
    if (file.exists(tree_file)) {
      candidate <- readRDS(tree_file)
      if (identical(candidate$run_config,cfg) &&
          length(candidate$replications)==cfg$R) archive <- candidate
    }
  }
  for (id in which(rowSums(completed)>0L)) {
    if (!is.null(archive)) {
      part <- archive$replications[[id]]; part$run_config <- cfg
    } else {
      path <- checkpoint_file(id); recover_rds_backup(path)
      if (!file.exists(path)) stop("Missing tree checkpoint: ",path,
        ". Restore it or start a fresh run without --resume.")
      part <- readRDS(path)
    }
    validate_sipp_tree_replication(part,cfg,id,which(completed[id,]))
    part$costs[which(!completed[id,])] <- rep(list(NULL),sum(!completed[id,]))
    tree_parts[[id]] <- part
  }
  pending <- which(rowSums(completed)<length(cfg$Cu))
  snapshot <- function() {
    x <- results[order(results$rep,results$Cu),,drop=FALSE]
    rownames(x) <- NULL; attr(x,"run_config") <- cfg
    x
  }
  finish_trees <- function() {
    for (id in seq_len(cfg$R))
      validate_sipp_tree_replication(tree_parts[[id]],cfg,id,seq_along(cfg$Cu))
    saved <- list(format="sipp-bootstrap-policy-trees-v1",run_config=cfg,
      covariates=sipp_covars,preprocessing=attr(dat,"preprocessing"),
      bootstrap=TRUE,tree_count=cfg$R*cfg$K*length(cfg$Cu)*length(cfg$estimators),
      replications=lapply(tree_parts,function(part) {part$run_config <- NULL; part}),
      notes=c("Index by replications[[rep]]$costs[[cost_index]]$rotations[[rotation]]$estimators[[estimator]].",
        "Each estimator entry contains tree, leaf_summary and held-out evaluation, including both constant policies.",
        "split$index maps bootstrap rows to sipp_data_dec.rds input rows; split$fold gives their fold.",
        "Leaf *_n counts bootstrap rows; *_people counts distinct input respondents.",
        "Use predict(tree, X)-1 for actions 0/1. All-person predictions can be reconstructed from the input data.",
        "Cutoffs use the analysis scale; preprocessing$source$scaling supplies original units.",
        "Rotations within a replication are dependent. No separate illustrative fits are made."))
    if (!identical(archive,saved)) atomic_rds_save(saved,tree_file)
    final_results <- snapshot()
    validate_sipp(final_results,cfg)
    stopifnot(identical(readRDS(file),final_results),
      identical(readRDS(tree_file),saved))
    message("SAVED/VERIFIED ",saved$tree_count," bootstrap trees: ",tree_file)
    remove_sipp_tree_checkpoints(out)
    invisible(saved)
  }
  if (!length(pending)) {
    finish_trees()
    message("REUSED complete SIPP result: ",file)
    return(invisible(snapshot()))
  }
  workers <- min(workers,length(pending))
  unsaved <- 0L
  last_progress <- nrow(results)
  started <- Sys.time()
  log_line <- function(message) {
    cat(sprintf("[%s] %s\n",format(Sys.time(),"%H:%M:%S"),message))
    flush.console()
  }
  persist <- function() {
    x <- snapshot(); validate_sipp(x,cfg,complete=FALSE)
    # Commit policies before the numeric rows that claim they are complete.
    # Only changed replication files are rewritten; the full archive is
    # assembled once at completion. All disk writes occur in this process.
    for (id in which(dirty_trees)) atomic_rds_save(tree_parts[[id]],checkpoint_file(id))
    atomic_rds_save(x,file)
    dirty_trees[] <<- FALSE
    unsaved <<- 0L
  }
  ## On a handled R error/interrupt, try to save already delivered results.
  ## A forced process termination cannot run this cleanup.
  on.exit(if (unsaved>0L) tryCatch(persist(),error=function(e)
    warning("Could not flush pending SIPP results: ",conditionMessage(e))),add=TRUE)
  accept_packet <- function(packet) {
    row <- packet$row
    validate_sipp(row,cfg,complete=FALSE)
    if (nrow(row)!=1L) stop("Expected one completed repetition-Cu pair.")
    id <- as.integer(row$rep); j <- match(row$Cu,cfg$Cu)
    validate_sipp_tree_cost(packet$cost,cfg)
    if (packet$cost$Cu!=row$Cu) stop("Tree cost and numeric result disagree.")
    for (e in cfg$estimators) {
      expected <- row[1,paste(e,sipp_metrics(cfg),sep=".")]
      actual <- colMeans(do.call(rbind,lapply(packet$cost$rotations,function(rotation)
        with(rotation$estimators[[e]]$evaluation,
          c(reduced_Rsup_estimate,treat_rate,constant0_Rsup_estimate,constant1_Rsup_estimate)))))
      if (any(abs(as.numeric(expected)-actual)>1e-12))
        stop("Saved tree evaluation and numeric result disagree.")
    }
    if (completed[id,j]) {
      prior <- plain_result(results[results$rep==id & results$Cu==cfg$Cu[j],,drop=FALSE])
      if (!identical(prior,plain_result(row)) ||
          !identical(tree_parts[[id]]$costs[[j]],packet$cost) ||
          !identical(tree_parts[[id]]$split,packet$split))
        stop("Conflicting duplicate SIPP result/tree.")
      return(invisible(NULL))
    }
    if (is.null(tree_parts[[id]])) tree_parts[[id]] <<- list(run_config=cfg,
      rep=id,split=packet$split,costs=vector("list",length(cfg$Cu)))
    if (!identical(tree_parts[[id]]$split,packet$split)) stop("Conflicting bootstrap split.")
    tree_parts[[id]]$costs[[j]] <<- packet$cost
    dirty_trees[id] <<- TRUE
    results <<- plain_result(rbind(results,row))
    completed[id,j] <<- TRUE
    unsaved <<- unsaved+1L
    ## Report total work independently so a long run does not stay silent
    ## while waiting for an entire cost to finish across all repetitions.
    if (nrow(results)-last_progress>=progress_every || nrow(results)==total) {
      log_line(sprintf("TOTAL PROGRESS: %d/%d pairs (%.1f%%) | Complete repetitions: %d/%d | Elapsed: %.1f min.",
        nrow(results),total,100*nrow(results)/total,
        sum(rowSums(completed)==length(cfg$Cu)),cfg$R,
        as.numeric(difftime(Sys.time(),started,units="mins"))))
      last_progress <<- nrow(results)
    }
    ## Print only once per cost, after every repetition has completed it.
    ## Checkpoint writes continue silently at the fixed save interval.
    if (all(completed[,j])) log_line(sprintf(
      "Cu=%.3f COMPLETE across all %d repetitions | Costs: %d/%d | Total work: %d/%d (%.1f%%) | Elapsed: %.1f min.",
      cfg$Cu[j],cfg$R,sum(colSums(completed)==cfg$R),length(cfg$Cu),
      nrow(results),total,100*nrow(results)/total,
      as.numeric(difftime(Sys.time(),started,units="mins"))))
    if (unsaved>=save_every || nrow(results)==total) persist()
    invisible(NULL)
  }
  resume_completed <- completed  # Frozen worker snapshot.
  previous_plan <- future::plan(); on.exit(future::plan(previous_plan),add=TRUE)
  doFuture::registerDoFuture()
  if (workers==1L) future::plan(future::sequential) else
    future::plan(future::multisession,workers=workers)
  if (!resume) {
    log_line("NEW RUN: replacing sipp_res.rds with a fresh checkpoint.")
    persist()
  }
  log_line(sprintf("SIPP: %d workers; %d/%d pairs already saved; save every %d pairs; total progress every %d pairs.",
    workers,nrow(results),total,save_every,progress_every))
  finished <- progressr::with_progress({
    progress <- progressr::progressor(steps=total-nrow(results))
    foreach(id=pending,.packages=c("nnet","policytree"),
      .options.future=list(chunk.size=1L),.options.RNG=cfg$seed) %dorng% {
      sipp_replication(id,dat,cfg,cost_indices=which(!resume_completed[id,]),
        on_cost=function(packet) progress(sipp_packet=packet))
    }
  },handlers=sipp_progress_handler(accept_packet),enable=TRUE,
    delay_stdout=FALSE,delay_conditions=character())
  ## Deduplicate the normal future return against live deliveries.
  for (part in finished) for (packet in part$packets)
    accept_packet(c(packet,list(split=part$split)))
  if (unsaved>0L) persist()
  answer <- snapshot(); validate_sipp(answer,cfg)
  stopifnot(identical(readRDS(file),answer))
  finish_trees()
  log_line(sprintf("FINISHED: verified %s; elapsed %.1f minutes.",file,
    as.numeric(difftime(Sys.time(),started,units="mins"))))
  invisible(answer)
}
if (!isTRUE(getOption("opl.sipp.functions_only",FALSE)))
  run_sipp(resume="--resume" %in% commandArgs(trailingOnly=TRUE))
