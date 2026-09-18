## Four-estimator December SIPP summaries. USER PATHS: see comments below.
## Outputs: sipp_trP.pdf, sipp_Rsup.pdf, sipp_Rsup.xlsx, sipp_summary.csv.
## To plot saved trees, choose a bootstrap replication and rotation, then call
## run_sipp_tree_plot(); no additional policy fitting is needed.
## Bootstrap SDs describe variation across redrawn partitions and within-fold
## resamples. This script loads sipp.r with model execution disabled.
sipp_plot_file <- local({
  files <- Filter(Negate(is.null),lapply(sys.frames(),function(x) x$ofile))
  arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
  if (length(files)) tail(files,1)[[1]] else if (length(arg))
    sub("^--file=","",arg[[1]]) else stop("Use source() or Rscript.")
})
local({
  saved <- options(opl.sipp.functions_only=TRUE)
  on.exit(options(saved))
  source(file.path(dirname(normalizePath(sipp_plot_file)),"sipp.r"),local=.GlobalEnv)
})
suppressPackageStartupMessages(library(ggplot2))
## USER PATH: replace the empty string with your absolute figure/table directory.
## Input results use sipp_output_dir() in sipp.r. Use forward slashes on Windows.
sipp_figure_dir <- function() require_user_path(getOption("opl.sipp.figure_dir", ""),
  "the output directory in sipp_figure_dir() in sipp_plot.r")
sipp_tree_plot_parameters <- list(Cu=c(.26,.29),width=10,height=3.8,
  estimator="orthogonal_smoothed",rep=1L,rotation=1L)
build_sipp_summary <- function(results) {
  cfg <- attr(results,"run_config")
  if (is.null(cfg) || !cfg$version %in%
      c("0903-sipp-december-v2-thirds","0903-sipp-december-v3-all-trees") ||
      !identical(cfg$estimators,names(sipp_labels))) stop("Unknown SIPP result configuration.")
  validate_sipp(results,cfg)
  do.call(rbind,lapply(cfg$Cu,function(cost) {
    x <- results[results$Cu==cost,]
    do.call(rbind,lapply(cfg$estimators,function(e) {
      risk <- x[[paste0(e,".Rsup_estimate")]]
      treatment <- x[[paste0(e,".treat_rate")]]
      data.frame(N=cfg$N,J=cfg$J,Cu=cost,estimator=e,R=cfg$R,
        Rsup_estimate=mean(risk),Rsup_boot_sd=sd(risk),
        treat_rate=mean(treatment),treat_rate_boot_sd=sd(treatment))
    }))
  }))
}
plot_sipp <- function(summary,metric) {
  stopifnot(metric %in% c("treat_rate","Rsup_estimate"))
  ## Broad panel: .05 grid. Detail panel: .01 grid, reusing the same estimates.
  full_grid <- seq(10L,60L,by=5L)/100
  full <- summary[vapply(summary$Cu,function(cost)
    any(abs(cost-full_grid)<1e-8),logical(1)),]
  full$view <- rep("Full range",nrow(full))
  zoom <- summary[summary$Cu>=.25-1e-8 & summary$Cu<=.35+1e-8,]
  zoom$view <- rep("Transition region",nrow(zoom))
  dat <- rbind(full,zoom)
  dat$view <- factor(dat$view,c("Full range","Transition region"))
  dat$estimator <- factor(dat$estimator,names(sipp_labels))
  dat$value <- dat[[metric]]
  line_dat <- dat[ave(dat$Cu,dat$view,dat$estimator,FUN=length)>1,]
  p <- ggplot(dat,aes(Cu,value,color=estimator,linetype=estimator,shape=estimator)) +
    geom_line(data=line_dat) + geom_point(size=2.15,stroke=.35) +
    facet_wrap(~view,nrow=1,scales="free_x") +
    scale_x_continuous(
      breaks=function(limits) if (diff(limits)<.2) seq(25L,35L,by=2L)/100 else seq(0,1,by=.1),
      labels=function(x) if (any(abs(x*10-round(x*10))>1e-8,na.rm=TRUE))
        sprintf("%.2f",x) else sprintf("%.1f",x)) +
    scale_color_manual(values=c(plugin="#0072B2",direct_if="#8E44AD",
      smoothed_plugin="#009E73",orthogonal_smoothed="#D62728"),labels=sipp_labels) +
    scale_linetype_manual(values=c(plugin="dashed",direct_if="dotdash",
      smoothed_plugin="dotted",orthogonal_smoothed="solid"),labels=sipp_labels) +
    scale_shape_manual(values=c(plugin=16,direct_if=18,smoothed_plugin=15,
      orthogonal_smoothed=17),labels=sipp_labels) +
    labs(x=expression(C[u]),y=if(metric=="treat_rate") "Treated proportion" else
      expression(hat(R)[sup]),color=NULL,linetype=NULL,shape=NULL) +
    theme_bw(base_size=13) + theme(
      text=element_text(size=13),
      axis.text=element_text(size=13),axis.title=element_text(size=13),
      strip.text=element_text(size=13),plot.title=element_text(size=13),
      plot.subtitle=element_text(size=13),plot.caption=element_text(size=13),
      plot.tag=element_text(size=13),
      legend.position="inside",legend.position.inside=c(1,1),
      legend.justification.inside=c(1,1),legend.direction="vertical",
      legend.background=element_blank(),
      legend.box.background=element_blank(),legend.key=element_blank(),
      legend.key.height=grid::unit(10,"pt"),
      legend.key.spacing.y=grid::unit(1,"pt"),
      legend.box.margin=margin(0,0,0,0),legend.margin=margin(4,4,4,4),
      legend.title=element_text(size=10),legend.text=element_text(size=10)) +
    guides(color=guide_legend(ncol=1),linetype=guide_legend(ncol=1),
      shape=guide_legend(ncol=1))
  if (metric=="treat_rate") p <- p + coord_cartesian(ylim=c(0,1))
  ## Attach the inside legend to the exact first panel cell. Using a fraction
  ## of the full plotting width would misalign it because of facet spacing.
  grob <- ggplotGrob(p)
  panel <- which(grob$layout$name=="panel-1-1")
  legend <- which(grob$layout$name=="guide-box-inside")
  stopifnot(length(panel)==1L,length(legend)==1L)
  grob$layout[legend,c("t","l","b","r")] <-
    grob$layout[panel,c("t","l","b","r")]
  grob
}
# Convert internal standardized cutoffs back to the dataset's original units.
# Binary nodes use No/Yes branches; numerical comparisons use Yes/No.
sipp_tree_split_label <- function(node,artifact) {
  variable <- artifact$covariates[node$split_variable]
  cutoff <- node$split_value
  binary <- c(female="Female",race_black="Black race",race_asian="Asian race",
    hispanic="Hispanic",married="Married",insured_any="Health insurance\ncoverage",
    work_limited="Work limitation")
  if (variable %in% names(binary) && cutoff>=0 && cutoff<1)
    return(list(label=paste0(binary[[variable]],"?"),left="No",right="Yes"))
  if (variable %in% c("age","povratio")) {
    scaling <- artifact$preprocessing$source$scaling
    j <- match(variable,scaling$variable)
    if (is.na(j)) stop("Missing original-unit scaling for ",variable)
    cutoff <- scaling$center[j]+scaling$scale[j]*cutoff
  }
  if (abs(cutoff-round(cutoff))<1e-8) cutoff <- round(cutoff)
  value <- format(signif(cutoff,6),trim=TRUE,scientific=FALSE)
  title <- switch(variable,age="Age (years)",povratio="Income-to-poverty ratio",
    educ_cat="Education level",variable)
  if (variable=="educ_cat" && cutoff>=1 && cutoff<4)
    value <- c("High school or less","Some college","Bachelor's degree")[floor(cutoff)]
  list(label=paste0(title,"\n\u2264 ",value,"?"),left="Yes",right="No")
}

sipp_tree_panel <- function(item,artifact,panel) {
  nodes <- item$tree$nodes
  positions <- list(); edges <- list()
  walk <- function(id,lo,hi,depth) {
    node <- nodes[[id]]
    leaf <- isTRUE(node$is_leaf)
    x <- (lo+hi)/2
    y <- if (leaf) {if (depth==0L) .50 else .18} else .79-.28*depth
    label <- if (leaf) NULL else sipp_tree_split_label(node,artifact)
    positions[[as.character(id)]] <<- list(id=id,x=x,y=y,leaf=leaf,
      width=if (leaf) .21 else if (depth==0L) .54 else .43,height=.145,label=label)
    if (!leaf) {
      walk(node$left_child,lo,x,depth+1L)
      walk(node$right_child,x,hi,depth+1L)
      edges[[length(edges)+1L]] <<- list(from=id,to=node$left_child,label=label$left)
      edges[[length(edges)+1L]] <<- list(from=id,to=node$right_child,label=label$right)
    }
  }
  walk(1L,0,1,0L)
  grid::grid.text(bquote(.(paste0("(",letters[panel],")"))~C[u]==.(sprintf("%.2f",item$Cu))),
    x=.5,y=.96,gp=grid::gpar(fontfamily="serif",fontsize=13))
  for (edge in edges) {
    from <- positions[[as.character(edge$from)]]
    to <- positions[[as.character(edge$to)]]
    y1 <- from$y-from$height/2; y2 <- to$y+to$height/2
    grid::grid.lines(x=c(from$x,to$x),y=c(y1,y2),
      arrow=grid::arrow(length=grid::unit(1.4,"mm"),type="closed"),
      gp=grid::gpar(col="#626262",fill="#626262",lwd=.8))
    mx <- (from$x+to$x)/2; my <- (y1+y2)/2
    grid::grid.rect(mx,my,width=.08,height=.052,
      gp=grid::gpar(fill="white",col=NA))
    grid::grid.text(edge$label,mx,my,gp=grid::gpar(fontfamily="serif",fontsize=9))
  }
  for (position in positions) {
    if (position$leaf) {
      row <- item$leaf_summary[item$leaf_summary$node_id==position$id,,drop=FALSE]
      if (nrow(row)!=1L || row$action!=nodes[[position$id]]$action-1L)
        stop("Tree and saved leaf summary disagree.")
      # Plotmath gives A mathematical italics and keeps text/numerals upright.
      label <- bquote(atop("Assign"~italic(A)==.(row$action),
        .(sprintf("%.1f%%",100*row$full_sample_share))))
      fill <- "white"
    } else {
      label <- position$label$label; fill <- "white"
    }
    grid::grid.roundrect(position$x,position$y,width=position$width,height=position$height,
      r=grid::unit(1.5,"mm"),gp=grid::gpar(fill=fill,col="#505050",lwd=.8))
    grid::grid.text(label,position$x,position$y,
      gp=grid::gpar(fontfamily="serif",fontsize=if(position$leaf) 10 else 10.5,lineheight=1.15))
  }
}

# Load with options(opl.sipp.plot_functions_only=TRUE) to call this separately;
# it reads fitted tree objects and never refits a model.
run_sipp_tree_plot <- function(expected_config=NULL) {
  settings <- sipp_tree_plot_parameters
  if (length(settings$Cu)!=2L || anyDuplicated(settings$Cu))
    stop("Choose two distinct utility thresholds for the tree panels.")
  input <- file.path(sipp_output_dir(),"sipp_bootstrap_trees.rds")
  use_bootstrap <- file.exists(input)
  if (!use_bootstrap) {
    if (!is.null(expected_config) &&
        identical(expected_config$version,"0903-sipp-december-v3-all-trees"))
      stop("Missing bootstrap tree archive: ",input)
    input <- file.path(sipp_output_dir(),"sipp_policy_trees.rds")
  }
  if (!file.exists(input)) stop("Missing saved policy trees: ",input)
  artifact <- readRDS(input)
  if (!is.null(expected_config) && !identical(artifact$run_config,expected_config))
    stop("Policy trees and bootstrap results belong to different runs.")
  if (use_bootstrap) {
    cfg <- artifact$run_config
    if (!identical(artifact$format,"sipp-bootstrap-policy-trees-v1") ||
        !settings$estimator %in% cfg$estimators || length(settings$rep)!=1L ||
        !settings$rep %in% seq_len(cfg$R) || length(settings$rotation)!=1L ||
        !settings$rotation %in% seq_len(cfg$K)) stop("Invalid bootstrap tree selection.")
    dat <- prepare_sipp(cfg)
    if (!identical(attr(dat,"preprocessing"),artifact$preprocessing))
      stop("Tree archive and current December input differ.")
    x <- as.matrix(dat[,artifact$covariates])
    part <- artifact$replications[[settings$rep]]
    artifact$Cu <- cfg$Cu
    artifact$estimator <- settings$estimator
    artifact$trees <- lapply(part$costs,function(cost) {
      item <- cost$rotations[[settings$rotation]]$estimators[[settings$estimator]]
      item$Cu <- cost$Cu
      node <- as.integer(predict(item$tree,x,type="node.id"))
      item$leaf_summary$full_sample_share <- vapply(item$leaf_summary$node_id,
        function(id) mean(node==id),numeric(1))
      item
    })
    message("Plotting saved bootstrap trees: replication ",settings$rep,
      ", rotation ",settings$rotation,", estimator ",settings$estimator,".")
  }
  if (!artifact$estimator %in% names(sipp_labels) ||
      !identical(artifact$covariates,sipp_covars)) stop("Unexpected policy-tree specification.")
  selected <- vapply(settings$Cu,function(cost) {
    hit <- which(abs(artifact$Cu-cost)<1e-10)
    if (length(hit)!=1L) stop("No unique saved tree for Cu=",cost)
    hit
  },integer(1))
  out <- sipp_figure_dir()
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  file <- file.path(out,"sipp_policy_trees.pdf")
  grDevices::cairo_pdf(file,width=settings$width,height=settings$height,family="serif")
  on.exit(grDevices::dev.off(),add=TRUE)
  grid::grid.newpage()
  for (panel in 1:2) {
    grid::pushViewport(grid::viewport(x=c(.25,.75)[panel],y=.50,width=.48,height=.96))
    sipp_tree_panel(artifact$trees[[selected[panel]]],artifact,panel)
    grid::popViewport()
  }
  message("Two-panel policy-tree PDF written to: ",file)
  invisible(file)
}

run_sipp_plot <- function() {
  input <- file.path(sipp_output_dir(),"sipp_res.rds")
  if (!file.exists(input)) stop("Run sipp.r first: ",input)
  results <- readRDS(input)
  summary <- build_sipp_summary(results)
  out <- sipp_figure_dir()
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  if (!requireNamespace("openxlsx",quietly=TRUE)) stop("Install openxlsx for the table.")
  ggsave(file.path(out,"sipp_trP.pdf"),plot_sipp(summary,"treat_rate"),width=8,height=3)
  ggsave(file.path(out,"sipp_Rsup.pdf"),plot_sipp(summary,"Rsup_estimate"),width=8,height=3)
  write.csv(summary,file.path(out,"sipp_summary.csv"),row.names=FALSE)
  selected <- summary[vapply(summary$Cu,function(x)
    any(abs(x-seq(.25,.3,.01))<1e-8),logical(1)),]
  notes <- data.frame(note=c("All values are in original units.",
    "Rsup is the reduced objective -mean(policy * estimated score).",
    "Bootstrap SD is resampling variability; no true Rsup/oracle is observed.",
    "N is the total bootstrap sample size; each estimation role uses N/3 rows.",
    "People are partitioned by treatment and outcome, then resampled within folds; partitions are redrawn each replication."))
  openxlsx::write.xlsx(list(Rsup_selected=selected,All_costs=summary,Notes=notes),
    file.path(out,"sipp_Rsup.xlsx"),overwrite=TRUE)
  if (identical(attr(results,"run_config")$version,"0903-sipp-december-v2-thirds"))
    run_sipp_tree_plot(expected_config=attr(results,"run_config")) else
    message("All-tree run: summary plots only. Choose a saved replication/rotation and call run_sipp_tree_plot() for a tree figure.")
  message("SIPP figures, workbook and summary written to: ",out)
  invisible(summary)
}
if (!isTRUE(getOption("opl.sipp.plot_functions_only",FALSE))) run_sipp_plot()
