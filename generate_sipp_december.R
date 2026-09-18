# Generate one December observation per eligible person from the original ZIP.
# Run with Rscript generate_sipp_december.R, or source this file in R.
# Requires haven, tidyselect, and PowerShell 7 (for the Deflate64 ZIP archive).
# USER PATHS: fill in the commented raw input and prepared output paths below.
# Running this optional step replaces the configured prepared data file.
script_file <- local({
  sourced <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  if (length(sourced)) tail(sourced,1)[[1]] else if (length(arg))
    sub('^--file=','',arg[[1]]) else stop('Run this file with source() or Rscript.')
})
generate_sipp_december <- function() {
  # USER PATHS: replace each empty string with an absolute FILE path.
  # output_file is the prepared .rds filename; archive includes the ZIP filename.
  # Use forward slashes on Windows. R options may override these paths.
  output_file <- path.expand(getOption('opl.sipp.data_file', ''))
  archive <- path.expand(getOption('opl.sipp.raw_zip', ''))
  for (path in c(output_file,archive)) {
    if (length(path)!=1L || is.na(path) || !nzchar(path) ||
        !grepl('^(/|[A-Za-z]:/)',gsub('\\\\','/',path)))
      stop('Set output_file and archive in generate_sipp_december.R to absolute file paths.')
  }
  if (!file.exists(archive)) stop('Raw ZIP not found: ',archive,
    '. Set archive in generate_sipp_december.R, or use the included prepared data.')
  archive <- normalizePath(archive,winslash='/',mustWork=TRUE)
  if (identical(tolower(output_file),tolower(archive)))
    stop('Raw input and prepared output must use different paths.')
  out <- dirname(output_file)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  protected <- archive
  before <- tools::md5sum(protected)
  stopifnot(requireNamespace('haven',quietly=TRUE),
            requireNamespace('tidyselect',quietly=TRUE))
  entries <- utils::unzip(archive,list=TRUE)
  stopifnot(nrow(entries)==1L,entries$Name=='pu2022.dta')
  scratch <- tempfile('december_extract_',tmpdir=out)
  dir.create(scratch)
  scratch <- normalizePath(scratch,winslash='/',mustWork=TRUE)
  stopifnot(startsWith(tolower(scratch),paste0(tolower(out),'/')))
  # Cleanup is confined to this newly created, checked directory.
  on.exit({
    if (dir.exists(scratch) && startsWith(tolower(normalizePath(scratch,winslash='/')),
                                        paste0(tolower(out),'/'))) {
      unlink(scratch,recursive=TRUE,force=TRUE)
    }
  },add=TRUE)
  pwsh <- unname(Sys.which('pwsh'))
  if (!nzchar(pwsh)) stop('PowerShell 7 (pwsh) must be available on PATH.')
  psquote <- function(x) paste0("'",gsub("'","''",x,fixed=TRUE),"'")
  command <- paste('$ErrorActionPreference =',psquote('Stop'),
                   '; Expand-Archive -LiteralPath',psquote(archive),
                   '-DestinationPath',psquote(scratch))
  message('Extracting the original archive; the temporary DTA will be removed afterward.')
  status <- system2(pwsh,c('-NoProfile','-NonInteractive','-Command',shQuote(command)))
  dta <- file.path(scratch,'pu2022.dta')
  if (status!=0L || !file.exists(dta) || file.info(dta)$size!=entries$Length)
    stop('Archive extraction failed or extracted size differs.')
  columns <- c('SSUID','PNUM','MONTHCODE','EOWN_THR401','EHLTSTAT','EMJOB_401',
    'TAGE','ESEX','ERACE','EORIGIN','EEDUC','THINCPOV','EMS','RHLTHMTH','EDISABL')
  message('Reading identifiers and required analysis variables.')
  raw <- as.data.frame(haven::read_dta(dta,col_select=tidyselect::all_of(columns)))
  stopifnot(setequal(names(raw),columns))
  raw <- raw[,columns]
  source_row <- seq_len(nrow(raw))
  december <- which(raw$MONTHCODE %in% 12)
  dec <- raw[december,]
  # Valid EMJOB_401 answers select retirement-account owners with a job in
  # December, with or without an account through their main employer or business.
  eligible <- which(dec$EMJOB_401 %in% c(1,2))
  d <- dec[eligible,]
  recode <- function(x,yes,no) ifelse(x %in% yes,1,ifelse(x %in% no,0,NA_real_))
  analysis <- data.frame(A=recode(d$EMJOB_401,1,2),Y=6-as.numeric(d$EHLTSTAT),
    age=as.numeric(d$TAGE),female=recode(d$ESEX,2,1),
    race_black=as.numeric(d$ERACE==2),race_asian=as.numeric(d$ERACE==3),
    hispanic=recode(d$EORIGIN,1,2),
    educ_cat=ifelse(d$EEDUC %in% 31:39,1,ifelse(d$EEDUC %in% 40:42,2,
      ifelse(d$EEDUC==43,3,ifelse(d$EEDUC %in% 44:46,4,NA_real_)))),
    povratio=as.numeric(d$THINCPOV),married=recode(d$EMS,c(1,2),3:6),
    insured_any=recode(d$RHLTHMTH,1,2),work_limited=recode(d$EDISABL,1,2))
  valid <- complete.cases(analysis) & analysis$Y %in% 1:5
  id_valid <- !is.na(d$SSUID) & !is.na(d$PNUM) &
    !as.character(d$SSUID) %in% c('','-999') & !as.character(d$PNUM) %in% c('','-999')
  keep <- which(valid & id_valid)
  analysis <- analysis[keep,]
  stopifnot(nrow(analysis)>1L,all(is.finite(as.matrix(analysis))),
            all(d$EOWN_THR401[keep] %in% 1))
  ids <- data.frame(SSUID=as.character(d$SSUID[keep]),PNUM=as.character(d$PNUM[keep]),
    MONTHCODE=as.integer(d$MONTHCODE[keep]),source_row=source_row[december[eligible[keep]]],
    stringsAsFactors=FALSE)
  ids$person_id <- paste(ids$SSUID,ids$PNUM,sep=':')
  if (anyDuplicated(ids$person_id)) stop('Duplicate people remain among eligible December records; no dataset saved.')
  # Standardize on the final December analytic sample, before any bootstrap.
  scaling <- data.frame(variable=c('age','povratio'),
    center=vapply(analysis[c('age','povratio')],mean,numeric(1)),
    scale=vapply(analysis[c('age','povratio')],sd,numeric(1)))
  stopifnot(all(is.finite(scaling$scale)),all(scaling$scale>0))
  original_units <- data.frame(age_years=analysis$age,income_to_poverty_ratio=analysis$povratio)
  analysis[c('age','povratio')] <- lapply(analysis[c('age','povratio')],
                                         function(x) as.numeric(scale(x)))
  rownames(analysis) <- NULL
  dataset <- cbind(ids,analysis,original_units)
  rownames(dataset) <- NULL
  metadata <- list(source_zip=basename(archive),source_zip_md5=unname(before[1]),
    selection='All December records with valid EMJOB_401, complete analysis variables, and valid person identifiers',
    person_key=c('SSUID','PNUM'),reference_year=2021L,
    outcome_coding='Y=6-EHLTSTAT: 1=Poor,...,5=Excellent (existing application coding)',
    treatment_coding='A=1 if EMJOB_401=1; A=0 if EMJOB_401=2',
    analysis_columns=names(analysis),scaling=scaling,random_sampling=FALSE,
    standardization_population='Final eligible complete-case December sample',
    created=as.character(Sys.Date()))
  metadata$sample_flow <- data.frame(stage=c('Original person-month records',
    'December records','December records with valid EMJOB_401',
    'Complete analysis variables and valid outcome',
    'Final records with valid person identifiers','Distinct people','Treated','Control'),
    count=c(nrow(raw),nrow(dec),nrow(d),sum(valid),nrow(dataset),
      length(unique(dataset$person_id)),sum(dataset$A==1),sum(dataset$A==0)))
  metadata$outcome_counts <- table(A=dataset$A,Y=dataset$Y)
  metadata$session_info <- capture.output(sessionInfo())
  metadata$data_dictionary <- paste0('https://www2.census.gov/programs-surveys/sipp/',
    'tech-documentation/data-dictionaries/2022/2022_SIPP_Data_Dictionary.pdf')
  after <- tools::md5sum(protected)
  stopifnot(identical(before,after))
  metadata$source_integrity <- data.frame(path=basename(names(before)),
    md5_before=unname(before),md5_after=unname(after),unchanged=unname(before==after))
  attr(dataset,'preprocessing') <- metadata
  # The only persistent output is this RDS. Identifiers remain aligned with
  # their records when the data frame is subsetted or resampled.
  # To obtain the 12 analysis columns after loading dat, use:
  # dat[, attr(dat, 'preprocessing')$analysis_columns, drop=FALSE]
  # Explicitly select covariates in models; do not use all remaining columns.
  # This script does not run policy models or write analysis results.
  saveRDS(dataset,output_file)
  saved <- readRDS(output_file)
  stopifnot(identical(saved,dataset),all(saved$MONTHCODE==12L),
    !anyDuplicated(saved$person_id),nrow(saved)==8319L,
    sum(saved$A==1)==6624L,sum(saved$A==0)==1695L)
  print(metadata$sample_flow,row.names=FALSE)
  message('Saved and verified sipp_data_dec.rds; original files unchanged.')
  invisible(dataset)
}

if (!isTRUE(getOption('opl.sipp.generate_functions_only',FALSE)))
  generate_sipp_december()
