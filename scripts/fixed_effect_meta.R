#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs")
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

required_args <- c("manifest", "analysis", "out", "qc-out")
args <- parse_args(commandArgs(trailingOnly = TRUE))
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0L) {
  stop("Missing arguments: ", paste(missing_args, collapse = ", "))
}

if (!args$analysis %in% c("european", "multi")) {
  stop("--analysis must be european or multi")
}

required_columns <- c("CHR", "POS", "EA", "OA", "BETA", "SE", "P", "EAF", "N")
valid_alleles <- c("A", "C", "G", "T")

is_palindromic <- function(ea, oa) {
  paste0(ea, oa) %in% c("AT", "TA", "CG", "GC")
}

complement <- function(allele) {
  map <- c(A = "T", T = "A", C = "G", G = "C")
  unname(map[allele])
}

read_cohort <- function(path, cohort, ancestry) {
  dat <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  absent <- setdiff(required_columns, names(dat))
  if (length(absent) > 0L) {
    stop(cohort, " is missing columns: ", paste(absent, collapse = ", "))
  }

  dat <- dat[required_columns]
  dat$EA <- toupper(dat$EA)
  dat$OA <- toupper(dat$OA)
  dat$COHORT <- cohort
  dat$ANCESTRY <- ancestry
  dat$KEY <- paste(dat$CHR, dat$POS, sep = ":")

  if (anyDuplicated(dat$KEY)) {
    stop(cohort, " contains duplicate chromosome-position keys")
  }
  dat
}

manifest <- read.delim(args$manifest, stringsAsFactors = FALSE, check.names = FALSE)
manifest_required <- c("cohort", "ancestry", "include_european", "path")
manifest_missing <- setdiff(manifest_required, names(manifest))
if (length(manifest_missing) > 0L) {
  stop("Manifest is missing columns: ", paste(manifest_missing, collapse = ", "))
}

if (args$analysis == "european") {
  keep <- toupper(as.character(manifest$include_european)) %in% c("TRUE", "T", "1")
  manifest <- manifest[keep, , drop = FALSE]
}
if (nrow(manifest) < 2L) {
  stop("At least two cohort inputs are required")
}

cohort_data <- lapply(seq_len(nrow(manifest)), function(i) {
  read_cohort(manifest$path[[i]], manifest$cohort[[i]], manifest$ancestry[[i]])
})
names(cohort_data) <- manifest$cohort

all_rows <- do.call(rbind, cohort_data)
reference_candidates <- all_rows[
  all_rows$EA %in% valid_alleles &
    all_rows$OA %in% valid_alleles &
    all_rows$EA != all_rows$OA,
]
reference <- reference_candidates[
  !duplicated(reference_candidates$KEY),
  c("KEY", "CHR", "POS", "EA", "OA")
]
names(reference)[names(reference) == "EA"] <- "REF_EA"
names(reference)[names(reference) == "OA"] <- "REF_OA"

harmonize_one <- function(dat) {
  x <- merge(dat, reference, by = c("KEY", "CHR", "POS"), all.x = TRUE, sort = FALSE)
  x$ACTION <- "mismatch_drop"
  x$BETA_H <- NA_real_
  x$EAF_H <- NA_real_

  numeric_ok <- is.finite(x$BETA) & is.finite(x$SE) & x$SE > 0 &
    is.finite(x$P) & x$P >= 0 & x$P <= 1 &
    is.finite(x$EAF) & x$EAF >= 0 & x$EAF <= 1 &
    is.finite(x$N) & x$N > 0
  allele_ok <- x$EA %in% valid_alleles & x$OA %in% valid_alleles & x$EA != x$OA
  pal <- is_palindromic(x$EA, x$OA)

  direct <- numeric_ok & allele_ok & !pal & x$EA == x$REF_EA & x$OA == x$REF_OA
  swapped <- numeric_ok & allele_ok & !pal & x$EA == x$REF_OA & x$OA == x$REF_EA
  comp_direct <- numeric_ok & allele_ok & !pal &
    complement(x$EA) == x$REF_EA & complement(x$OA) == x$REF_OA
  comp_swapped <- numeric_ok & allele_ok & !pal &
    complement(x$EA) == x$REF_OA & complement(x$OA) == x$REF_EA

  x$ACTION[direct] <- "direct"
  x$BETA_H[direct] <- x$BETA[direct]
  x$EAF_H[direct] <- x$EAF[direct]

  x$ACTION[swapped] <- "swapped"
  x$BETA_H[swapped] <- -x$BETA[swapped]
  x$EAF_H[swapped] <- 1 - x$EAF[swapped]

  x$ACTION[comp_direct] <- "strand_complement"
  x$BETA_H[comp_direct] <- x$BETA[comp_direct]
  x$EAF_H[comp_direct] <- x$EAF[comp_direct]

  x$ACTION[comp_swapped] <- "strand_complement_swapped"
  x$BETA_H[comp_swapped] <- -x$BETA[comp_swapped]
  x$EAF_H[comp_swapped] <- 1 - x$EAF[comp_swapped]

  x$ACTION[pal & numeric_ok & allele_ok] <- "palindromic_drop"
  x$ACTION[!numeric_ok] <- "invalid_numeric_drop"
  x$ACTION[!allele_ok] <- "invalid_allele_drop"
  x
}

harmonized <- do.call(rbind, lapply(cohort_data, harmonize_one))

qc <- aggregate(
  list(N_ROWS = rep(1L, nrow(harmonized))),
  by = list(COHORT = harmonized$COHORT, ACTION = harmonized$ACTION),
  FUN = sum
)
qc <- qc[order(qc$COHORT, qc$ACTION), ]

usable <- harmonized[!is.na(harmonized$BETA_H), ]
variant_groups <- split(usable, usable$KEY)

meta_one <- function(x) {
  if (nrow(x) < 2L) {
    return(NULL)
  }
  w <- 1 / (x$SE ^ 2)
  beta <- sum(w * x$BETA_H) / sum(w)
  se <- sqrt(1 / sum(w))
  z <- beta / se
  p <- 2 * pnorm(-abs(z))
  q <- sum(w * (x$BETA_H - beta) ^ 2)
  df <- nrow(x) - 1L
  q_p <- pchisq(q, df = df, lower.tail = FALSE)
  i2 <- if (q > 0) max(0, (q - df) / q) else 0
  eaf <- weighted.mean(x$EAF_H, w = x$N)

  data.frame(
    CHR = x$CHR[[1]],
    POS = x$POS[[1]],
    EA = x$REF_EA[[1]],
    OA = x$REF_OA[[1]],
    K = nrow(x),
    N = sum(x$N),
    EAF = eaf,
    BETA = beta,
    SE = se,
    P = p,
    OR = exp(beta),
    CI_LOWER = exp(beta - 1.96 * se),
    CI_UPPER = exp(beta + 1.96 * se),
    Q = q,
    Q_P = q_p,
    I2 = i2,
    stringsAsFactors = FALSE
  )
}

meta_list <- lapply(variant_groups, meta_one)
meta_list <- meta_list[!vapply(meta_list, is.null, logical(1))]
if (length(meta_list) == 0L) {
  stop("No variants remained in at least two cohorts after harmonization")
}
meta <- do.call(rbind, meta_list)
meta <- meta[order(as.numeric(meta$CHR), meta$POS), ]

dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args[["qc-out"]]), recursive = TRUE, showWarnings = FALSE)
write.table(meta, args$out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qc, args[["qc-out"]], sep = "\t", quote = FALSE, row.names = FALSE)

message("Analysis: ", args$analysis)
message("Cohorts: ", paste(manifest$cohort, collapse = ", "))
message("Meta-analyzed variants: ", nrow(meta))
message("Wrote: ", args$out)
message("Wrote: ", args[["qc-out"]])
