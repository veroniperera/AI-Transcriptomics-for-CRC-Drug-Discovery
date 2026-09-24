options(stringsAsFactors = FALSE)
SEED <- 42
set.seed(SEED)

DATA <- file.path(getwd(), "data")
RES  <- file.path(getwd(), "results")
for (d in c(file.path(DATA, "raw"), file.path(DATA, "processed"),
            file.path(RES, "figures"), file.path(RES, "tables"), file.path(RES, "rds"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOGFC_CUT      <- 2
FDR_CUT        <- 0.05
QUANTILE_NORM  <- TRUE
EBAYES_TREND   <- FALSE
WGCNA_TOP_MAD  <- NULL
SOFT_POWER     <- 16
MM_CUT         <- 0.8
GS_CUT         <- 0.2
ME_TRAIT_R     <- 0.3

HUB_TOP10 <- c("ZZZ3", "SSBP4", "IL12RB2", "UCHL3", "SIX5",
               "RNF19A", "SYVN1", "SHMT2", "CD2", "DDC")
