# FIGURE 5, PANEL b - step 1 | Signature reversal in the Beltran NEPC cohort (branch B of the pipeline).
#
# An earlier version of this arm ran on TCGA PRAD, SU2C and an ARCHS4 adipose compendium. TCGA PRAD is primary,
# hormone-naive, overwhelmingly adenocarcinoma prostate tissue: it is the wrong cohort for a question about
# neuroendocrine disease. This runs instead on the Beltran metastatic cohort (RR1384), restricted to its
# NEPC-like biopsies, with the AR-altered biopsies of the SAME cohort as a contrast - so a compound that
# reverses the programme in NEPC-like disease specifically can be told apart from one that reverses it in any
# metastatic prostate tissue.
#
#   signature  top and bottom 250 genes correlated (Spearman) with MDK and with TNFSF12 within each disease
#              group, submitted to LINCS L1000FWD; the signatures that OPPOSE that programme are collected.
#
# Inputs  RR1384 counts ($F5_RR1384), the disease-group assignment from R/fig5a_01_bmi_composition.R
# Output  $F5_DATA/compound_beltran/l1000_opposite_beltran.csv  (also shipped as data/l1000_opposite_beltran.csv)
# Usage   Rscript R/fig5b_01_signature_reversal.R
suppressPackageStartupMessages({library(dplyr); library(jsonlite); library(httr); library(tibble)})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
OUT <- f5_data("compound_beltran"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42); options(width = 200)
L1000 <- "https://maayanlab.cloud/L1000FWD/"

# ---------------------------------------------------------------- the cohort
sc  <- read.csv(f5_data("wcm_bmi_composition", "wcm_sample_scores.csv"))
raw <- read.csv(file.path(F5_RR1384, "RR1384_all_PM_cases_counts.csv"), check.names = FALSE)
g <- toupper(raw[[1]]); cnt <- rowsum(as.matrix(raw[, -1]), group = g); storage.mode(cnt) <- "double"
cnt <- cnt[rowSums(cnt > 0) >= 0.2 * ncol(cnt), , drop = FALSE]
cpm <- log2(sweep(cnt, 2, pmax(colSums(cnt), 1), "/") * 1e6 + 1)
grp <- list(NEPC_like  = intersect(sc$sample[sc$category == "NEPC_like"],  colnames(cpm)),
            AR_altered = intersect(sc$sample[sc$category == "AR_altered"], colnames(cpm)))
msg("Beltran RR1384: %d NEPC-like and %d AR-altered biopsies", length(grp$NEPC_like), length(grp$AR_altered))

sig_for <- function(cols, n = 250) {
  M <- cpm[, cols, drop = FALSE]; M <- M[apply(M, 1, sd) > 0, , drop = FALSE]
  out <- list()
  for (gene in c("MDK", "TNFSF12")) {
    if (!gene %in% rownames(M)) { warn_msg("%s absent", gene); next }
    r <- as.numeric(cor(t(M), M[gene, ], method = "spearman")); names(r) <- rownames(M)
    r <- r[names(r) != gene & is.finite(r)]
    out[[gene]] <- list(up = names(sort(r, decreasing = TRUE))[seq_len(n)], down = names(sort(r))[seq_len(n)])
  }
  list(up = unique(unlist(lapply(out, `[[`, "up"))), down = unique(unlist(lapply(out, `[[`, "down"))))
}
SIG <- lapply(grp, sig_for)
for (k in names(SIG)) msg("%s signature: %d up, %d down", k, length(SIG[[k]]$up), length(SIG[[k]]$down))

# ---------------------------------------------------------------- L1000FWD
l1000 <- function(up, down) {
  r <- tryCatch(httr::RETRY("POST", paste0(L1000, "sig_search"), body = list(up_genes = up, down_genes = down),
                            encode = "json", times = 3, pause_min = 1, httr::timeout(90)), error = function(e) NULL)
  if (is.null(r) || status_code(r) != 200) { warn_msg("sig_search failed"); return(NULL) }
  rid <- fromJSON(content(r, "text", encoding = "UTF-8"))$result_id
  if (is.null(rid)) return(NULL)
  Sys.sleep(2)
  t <- tryCatch(httr::RETRY("GET", paste0(L1000, "result/topn/", rid), times = 4, pause_min = 2, httr::timeout(90)),
                error = function(e) NULL)
  if (is.null(t) || status_code(t) != 200) return(NULL)
  js <- fromJSON(gsub("NaN|(-?)Infinity", "null", content(t, "text", encoding = "UTF-8")))
  opp <- js$opposite
  if (is.null(opp) || !length(opp)) return(NULL)
  sigs <- if (is.data.frame(opp)) opp else bind_rows(opp)
  ids <- sigs$sig_id %||% sigs[[1]]
  bind_rows(lapply(ids, function(s) {
    rr <- tryCatch(httr::GET(paste0(L1000, "sig/", URLencode(s, reserved = TRUE)), httr::timeout(45)),
                   error = function(e) NULL)
    if (is.null(rr) || status_code(rr) != 200) return(NULL)
    j <- tryCatch(fromJSON(gsub("NaN|(-?)Infinity", "null", content(rr, "text", encoding = "UTF-8"))),
                  error = function(e) NULL)
    if (is.null(j)) return(NULL)
    tibble(sig_id = s, drug = toupper(j$pert_desc %||% NA_character_),
           pert_id = j$pert_id %||% NA_character_, cell = j$cell_id %||% NA_character_)
  }))
}
CACHE <- file.path(OUT, "l1000_raw.rds")
REV <- if (file.exists(CACHE)) readRDS(CACHE) else list()
for (k in setdiff(names(SIG), names(REV))) {
  msg("L1000FWD: querying %s ...", k)
  REV[[k]] <- l1000(SIG[[k]]$up, SIG[[k]]$down)
  msg("   %s: %d opposing signatures", k, if (is.null(REV[[k]])) 0L else nrow(REV[[k]]))
}
saveRDS(REV, CACHE)

# "-666" and the BRD-xxxxx identifiers are L1000's placeholders for a compound with no released name; they
# cannot be sourced or interpreted, so they are dropped rather than counted.
L <- bind_rows(lapply(names(REV), function(k) if (is.null(REV[[k]])) NULL else REV[[k]] %>% mutate(cohort = k))) %>%
  filter(!is.na(drug), !grepl("CONFIDENTIAL", drug), drug != "-666", !grepl("^BRD-", drug))
write.csv(L, file.path(OUT, "l1000_opposite_beltran.csv"), row.names = FALSE)

nep <- unique(L$drug[L$cohort == "NEPC_like"]); ar <- unique(L$drug[L$cohort == "AR_altered"])
msg("%d named compounds oppose the programme in the NEPC-like biopsies, %d of them not in the AR-altered contrast",
    length(nep), length(setdiff(nep, ar)))
cat("\n", paste(sort(nep), collapse = ", "), "\n", sep = "")
msg("written %s", file.path(OUT, "l1000_opposite_beltran.csv"))
