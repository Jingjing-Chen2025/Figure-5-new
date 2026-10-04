# FIGURE 5, PANEL b - step 2 | The compound-nomination pipeline: the list screened in panel b.
#
# WHAT THE ORIGINAL PIPELINE DID, AND WHY IT WAS REPLACED.
#
#   Branch A step 3 intersected a promoter-derived transcription-factor list with
#   `cluster.markers.sig.csv %>% filter(cluster == "c5")`. That file is the WHOLE-TISSUE clustering of the
#   single mouse dataset GSE235036, and in that numbering c5 is the ENTIRE stromal compartment - not a state.
#   Its markers therefore separate stroma from epithelium, in one mouse, and cannot distinguish PI16+ from
#   POSTN+ because the two were never resolved. It is also unrelated to atlas c5, the POSTN-high mesenchymal
#   state of the mouse mesenchymal atlas, with which it is easily confused. For the record: 1,800 c5 markers,
#   96 of them transcription factors, only 15 of which survive the state analysis below.
#
#   The branch also began with a HOCOMOCO motif scan of the MDK promoter and an active-chromatin filter. That
#   step discarded transcription factors on motif presence at one of the two ligand genes - TNFSF12 was never
#   scanned the same way - and does not discriminate regardless: 1,116 TFs have a hit in that window.
#
# WHAT REPLACES IT. The PI16+ and POSTN+ per-cell states of Figures 1 and 2, in the 15 single-cell datasets
# those figures use (11 human, 4 mouse). Within the fibroblast lineage of each dataset:
#     PI16+   per-cell call on the 200-gene progenitor signature (PI16, DPT, MFAP5, CLEC3B, CD34, PDGFRA ...)
#     POSTN+  non-PI16 lineage cell in the top quintile of the 100-gene POSTN-high signature
#     rest    the remainder of the lineage
# A transcription factor enters the list when it is differentially expressed in one of those two states, in the
# same direction, in at least three datasets of a species (tier 1), or when one of the supporting analyses run
# on the same cells returns it (tier 2).
#
# THE REST OF THE ARCHITECTURE
#   upstream  kinases and receptors of those TFs, from OmniPath
#   drugs     Open Targets, kept where maximum clinical stage = approval
#   branch B  LINCS L1000 compounds opposing the MDK/TNFSF12 programme in the Beltran NEPC-like biopsies, with
#             the AR-altered biopsies as a contrast (R/fig5b_01_signature_reversal.R)
#   ranking   cross-evidence over those three axes
#
# DECLARED SUBSTITUTIONS: upstream kinases come from OmniPath because the full PhosphoSitePlus kinase-substrate
# dataset is not redistributable; approved drugs come from Open Targets because Selleck Chemicals publishes no
# programmatic catalogue. Every compound can be looked up in the Selleck catalogue by name.
#
# Inputs  the per-state differential-expression tables under $F5_DATA (see README, "Upstream dependencies"),
#         OmniPath interaction and receptor tables, data/l1000_opposite_beltran.csv
# Output  $F5_DATA/compound_pipeline/compound_list.csv  (also shipped as data/compound_list.csv)
# Usage   Rscript R/fig5b_02_compound_pipeline.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(writexl); library(jsonlite); library(httr); library(tibble)
  library(org.Hs.eg.db); library(AnnotationDbi); library(data.table)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
OUT <- f5_data("compound_pipeline"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42); options(width = 200)
TESTED <- c("ATORVASTATIN","BETAMETHASONE","FULVESTRANT","MIRABEGRON","PAROXETINE",
            "SALBUTAMOL","SILDENAFIL","SORAFENIB","TAMOXIFEN","TRAMETINIB")
rd <- function(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL

# ================================================================ BRANCH A: TFs of the Figure 1/2 states
dep <- rd(f5_data("state_tf_network", "depth_by_state.csv"))
msg("Figure 1/2 cells: %d datasets (%d human, %d mouse); PI16+ %s cells, POSTN+ %s, rest of lineage %s",
    n_distinct(dep$dataset), n_distinct(dep$dataset[dep$species == "human"]),
    n_distinct(dep$dataset[dep$species == "mouse"]),
    format(sum(dep$n[dep$state == "ASPC"]), big.mark = ","),
    format(sum(dep$n[dep$state == "C5_high"]), big.mark = ","),
    format(sum(dep$n[dep$state == "Fib_other"]), big.mark = ","))

tiers <- list()
# tier 1 - the declared primary criterion: recurrent in >= 3 datasets of a species, same state, same direction
x <- rd(f5_data("state_tf_network", "tf_recurrent.csv"))
if (!is.null(x)) tiers$rec <- x %>% filter(state != "Rest_of_lineage") %>%
  transmute(gene = toupper(gene), state, species, tier = 1L,
            evidence = sprintf("recurrent in %d/%d %s datasets", n_hit, n_tested, species))
# tier 1 - the same cells, three stricter readings of the same contrast
x <- rd(f5_data("nepc_tf_replication", "replicated_TFs.csv"))
if (!is.null(x)) tiers$rep <- x %>% filter(state != "Rest_of_lineage") %>%
  transmute(gene = toupper(gene), state, species = "both", tier = 1L, evidence = "replicated in human and mouse")
x <- rd(f5_data("collectri_per_dataset", "shared_hubs.csv"))
if (!is.null(x)) tiers$hub <- x %>% filter(state != "Rest_of_lineage", human_datasets >= 1, mouse_datasets >= 1) %>%
  transmute(gene = toupper(gene), state, species = "both", tier = 1L, evidence = "CollecTRI network hub in both species")
x <- rd(f5_data("nepc_tf_consistency", "tf_sample_counts.csv"))
if (!is.null(x)) tiers$con <- x %>% filter(state != "Rest_of_lineage") %>%
  group_by(gene, state) %>% summarise(h = any(species == "human" & n_up / n >= 6/7),
                                      m = any(species == "mouse" & n_up / n >= .85), .groups = "drop") %>%
  filter(h, m) %>% transmute(gene = toupper(gene), state, species = "both", tier = 1L,
                             evidence = "consistent across samples in both species")
# tier 2 - the same cells, single-species or single-pass differential expression
x <- rd(f5_data("human_tf_strategy", "strategy_comparison.csv"))
if (!is.null(x)) tiers$hde <- x %>% filter(strategy == "ALL human, all cells", fdr < 0.05, lfc > 0,
                                           state != "Rest_of_lineage") %>%
  transmute(gene = toupper(gene), state, species = "human", tier = 2L, evidence = "all-human differential expression")
for (sp in c("human", "mouse")) {
  x <- rd(f5_data("nepc_tf_expression", sprintf("%s_conventional.csv", sp)))
  if (!is.null(x)) tiers[[sp]] <- x %>% filter(fdr_wilcox < 0.05, lfc > 0, state != "Rest_of_lineage") %>%
    transmute(gene = toupper(gene), state, species = sp, tier = 2L, evidence = sprintf("%s, conventional FDR", sp))
}
STATE <- bind_rows(tiers) %>% group_by(gene, state) %>%
  summarise(tier = min(tier), species = paste(sort(unique(species)), collapse = "/"),
            evidence = paste(sort(unique(evidence)), collapse = "; "), .groups = "drop") %>%
  arrange(tier, state, gene)
write.csv(STATE, file.path(OUT, "state_TFs_fig1_2.csv"), row.names = FALSE)
msg("state TFs: %d (%d PI16+, %d POSTN+; %d tier 1, of which %d recurrent)",
    n_distinct(STATE$gene), sum(STATE$state == "PI16_progenitor"), sum(STATE$state == "POSTN_activated"),
    sum(STATE$tier == 1), sum(grepl("recurrent", STATE$evidence)))

# ================================================================ upstream kinases and receptors (OmniPath)
omni_file <- f5_data("omnipath", "interactions.tsv"); rec_file <- f5_data("omnipath", "intercell_receptors.tsv")
if (!file.exists(omni_file)) {
  dir.create(dirname(omni_file), showWarnings = FALSE, recursive = TRUE)
  msg("downloading OmniPath ...")
  download.file("https://omnipathdb.org/interactions?genesymbols=1&fields=is_directed", omni_file, quiet = TRUE)
  download.file("https://omnipathdb.org/intercell?categories=receptor&genesymbols=1", rec_file, quiet = TRUE)
}
omni <- fread(omni_file); rec <- fread(rec_file)
RECEPTORS <- unique(rec$genesymbol[rec$consensus_score >= 5])
# PhosphoSitePlus is licence-restricted and is not redistributed here; a kinase is recognised either from a
# local copy of its kinase-substrate table, if the user has one, or from the gene-family pattern below.
psp_file <- f5_data("PhosphoSitePlus_Kinase_Substrate_List.xlsx")
psp <- if (file.exists(psp_file) && requireNamespace("readxl", quietly = TRUE))
  toupper(trimws(suppressMessages(readxl::read_excel(psp_file))$KINASE)) else character(0)
KIN_RE <- "^(MAP[0-9K]|CDK|CSNK|GSK|PRKA|PRKC|PRKD|AKT|SRC|JAK|RPS6K|CAMK|PLK|AURK|ROCK|IKBK|CHUK|MTOR|PIK3|ABL|EGFR|ERBB|FGFR|PDGFR|KDR|FLT|KIT|MET|RET|ALK|BRAF|RAF1|ARAF|SYK|BTK|LCK|FYN|TGFBR|ACVR|NTRK|INSR|IGF1R)"
up <- omni[is_directed == TRUE & target_genesymbol %in% STATE$gene,
           .(regulator = source_genesymbol, tf = target_genesymbol)] %>%
  distinct() %>% mutate(class = ifelse(regulator %in% RECEPTORS, "receptor",
                                ifelse(regulator %in% psp | grepl(KIN_RE, regulator), "kinase", "other"))) %>%
  filter(class != "other")
msg("upstream: %d kinases/receptors for %d of the state TFs", n_distinct(up$regulator), n_distinct(up$tf))

TARGETS <- bind_rows(
  STATE %>% transmute(gene, role = "state TF", tier, detail = state),
  up %>% transmute(gene = regulator, role = paste0("upstream ", class), tier = 3L, detail = paste0("regulates ", tf))) %>%
  group_by(gene) %>% summarise(role = paste(sort(unique(role)), collapse = "; "), tier = min(tier),
                               detail = paste(sort(unique(detail)), collapse = "; "), .groups = "drop")
ens <- AnnotationDbi::select(org.Hs.eg.db, keys = TARGETS$gene, keytype = "SYMBOL", columns = "ENSEMBL")
TARGETS <- TARGETS %>%
  left_join(ens %>% filter(!is.na(ENSEMBL)) %>% group_by(SYMBOL) %>% dplyr::slice(1) %>% ungroup(),
            by = c("gene" = "SYMBOL")) %>% filter(!is.na(ENSEMBL)) %>% arrange(tier, gene)
write.csv(TARGETS, file.path(OUT, "target_genes.csv"), row.names = FALSE)
msg("targets to query: %d (%d state TFs, %d upstream)", nrow(TARGETS),
    sum(grepl("state TF", TARGETS$role)), sum(!grepl("state TF", TARGETS$role)))

# ================================================================ approved drugs: Open Targets
CACHE_OT <- file.path(OUT, "opentargets.rds")
Q <- '{ target(ensemblId:"%s"){ approvedSymbol drugAndClinicalCandidates{ rows{ drug{ name drugType maximumClinicalStage mechanismsOfAction{ rows{ actionType } } } } } } }'
ot <- function(e) {
  r <- tryCatch(httr::POST("https://api.platform.opentargets.org/api/v4/graphql", httr::content_type_json(),
                           body = toJSON(list(query = sprintf(Q, e)), auto_unbox = TRUE), httr::timeout(45)),
                error = function(x) NULL)
  if (is.null(r) || httr::status_code(r) != 200) return(NULL)
  j <- tryCatch(fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE), error = function(x) NULL)
  rows <- j$data$target$drugAndClinicalCandidates$rows
  if (is.null(rows) || !length(rows)) return(NULL)
  bind_rows(lapply(rows, function(z) data.frame(
    target_symbol = j$data$target$approvedSymbol %||% NA_character_, drug = z$drug$name %||% NA_character_,
    drug_type = z$drug$drugType %||% NA_character_, max_stage = z$drug$maximumClinicalStage %||% NA_character_,
    action = paste(unique(unlist(lapply(z$drug$mechanismsOfAction$rows, `[[`, "actionType"))), collapse = ";"),
    stringsAsFactors = FALSE)))
}
CA <- if (file.exists(CACHE_OT)) readRDS(CACHE_OT) else list()
todo <- setdiff(TARGETS$ENSEMBL, names(CA))
msg("Open Targets: %d cached, %d to fetch", length(CA), length(todo))
for (i in seq_along(todo)) {
  CA[[todo[i]]] <- ot(todo[i]) %||% data.frame()
  if (i %% 50 == 0) { msg("   %d/%d", i, length(todo)); saveRDS(CA, CACHE_OT) }
}
saveRDS(CA, CACHE_OT)
DR <- bind_rows(CA[TARGETS$ENSEMBL]) %>% filter(max_stage == "APPROVAL") %>%
  left_join(TARGETS %>% dplyr::select(gene, role, tier), by = c("target_symbol" = "gene"))
msg("approved drug-target rows: %d over %d targets", nrow(DR), n_distinct(DR$target_symbol))

# ================================================================ BRANCH B: L1000 reversal, Beltran NEPC
LB <- rd(f5_data("compound_beltran", "l1000_opposite_beltran.csv")) %||%
      rd(f5_repo("data", "l1000_opposite_beltran.csv"))
NEPC_REV <- LB %>% filter(cohort == "NEPC_like") %>% group_by(drug) %>%
  summarise(n_sig_NEPC = n(), cells_NEPC = paste(sort(unique(cell)), collapse = ";"), .groups = "drop")
AR_REV   <- LB %>% filter(cohort == "AR_altered") %>% group_by(drug) %>%
  summarise(n_sig_AR = n(), .groups = "drop")
msg("branch B: %d compounds oppose the programme in the NEPC-like biopsies, %d of them not in AR-altered",
    nrow(NEPC_REV), sum(!NEPC_REV$drug %in% AR_REV$drug))

# ================================================================ cross-evidence ranking
# Three independent axes: acts on a state TF; acts on an upstream kinase or receptor; reverses the ligand
# programme. A compound supported on more than one axis ranks above one supported on a single axis.
CE <- DR %>% group_by(drug, drug_type) %>%
  summarise(state_TF_targets = paste(sort(unique(target_symbol[grepl("state TF", role)])), collapse = ";"),
            upstream_targets = paste(sort(unique(target_symbol[!grepl("state TF", role)])), collapse = ";"),
            n_state_TF = n_distinct(target_symbol[grepl("state TF", role)]),
            n_upstream = n_distinct(target_symbol[!grepl("state TF", role)]),
            actions = paste(sort(unique(unlist(strsplit(action, ";")))), collapse = ";"),
            best_tier = min(tier), .groups = "drop") %>%
  mutate(DRUG = toupper(drug)) %>%
  full_join(NEPC_REV %>% dplyr::rename(DRUG = drug), by = "DRUG") %>%
  left_join(AR_REV   %>% dplyr::rename(DRUG = drug), by = "DRUG") %>%
  mutate(across(c(n_state_TF, n_upstream, n_sig_NEPC, n_sig_AR), ~ coalesce(.x, 0)),
         drug = coalesce(drug, DRUG),
         axis_state_TF = n_state_TF > 0, axis_upstream = n_upstream > 0, axis_L1000 = n_sig_NEPC > 0,
         NEPC_specific = n_sig_NEPC > 0 & n_sig_AR == 0,
         n_axes = axis_state_TF + axis_upstream + axis_L1000,
         score = 3 * n_state_TF + n_upstream + 2 * n_sig_NEPC + 2 * NEPC_specific,
         tested_in_Fig5b = vapply(DRUG, function(x) any(vapply(TESTED, function(t) grepl(t, x, fixed = TRUE),
                                                               logical(1))), logical(1))) %>%
  arrange(desc(n_axes), desc(score), drug)
write.csv(CE, file.path(OUT, "compound_list.csv"), row.names = FALSE)

acct <- bind_rows(lapply(TESTED, function(t) {
  i <- grep(t, CE$DRUG, fixed = TRUE)
  if (!length(i)) return(tibble(compound = t, in_pipeline = "NO", route = "-", targets = "-"))
  r <- CE[i[1], ]
  tibble(compound = t, in_pipeline = "yes",
         route = paste(c("drug-target (state TF)", "drug-target (upstream)", "L1000 reversal, Beltran NEPC")[
                         c(r$axis_state_TF, r$axis_upstream, r$axis_L1000)], collapse = " + "),
         targets = paste(c(r$state_TF_targets, r$upstream_targets)[
                           nzchar(c(r$state_TF_targets, r$upstream_targets))], collapse = " / "))
}))
write.csv(acct, file.path(OUT, "ten_tested_compounds.csv"), row.names = FALSE)
short <- CE %>% filter(axis_state_TF, grepl("ANTAGONIST|INHIBITOR|DEGRADER", actions)) %>% arrange(desc(score), drug)
write_xlsx(list(`Compound list` = as.data.frame(CE),
                `Ten tested compounds` = as.data.frame(acct),
                `Shortlist (antagonists)` = as.data.frame(short),
                `State TFs (Fig 1 and 2)` = as.data.frame(STATE),
                `Target genes` = as.data.frame(TARGETS)),
           file.path(OUT, "Compound list.xlsx"))

cat(sprintf("\n================ COMPOUND LIST: %d compounds\n", nrow(CE)))
cat(sprintf("   FDA-approved, from drug-target mapping : %d\n", sum(!is.na(CE$drug_type))))
cat(sprintf("   acting on a PI16+/POSTN+ state TF      : %d\n", sum(CE$axis_state_TF)))
cat(sprintf("   acting on an upstream kinase/receptor  : %d\n", sum(CE$axis_upstream)))
cat(sprintf("   reversing the Beltran NEPC signature   : %d  (NEPC-specific %d)\n",
            sum(CE$axis_L1000), sum(CE$NEPC_specific)))
cat(sprintf("   supported on 2+ axes                   : %d\n", sum(CE$n_axes >= 2)))
cat("\n-- top 25 by cross-evidence score\n")
print(as.data.frame(CE %>% transmute(drug, n_axes, score, state_TF = substr(state_TF_targets, 1, 20),
      upstream = substr(upstream_targets, 1, 24), L1000 = n_sig_NEPC, tested = tested_in_Fig5b) %>% head(25)),
      row.names = FALSE)
cat("\n-- the ten compounds tested in panel b\n")
print(as.data.frame(acct), row.names = FALSE)
cat("\n-- the state transcription factors that drive branch A\n")
for (s in c("PI16_progenitor", "POSTN_activated")) {
  g <- STATE %>% filter(state == s)
  cat(sprintf("  %-16s tier 1 (%2d): %s\n", s, sum(g$tier == 1), paste(sort(g$gene[g$tier == 1]), collapse = ", ")))
  cat(sprintf("  %-16s tier 2 (%2d)\n", "", sum(g$tier == 2)))
}
msg("written %s", file.path(OUT, "compound_list.csv"))
