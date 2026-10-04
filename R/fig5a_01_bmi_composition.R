# FIGURE 5, PANEL a - step 1 | Does host adiposity change the COMPOSITION of the stroma, or only its amount?
#
# Beltran (RR1384) bulk RNA-seq joined to the clinical record (BMI, hormone therapy, histology).
#
# The problem this is built to avoid: every stromal signature scores high wherever there is more stroma, so the
# PI16+ and POSTN+ scores are collinear in bulk (Spearman ~0.9). A raw "the PI16+ score rises with BMI" result is
# therefore indistinguishable from "stromal content rises with BMI". The primary outcome is a CONTRAST,
#     delta = z(PI16+ specific) - z(POSTN+ specific),
# which is invariant to a shared multiplicative content term: a biopsy with twice the stroma raises both scores
# and leaves delta unchanged. The question becomes whether adiposity changes the composition of the stroma, not
# its amount. The total stromal score is reported alongside and entered as a covariate, so a content effect
# cannot masquerade as a composition effect.
#
# Three mutually exclusive gene sets (verified disjoint before scoring):
#   PI16+ specific  = the 200-gene progenitor signature, minus POSTN+ genes, minus the pan-fibroblast signature
#   POSTN+ specific = the 100-gene POSTN-high signature, minus PI16+ genes, minus the pan-fibroblast signature
#   rest of stroma  = the pan-fibroblast signature, minus both
#
# PRE-SPECIFIED (written before the run, 2026-10-01):
#   primary    delta ~ BMI x category; the interaction P is the test. Supported only if the interaction holds
#              after adjustment for the total stromal score, and the two categories differ in slope direction.
#   secondary  each score separately, the total stromal score, and the androgen-response programme, against BMI.
#   Repeat biopsies per patient violate independence, so a mixed model with (1|record_id) is also reported and
#   the number of patients contributing more than one sample is stated.
# EXPLORATORY. With roughly 25 samples per category a null result is uninformative and is reported as such.
#
# Inputs  RR1384 counts + sample key ($F5_RR1384), the clinical table ($F5_CLINICAL), the signatures in data/
# Output  $F5_DATA/wcm_bmi_composition/{wcm_sample_scores.csv, wcm_bmi_interaction_stats.csv}
# Usage   Rscript R/fig5a_01_bmi_composition.R
suppressPackageStartupMessages({library(readxl); library(dplyr); library(tidyr)})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
OUT <- f5_data("wcm_bmi_composition"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42); options(width = 200)

# ---------------- gene sets
aspc <- read_set(f5_repo("data", "signatures", "PI16_PROGENITOR_200.txt"))
c5   <- read_set(f5_repo("data", "signatures", "POSTN_HIGH_100.txt"))
fib  <- read_set(f5_repo("data", "signatures", "FIBROBLAST_STROMA.txt"))
SIG <- list(PI16_specific  = setdiff(aspc, union(c5, fib)),
            POSTN_specific = setdiff(c5,   union(aspc, fib)),
            Rest_of_stroma = setdiff(fib,  union(aspc, c5)),
            Total_stroma   = unique(c(aspc, c5, fib)),
            AR             = HALLMARK_ANDROGEN_RESPONSE,
            NEPC           = NEPC_BELTRAN_CUSTOM_UP)
stopifnot(length(intersect(SIG$PI16_specific, SIG$POSTN_specific)) == 0,
          length(intersect(SIG$PI16_specific, SIG$Rest_of_stroma)) == 0,
          length(intersect(SIG$POSTN_specific, SIG$Rest_of_stroma)) == 0)
msg("gene sets: PI16+ %d | POSTN+ %d | rest %d | total stroma %d (disjoint verified)",
    length(SIG$PI16_specific), length(SIG$POSTN_specific), length(SIG$Rest_of_stroma), length(SIG$Total_stroma))

# ---------------- expression: counts -> log2 CPM, gene-wise z across samples
raw <- read.csv(file.path(F5_RR1384, "RR1384_all_PM_cases_counts.csv"), check.names = FALSE)
g <- toupper(raw[[1]]); raw <- as.matrix(raw[, -1])
cnt <- rowsum(raw, group = g); storage.mode(cnt) <- "double"   # duplicate symbols are summed
msg("counts: %d rows -> %d unique symbols", length(g), nrow(cnt))
cnt <- cnt[rowSums(cnt > 0) >= 0.2 * ncol(cnt), , drop = FALSE]
cpm <- log2(sweep(cnt, 2, pmax(colSums(cnt), 1), "/") * 1e6 + 1)
z   <- t(scale(t(cpm)))                                        # per gene across samples
msg("expression: %d genes x %d samples after filtering", nrow(z), ncol(z))
score <- sapply(SIG, function(gg) { gg <- intersect(gg, rownames(z)); colMeans(z[gg, , drop = FALSE]) })
sc <- as.data.frame(score); sc$sample <- rownames(sc)
for (n in names(SIG)) msg("  %-15s %3d of %3d genes present", n, length(intersect(SIG[[n]], rownames(z))), length(SIG[[n]]))

# ---------------- sample -> patient -> clinical
key <- read_excel(file.path(F5_RR1384, "Beltran2020_sample_key.xlsx")) %>%
  transmute(sample = `Sample name on matrix`, patient = Patient.ID, diagnosis = Diagnosis, site = Primary.site)
cl <- read_excel(F5_CLINICAL, guess_max = 5000)
num <- function(x) suppressWarnings(as.numeric(x))
hcat <- grep("drugs_hormonecate?o?gry___", names(cl), value = TRUE)
clin <- cl %>%
  transmute(record_id = as.character(record_id), histology = as.character(histology),
            bmi = num(bmi_measure),
            n_hormone = rowSums(sapply(cl[hcat], function(v) num(v) %in% 1), na.rm = TRUE)) %>%
  mutate(bmi = ifelse(is.finite(bmi) & bmi > 10 & bmi < 70, bmi, NA_real_))
# the clinical export carries many longitudinal rows per patient, so collapse to one row per record_id before
# joining; otherwise every biopsy is duplicated and the P values are meaningless
clin <- clin %>% group_by(record_id) %>%
  summarise(histology = first(na.omit(histology)), bmi = suppressWarnings(median(bmi, na.rm = TRUE)),
            n_hormone = max(n_hormone, na.rm = TRUE), .groups = "drop") %>%
  mutate(bmi = ifelse(is.finite(bmi), bmi, NA_real_))
msg("clinical records after collapsing: %d (BMI on %d)", nrow(clin), sum(is.finite(clin$bmi)))
# Sample_ID_1/2/3 hold the sample name in three spellings (PM0_LM1, PM0_LM1_RNA, PM0_LM1_Case_RNASeq), so the
# join is on sample, not on patient; the suffixes are stripped to make the three forms agree.
strip <- function(x) sub("_(RNA|CASE_RNASEQ)$", "", toupper(trimws(x)))
map <- cl %>% transmute(record_id = as.character(record_id), s1 = as.character(Sample_ID_1),
                        s2 = as.character(Sample_ID_2), s3 = as.character(Sample_ID_3)) %>%
  pivot_longer(c(s1, s2, s3), values_to = "sid") %>% filter(!is.na(sid), sid != "NA") %>%
  mutate(sid = strip(sid)) %>% distinct(record_id, sid)
dat <- sc %>% left_join(key, by = "sample") %>% mutate(sid = strip(sample)) %>%
  left_join(map, by = "sid") %>%
  left_join(clin %>% dplyr::select(record_id, histology, bmi, n_hormone), by = "record_id")
msg("matched to a clinical record: %d of %d samples | with BMI: %d",
    sum(!is.na(dat$record_id)), nrow(dat), sum(is.finite(dat$bmi)))

# ---------------- category: histology first, then a score-based fallback
NEPAT <- "(?i)(\\bnepc\\b|\\bned\\b|neuroendocrine|small\\s*-?\\s*cell)"
dat <- dat %>% mutate(
  category = case_when(grepl(NEPAT, diagnosis) | grepl(NEPAT, histology) ~ "NEPC_like",
                       grepl("(?i)adeno", diagnosis) | grepl("(?i)adeno", histology) ~ "AR_altered",
                       NEPC >= median(NEPC, na.rm = TRUE) & AR < median(AR, na.rm = TRUE) ~ "NEPC_like",
                       NEPC <  median(NEPC, na.rm = TRUE) & AR >= median(AR, na.rm = TRUE) ~ "AR_altered",
                       TRUE ~ "other"),
  delta = PI16_specific - POSTN_specific)
print(table(dat$category, useNA = "ifany"))

# ---------------- the test
an <- dat %>% filter(category %in% c("NEPC_like", "AR_altered"), is.finite(bmi))
rep_n <- an %>% count(record_id) %>% filter(n > 1) %>% nrow()
msg("analysis set: %d samples (%s) | patients with >1 sample: %d", nrow(an),
    paste(names(table(an$category)), table(an$category), sep = "=", collapse = " "), rep_n)
write.csv(dat, file.path(OUT, "wcm_sample_scores.csv"), row.names = FALSE)

test_y <- function(y) {
  d <- an %>% mutate(.y = .data[[y]])
  f0 <- lm(.y ~ bmi + category, d); f1 <- lm(.y ~ bmi * category, d)
  p_int <- anova(f0, f1)$`Pr(>F)`[2]
  f1a <- lm(.y ~ bmi * category + Total_stroma, d)
  p_adj <- tryCatch(anova(lm(.y ~ bmi + category + Total_stroma, d), f1a)$`Pr(>F)`[2], error = function(e) NA)
  per <- d %>% group_by(category) %>%
    summarise(n = n(), rho = cor(.y, bmi, method = "spearman"),
              p = cor.test(.y, bmi, method = "spearman", exact = FALSE)$p.value, .groups = "drop")
  data.frame(y = y, category = per$category, n = per$n, rho = round(per$rho, 3), p = signif(per$p, 3),
             interaction_p = signif(p_int, 3), interaction_p_adj_total_stroma = signif(p_adj, 3))
}
res <- bind_rows(lapply(c("delta", "PI16_specific", "POSTN_specific", "Rest_of_stroma", "Total_stroma"), test_y))
write.csv(res, file.path(OUT, "wcm_bmi_interaction_stats.csv"), row.names = FALSE)
cat("\n=== score vs BMI, by category (primary outcome is delta) ===\n"); print(as.data.frame(res), row.names = FALSE)

if (requireNamespace("lmerTest", quietly = TRUE) && rep_n > 0) {
  m <- lmerTest::lmer(delta ~ bmi * category + (1 | record_id), data = an)
  cat("\n=== mixed model, delta ~ BMI x category + (1|record_id) ===\n"); print(coef(summary(m)))
}
cat("\n=== collinearity of the two state scores (the reason for the contrast) ===\n")
print(round(cor(an[, c("PI16_specific", "POSTN_specific", "Rest_of_stroma", "Total_stroma", "delta")],
                method = "spearman"), 2))
msg("written to %s", OUT)
