# Schematic of the compound-nomination pipeline, with the counts each step actually produced.
#
# Two branches converge on one approved-drug list. Both inputs have now been replaced relative to the original
# diagram (Figure 5/Compound screen/compound pipeline.pdf):
#
#   BRANCH A  the original step 3 intersected the promoter TF list with markers of cluster c5 in GSE235036.
#             In that dataset's whole-tissue numbering c5 is the ENTIRE stromal compartment of one mouse, so
#             its markers separate stroma from epithelium and cannot distinguish PI16+ from POSTN+; it is also
#             unrelated to atlas c5, the POSTN-high mesenchymal state. The branch now starts from the PI16+ and
#             POSTN+ per-cell states of Figures 1 and 2. The MDK-promoter HOCOMOCO scan and its chromatin
#             filter were removed earlier and do not return: 1,116 TFs have a hit in that window.
#
#   BRANCH B  ran on TCGA PRAD, SU2C and ARCHS4 adipose, none of which is neuroendocrine disease. It now runs
#             on the Beltran RR1384 metastatic cohort, NEPC-like biopsies, with the AR-altered biopsies of the
#             same cohort as a contrast.
#
# Every number on the diagram is taken from a run, not an estimate.
#
# Usage: Rscript R/fig5b_03_schematic.R
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tibble)})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
OUT <- F5_OUT
# palette comes from common.R
FILL_A <- "#EAF0F7"; FILL_B <- "#F5EEEA"; FILL_C <- "#EFEFEF"; FILL_NEW <- "#DCE8F5"; FILL_OUT <- "#F4F4F4"

box <- tribble(
  ~id, ~x0, ~x1, ~y0, ~y1, ~fill, ~border, ~title, ~body,
  # ---- branch A
  "A1", 3, 47, 84, 96, FILL_NEW, BLUE, "PI16+ and POSTN+ cells  —  Figures 1 and 2",
    "fibroblast lineage of 15 single-cell datasets, 11 human + 4 mouse\nPI16+ 10,322 cells  ·  POSTN+ 6,356  ·  rest of lineage 25,394",
  "A2", 3, 47, 68, 81, FILL_NEW, BLUE, "State transcription factors",
    "per-dataset differential expression, then recurrence:\n≥ 3 datasets of a species, or replicated across species,\nor a network hub in both, or consistent across samples\n138 TFs  ·  50 PI16+  ·  89 POSTN+  ·  37 tier 1",
  "A3", 3, 47, 57, 65, FILL_A, BLUE, "Upstream regulators",
    "kinases and receptors from OmniPath\n100 regulators of 74 state TFs",
  # ---- branch B
  "B1", 53, 97, 84, 96, FILL_B, RED, "Bulk mRNA cohort",
    "Beltran RR1384 metastatic biopsies\n34 NEPC-like  ·  71 AR-altered as contrast",
  "B2", 53, 97, 68, 81, FILL_B, RED, "Correlation with the ligands",
    "top 250 positively and negatively correlated genes\nwith MDK and with TNFSF12, within each disease group\n500 up  ·  403 down in the NEPC-like biopsies",
  "B3", 53, 97, 57, 65, FILL_B, RED, "Signature reversal  (LINCS L1000)",
    "compounds whose signature opposes the programme\n26 in NEPC-like  ·  19 of them not in AR-altered",
  # ---- converge
  "C1", 20, 80, 46, 54, FILL_C, INK, "Approved-drug retrieval",
    "Open Targets Platform  ·  234 targets queried\nkept where maximum clinical stage = approval",
  "C2", 10, 90, 33, 43, "#FFFFFF", INK, "343 compounds nominated",
    "318 FDA-approved from drug-target mapping  ·  163 act on a PI16+/POSTN+ state TF  ·  155 on an upstream kinase or receptor\n26 from signature reversal  ·  perhexiline is the only compound supported on two axes",
  "D1", 3, 47, 18, 29, FILL_NEW, BLUE, "Target-matched shortlist",
    "ESR1   fulvestrant, tamoxifen, elacestrant\nNR3C1  betamethasone, mifepristone, vamorolone\nAR     enzalutamide, abiraterone, darolutamide, +5",
  "D2", 53, 97, 18, 29, "#FFFFFF", MUTED, "Tested in Figure 5b  (n = 10)",
    "recovered, branch A: fulvestrant, tamoxifen (ESR1),\nbetamethasone (NR3C1), sorafenib, trametinib\nnot recovered: atorvastatin, mirabegron, paroxetine,\nsalbutamol, sildenafil",
  # ---- what was removed
  "X1", 3, 97, 1.5, 14, FILL_OUT, GREY, "Removed from the original pipeline",
    "GSE235036 cluster c5  —  whole-tissue numbering, in which c5 is the entire stromal compartment of one mouse dataset: 1,800 markers, 96 of them transcription\nfactors, only 15 of which survive the state analysis. It cannot separate PI16+ from POSTN+, and it is not atlas c5, the POSTN-high mesenchymal state.\nMDK-promoter HOCOMOCO motif scan and active-chromatin filter  —  1,116 TFs have a hit in that window, and TNFSF12 was never scanned the same way.\nTCGA PRAD, SU2C and ARCHS4 adipose  —  not neuroendocrine disease; the four compounds they alone nominated are no longer supported.")

arrow_seg <- tribble(~x, ~xend, ~y, ~yend,
  25, 25, 84, 81.4,   25, 25, 68, 65.4,   25, 25, 57, 54.4,
  75, 75, 84, 81.4,   75, 75, 68, 65.4,   75, 75, 57, 54.4,
  25, 40, 54.2, 54.2, 75, 60, 54.2, 54.2,
  50, 50, 46, 43.4,   25, 25, 33, 29.4,   75, 75, 33, 29.4)

p <- ggplot() +
  geom_rect(data = box, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill, colour = border),
            linewidth = .45, show.legend = FALSE) +
  scale_fill_identity() + scale_colour_identity() +
  geom_segment(data = arrow_seg, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"), linewidth = .35, colour = MUTED) +
  geom_text(data = box, aes(x = x0 + 1.6, y = y1 - 1.9, label = title),
            hjust = 0, vjust = 1, size = 2.3, fontface = "bold", colour = INK) +
  geom_text(data = box, aes(x = x0 + 1.6, y = y1 - 4.3, label = body),
            hjust = 0, vjust = 1, size = 1.85, colour = MUTED, lineheight = 1.3) +
  annotate("text", x = 25, y = 98.6, label = "BRANCH A   PI16+ / POSTN+ stromal-state transcription factors",
           size = 2.3, fontface = "bold", colour = BLUE) +
  annotate("text", x = 75, y = 98.6, label = "BRANCH B   ligand co-expression and signature reversal",
           size = 2.3, fontface = "bold", colour = RED) +
  annotate("text", x = 50, y = -3.2, hjust = .5, size = 1.8, colour = MUTED, lineheight = 1.3,
           label = paste("Substitutions, both forced by data availability: upstream kinases come from OmniPath because the full PhosphoSitePlus kinase-substrate dataset is not available locally;",
                         "approved drugs come from Open Targets because Selleck Chemicals publishes no programmatic catalogue. Every compound can be looked up in the Selleck catalogue by name.", sep = "\n")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(-6.5, 100), expand = FALSE) +
  theme_void(base_family = "Helvetica") + theme(plot.margin = margin(4, 4, 4, 4))

f <- file.path(OUT, "Compound pipeline schematic.pdf")
ggsave(f, p, width = 7.6, height = 7.2, device = cairo_pdf)
cat("written", f, "\n")
