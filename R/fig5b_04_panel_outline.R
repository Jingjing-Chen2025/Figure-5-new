# Figure 5, new panel b | Outline of the compound-nomination pipeline.
#
# Sits between the BMI panel (a) and the compound screen (now c), so the reader sees where the tested compounds
# came from before seeing what they did. A condensed version of the full schematic in
# Figure 5/Compound screen/Compound pipeline schematic.pdf, which keeps the method detail and the caveats.
#
# Drawn as a horizontal strip rather than a vertical flow: at main-figure size a vertical diagram would have to
# be squeezed into roughly 190 pt of width, where the text becomes unreadable. Full width keeps 7-point type.
#
# Every count is from a run (R/fig5b_02_compound_pipeline.R): the PI16+ and POSTN+ per-cell states of
# Figures 1 and 2 across 15 datasets; 138 state transcription factors; 100 upstream kinases and receptors;
# 318 approved compounds, 163 of them against a state TF; 26 further compounds from signature reversal in the
# Beltran NEPC-like biopsies; 10 tested, of which 5 are recovered.
#
# Usage: Rscript R/fig5b_04_panel_outline.R
suppressPackageStartupMessages({library(ggplot2); library(tibble)})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
OUT <- F5_OUT
# palette comes from common.R
FILL_A <- "#EAF0F7"; FILL_B <- "#F5EEEA"; FILL_HI <- "#DCE8F5"; FILL_N <- "#FFFFFF"

b <- tribble(
  ~x0,  ~x1, ~y0, ~y1, ~fill,   ~col,  ~head,                    ~body,
   1.0, 23.0, 42,  84, FILL_HI, BLUE,  "PI16+ / POSTN+ cells",   "per-cell states of Fig. 1 and 2,\n15 datasets, 11 human + 4 mouse\n10,322 PI16+  ·  6,356 POSTN+",
  26.5, 46.0, 42,  84, FILL_HI, BLUE,  "State TFs",              "recurrent across datasets\nand species\n138 TFs  ·  50 PI16+  ·  89 POSTN+",
  49.5, 64.0, 46,  80, FILL_A,  BLUE,  "Upstream regulators",    "kinases and receptors\n100 for 74 TFs",
  67.5, 84.0, 42,  84, FILL_N,  INK,   "Approved drugs",         "Open Targets\n318 FDA-approved\n163 against a state TF",
  87.5, 99.0, 46,  80, FILL_N,  MUTED, "Screened (c)",           "10 compounds\n5 from this pipeline",
   1.0, 64.0,  6,  30, FILL_B,  RED,   "Second input",           "Beltran NEPC-like biopsies  →  top 250 genes correlated with MDK and TNFSF12  →  LINCS L1000\ncompounds opposing that programme: 26, of which 19 do not oppose it in AR-altered biopsies")

arr <- tribble(~x, ~xend, ~y, ~yend,
  23.0, 25.9, 63, 63,   46.0, 48.9, 63, 63,   64.0, 66.9, 63, 63,   84.0, 86.9, 63, 63,
  64.0, 72.0, 18, 18,   72.0, 72.0, 18, 41.4)

p <- ggplot() +
  geom_rect(data = b, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill, colour = col),
            linewidth = .4, show.legend = FALSE) +
  scale_fill_identity() + scale_colour_identity() +
  geom_segment(data = arr, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(1.4, "mm"), type = "closed"), linewidth = .32, colour = MUTED) +
  geom_text(data = b, aes(x = x0 + 1.1, y = y1 - 4.5, label = head),
            hjust = 0, vjust = 1, size = 2.1, fontface = "bold", colour = INK) +
  geom_text(data = b, aes(x = x0 + 1.1, y = y1 - 12, label = body),
            hjust = 0, vjust = 1, size = 1.8, colour = MUTED, lineheight = 1.3) +
  annotate("text", x = 1, y = 95, hjust = 0, size = 2.6, fontface = "bold", colour = INK,
           label = "b  How the screened compounds were nominated") +
  annotate("text", x = 1, y = 88.5, hjust = 0, size = 1.85, colour = MUTED, lineheight = 1.25,
           label = "Transcription factors of the PI16+ and POSTN+ stromal states, with their upstream kinases and receptors, were used to retrieve approved drugs.") +
  coord_cartesian(xlim = c(0, 100), ylim = c(2, 99), expand = FALSE) +
  theme_void(base_family = "Helvetica") + theme(plot.margin = margin(2, 2, 2, 2))

f <- file.path(OUT, "Figure 5b - compound pipeline outline.pdf")
ggsave(f, p, width = 7.3, height = 1.95, device = cairo_pdf)
cat("written", f, "\n")
