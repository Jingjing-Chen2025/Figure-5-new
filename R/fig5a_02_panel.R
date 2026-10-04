# FIGURE 5, PANEL a - step 2 | Host adiposity shifts stromal composition in patient tumours.
#
# Print-ready vector panel, sized for the top-left slot of Figure 5 (3.6 x 2.4 in).
# Reads the per-sample scores written by R/fig5a_01_bmi_composition.R.
#
# The outcome is a CONTRAST, delta = z(PI16+ specific) - z(POSTN+ specific), because the two state scores are
# collinear in bulk and each tracks stromal content. The contrast is invariant to a shared multiplicative
# content term, so it asks whether adiposity changes the COMPOSITION of the stroma rather than its amount.
# Total stromal content does not vary with BMI in either group (shown in Extended Data).
#
# Usage: Rscript R/fig5a_02_panel.R
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
V <- f5_data("wcm_bmi_composition")

d <- read.csv(file.path(V, "wcm_sample_scores.csv")) %>%
  filter(category %in% c("NEPC_like", "AR_altered"), is.finite(bmi)) %>%
  mutate(grp = factor(ifelse(category == "NEPC_like", "NEPC-like", "AR-altered"), c("NEPC-like", "AR-altered")))
st <- read.csv(file.path(V, "wcm_bmi_interaction_stats.csv")) %>% filter(y == "delta")
p_int <- signif(st$interaction_p[1], 2); p_adj <- signif(st$interaction_p_adj_total_stroma[1], 2)
lab <- st %>% mutate(grp = ifelse(category == "NEPC_like", "NEPC-like", "AR-altered"),
                     t = sprintf("%s  rho = %+.2f, P = %.3g (n = %d)", grp, rho, p, n)) %>% arrange(grp)

p <- ggplot(d, aes(bmi, delta, colour = grp, fill = grp)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.5, alpha = 0.13) +
  geom_point(size = 1.3, stroke = 0, alpha = 0.9) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.5, size = 1.95, colour = INK,
           label = paste(lab$t, collapse = "\n")) +
  annotate("text", x = -Inf, y = -Inf, hjust = -0.06, vjust = -0.7, size = 1.95, colour = MUTED,
           label = sprintf("interaction P = %s; adjusted for total stromal score, P = %s", p_int, p_adj)) +
  scale_colour_manual(values = c(`NEPC-like` = RED, `AR-altered` = BLUE), name = NULL,
                      aesthetics = c("colour", "fill")) +
  labs(title = "a  Host adiposity shifts stromal composition",
       subtitle = "PI16+ minus POSTN+ score; invariant to total stromal content",
       x = expression("BMI (kg m"^-2*")"), y = "PI16+ − POSTN+  (z)") +
  theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.3, colour = MUTED),
        axis.ticks = element_line(linewidth = 0.3, colour = MUTED),
        axis.text = element_text(colour = INK), axis.title = element_text(size = 6.4, colour = MUTED),
        plot.title = element_text(face = "bold", size = 8), plot.subtitle = element_text(size = 5.8, colour = MUTED),
        legend.position = "none", plot.title.position = "plot", plot.margin = margin(4, 6, 3, 3))

f <- file.path(F5_OUT, "Figure 5a - BMI stromal composition.pdf")
ggsave(f, p, width = 3.6, height = 2.4, device = cairo_pdf)
cat("written", f, "\n")
