# ==============================================================================
# Script: chrombpnet_variant_scores_analysis.R
#
# Purpose:
#   Average ChromBPNet variant prediction scores across multiple folds.
#   Generate volcano plots, statistical comparisons, and identify variants
#   with significant predicted regulatory effects.
#
# Workflow:
#   1. Load variant scores from 5 fold predictions
#   2. Calculate mean scores across folds
#   3. Generate volcano plot (logFC vs p-value)
#   4. Identify and export significant variants
#   5. Compare regulatory effects between cohorts (GEL vs Common)
#
# Input:
#   results_fold_0/5: ChromBPNet output directories with variant_scores.tsv files
#   mean_variant_scores.tsv: averaged scores across folds
#
# Output:
#   mean_variant_scores.tsv: averaged scores from all 5 folds
#   significant_variants.tsv: variants passing p-value and logFC thresholds
#   Figure1_volcano_accessibility_p0.01_fc0.5.pdf: volcano plot visualisation
#   Figure2_GEL_vs_Common_enrichment.pdf: cohort comparison plots
#   Figure3_Top_GEL_variants_directional.pdf: top variants by effect direction
#
# Notes:
#   - Variant scores are averaged across 5-fold cross-validation
#   - Thresholds: p < 0.01, |logFC| >= 0.5
#   - Cohort analysis requires "Cohort" column in data
# ==============================================================================

# ==============================================================================
# SECTION 1: LOAD PACKAGES
# ==============================================================================
library(tidyverse)
library(ggrepel)
library(patchwork)
library(dplyr)
library(ggplot2)

# ==============================================================================
# SECTION 2: AVERAGE VARIANT SCORES ACROSS 5 FOLDS
# ==============================================================================
# After running variant_prediction.sh for all 5 folds, average the logFC and
# p-values to get a single prediction per variant

base_dir <- "/rds/general/user/tv722/projects/tamara_vujic_phd/live/ChromBPNet/chrombpnet0.5"
name <- "chr21_dbSnp155Common_lymphoedema"

# File naming pattern: name_variant_scores.variant_scores.tsv
file_stem <- paste0(name, "_variant_scores.variant_scores.tsv")

# Load scores from all 5 folds
scores <- list()
for (fold in 0:4) {
  f <- file.path(base_dir, paste0("results_fold_", fold), file_stem)
  if (!file.exists(f)) stop("Missing file: ", f)
  scores[[paste0("fold_", fold)]] <- read.table(f, header = TRUE, sep = "\t")
}

# Keep only variants present in all 5 folds
common_ids <- Reduce(intersect, lapply(scores, \(x) x$variant_id))

# Start with fold_0, filter to common variants
base <- scores[["fold_0"]] %>%
  filter(variant_id %in% common_ids) %>%
  arrange(variant_id)

# Calculate mean scores across all folds for each variant
means <- base %>%
  select(chr, pos, allele1, allele2, variant_id) %>%
  mutate(
    mean_logfc = rowMeans(sapply(scores, \(x)
                                 x %>% filter(variant_id %in% common_ids) %>% 
                                 arrange(variant_id) %>% pull(logfc)
    )),
    mean_logfc_pval = rowMeans(sapply(scores, \(x)
                                      x %>% filter(variant_id %in% common_ids) %>% 
                                      arrange(variant_id) %>% pull(logfc.pval)
    )),
    mean_allele1_pred_counts = rowMeans(sapply(scores, \(x)
                                               x %>% filter(variant_id %in% common_ids) %>% 
                                               arrange(variant_id) %>% pull(allele1_pred_counts)
    )),
    mean_allele2_pred_counts = rowMeans(sapply(scores, \(x)
                                               x %>% filter(variant_id %in% common_ids) %>% 
                                               arrange(variant_id) %>% pull(allele2_pred_counts)
    ))
  )

# Save averaged scores
out_file <- file.path(base_dir, "mean_variant_scores.tsv")
write.table(means, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved averaged scores to:", out_file, "\n")

# ==============================================================================
# SECTION 3: GENERATE VOLCANO PLOT
# ==============================================================================
# Volcano plot: visualise variants with significant predicted regulatory effects
# X-axis: log fold-change (logFC)
# Y-axis: -log10(p-value)
# Color: accessibility change direction

setwd(base_dir)
m <- read.table("mean_variant_scores.tsv", header = TRUE, sep = "\t")

# Define significance thresholds, can change as desired
p_cutoff  <- 0.01
fc_cutoff <- 0.5

# Categorise variants by threshold criteria
m <- m %>%
  mutate(
    y = -log10(mean_logfc_pval),
    category = case_when(
      mean_logfc_pval < p_cutoff & mean_logfc <= -fc_cutoff ~ "Decreased accessibility",
      mean_logfc_pval < p_cutoff & mean_logfc >=  fc_cutoff ~ "Increased accessibility",
      TRUE ~ "Insignificant"
    )
  )

# Select top 12 most significant variants to label on plot
top_to_label <- m %>%
  filter(category != "Insignificant") %>%
  arrange(mean_logfc_pval) %>%
  head(12)

# Force symmetric x-axis for better visualization
xmax <- max(abs(m$mean_logfc), na.rm = TRUE)
xlim_sym <- ceiling(xmax * 10) / 10  # round up to 0.1

# Header panel: adds title and thresholds above plot
header <- ggplot() +
  annotate("text", x = -0.6, y = 0.75,
           label = "Decreased accessibility",
           color = "dodgerblue3", fontface = "bold", size = 5.5, hjust = 0.5) +
  annotate("text", x =  0.6, y = 0.75,
           label = "Increased accessibility",
           color = "firebrick2", fontface = "bold", size = 5.5, hjust = 0.5) +
  annotate("text", x = 0, y = 0.12,
           label = "P value < 0.01 and logFC >= 0.5",
           size = 3, hjust = -0.5) +
  xlim(-1, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(t = 10, r = 20, b = 0, l = 20))

# Main volcano plot
volcano <- ggplot() +
  # Insignificant variants (grey background)
  geom_point(data = filter(m, category == "Insignificant"),
             aes(mean_logfc, y),
             color = "grey80", alpha = 0.8, size = 2.6) +
  # Decreased accessibility (blue)
  geom_point(data = filter(m, category == "Decreased accessibility"),
             aes(mean_logfc, y),
             color = "dodgerblue3", alpha = 0.95, size = 3.0) +
  # Increased accessibility (red)
  geom_point(data = filter(m, category == "Increased accessibility"),
             aes(mean_logfc, y),
             color = "firebrick2", alpha = 0.95, size = 3.0) +
  # Threshold lines: vertical for logFC, horizontal for p-value
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.6) +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", linewidth = 0.6) +
  # Label top significant variants
  ggrepel::geom_text_repel(
    data = top_to_label,
    aes(mean_logfc, y, label = variant_id),
    size = 4.2,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  labs(
    x = "LogFC",
    y = expression(-log[10]("P value"))
  ) +
  coord_cartesian(
    xlim = c(-xlim_sym, xlim_sym),
    ylim = c(0, 4),
    clip = "off"
  ) +
  scale_y_continuous(breaks = 0:4) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 14),
    plot.margin = margin(t = 0, r = 20, b = 15, l = 20)
  )

# Combine header + volcano plot
p_final <- header / volcano + plot_layout(heights = c(1.2, 10))
print(p_final)

# Save volcano plot
ggsave("Figures1_volcano_accessibility_p0.01_fc0.5.pdf", plot = p_final, width = 8.5, height = 9.5)

# ==============================================================================
# SECTION 4: EXPORT SIGNIFICANT VARIANTS
# ==============================================================================
# Save table of variants meeting significance thresholds
# Change thresholds as desired

m <- read.table("mean_variant_scores.tsv", header = TRUE, sep = "\t")

sig <- m %>%
  filter(mean_logfc_pval < 0.01, abs(mean_logfc) >= 0.5) %>%
  arrange(mean_logfc_pval, desc(abs(mean_logfc)))

write.table(sig,
            file = "significant_variants.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("Saved", nrow(sig), "significant variants to significant_variants.tsv\n")

# ==============================================================================
# SECTION 5: COMPARE COHORTS (GEL vs COMMON VARIANTS)
# ==============================================================================
# If data contains multiple cohorts, compare regulatory impact distributions
# NOTE: Requires "Cohort" column in mean_variant_scores.tsv

setwd('~/Desktop/ChrombpNET')
m <- read.table("mean_variant_scores.tsv", header = TRUE, sep = "\t")

# Statistical test: Wilcoxon rank-sum test
w <- wilcox.test(abs(mean_logfc) ~ Cohort, data = m)

# Visualisation: Violin plot with overlaid boxplot
p_fig2 <- ggplot(m, aes(x = Cohort, y = abs(mean_logfc), fill = Cohort)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = c("Common" = "grey75", "GEL" = "mediumpurple3")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 16) +
  labs(
    y = "Mean logFC",
    x = "",
    title = "GEL variants show stronger predicted regulatory effects"
  ) +
  annotate("text", x = 1.5, y = 0.28,
           label = paste0("Wilcoxon p = ", signif(w$p.value, 3)),
           size = 5)

print(p_fig2)
ggsave("Figure2_GEL_vs_Common_enrichment.pdf", p_fig2, width = 6, height = 6)

# ==============================================================================
# SECTION 6: VISUALISE TOP VARIANTS BY DIRECTION
# ==============================================================================
# Create barplot showing top variants with increased vs decreased accessibility

topN <- 10

gel <- m %>%
  filter(Cohort == "GEL") %>%
  mutate(direction = ifelse(mean_logfc > 0,
                            "Increased accessibility",
                            "Decreased accessibility"),
         abs_logfc = abs(mean_logfc))

# Select top N variants in each direction
gel_top <- bind_rows(
  gel %>% filter(direction == "Increased accessibility") %>%
    arrange(desc(mean_logfc)) %>% slice_head(n = topN),
  gel %>% filter(direction == "Decreased accessibility") %>%
    arrange(mean_logfc) %>% slice_head(n = topN)
)

# Barplot of top variants
p_fig3 <- ggplot(gel_top,
                 aes(x = reorder(variant_id, mean_logfc),
                     y = mean_logfc,
                     fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Increased accessibility" = "firebrick2",
                               "Decreased accessibility" = "dodgerblue3")) +
  theme_classic(base_size = 16) +
  labs(
    x = "GEL variants",
    y = "Mean logFC",
    title = "Top GEL variants by direction of regulatory effect",
    fill = ""
  )

print(p_fig3)
ggsave("Figure3_Top_GEL_variants_directional.pdf", p_fig3, width = 7, height = 7)

cat("Analysis complete! All figures saved.\n")
