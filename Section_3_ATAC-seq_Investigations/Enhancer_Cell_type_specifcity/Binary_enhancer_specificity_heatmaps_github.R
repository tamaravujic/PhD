# ==============================================================================
# Script: binary_enhancer_specificity_heatmaps.R
#
# Purpose:
#   Visualises the presence/absence of ATAC-seq peaks at putative cis-regulatory
#   elements (CREs/enhancers) across different endothelial cell types. Creates
#   comparative heatmaps showing enhancer accessibility patterns in:
#     (1) In-house endothelial cell (EC) ATAC-seq data
#     (2) ENCODE project endothelial cell ATAC-seq data (averaged signal)
#
# Workflow:
#   1. Load in-house and ENCODE public EC ATAC-seq peak data
#   2. Transform ENCODE data to binary (peak present/absent) using threshold
#   3. Reshape data for visualization (melting to long format)
#   4. Generate heatmap for in-house EC data
#   5. Generate heatmap for ENCODE data
#   6. Export both heatmaps as high-resolution PDFs
#   7. Save R session information for reproducibility
#
# Input Data:
#   LEC_specificity.xlsx
#     - Excel workbook with two sheets:
#       Sheet2 (ECs): In-house ATAC-seq peak data
#         Columns: Enhancer_region 
#         row names: cell_type_1, cell_type_2, ...
#         Values: Numeric accessibility/peak intensity scores (0-1 scale)
#       Sheet3 (ENCODE): Public ENCODE ATAC-seq data
#         Columns: Enhancer_region 
#         row names: cell_type_1, cell_type_2, ...
#         Values: Numeric peak intensity scores (0-1 scale)
#
# Outputs:
#   ECs_enhancer_specificity_plot.pdf
#     - Heatmap showing ATAC peak accessibility across in-house EC types
#   ENCODE_enhancer_specificity_plot_averages.pdf
#     - Heatmap showing ATAC peak accessibility across ENCODE EC types
#   Enhancer_specificity_heatmaps_session_info
#     - Text file containing R version, package versions, and system info
#       (useful for reproducibility and debugging)
#
# Notw
#   - Input data are numeric accessibility scores (0-1 scale, where 1 = peak present)
#   - ENCODE data represents averaged signal across replicates
#   - A threshold of >= 0.5 is used to convert ENCODE continuous scores to binary
# ==============================================================================

# ==============================================================================
# SECTION 1: ENVIRONMENT SETUP AND DEPENDENCIES
# ==============================================================================
# Set working directory to the folder containing input data files.
setwd("/Users/tamaravujic/Desktop/BHF 1+3/Enhancer_specificty")

# Load required R packages.
library(readxl)
library(ggplot2)
library(tidyr)
library(gplots)

# ==============================================================================
# SECTION 2: LOAD IN-HOUSE ENDOTHELIAL CELL ATAC-SEQ DATA
# ==============================================================================
# Define sheet names for readability. These correspond to specific sheets
# in the LEC_specificity.xlsx workbook.
ECs <- "ECs"
ENCODE <- "ENNCODE"

# Load in-house ATAC-seq peak data f
enhancers <- read_excel("/Users/tamaravujic/Desktop/BHF 1+3/Enhancer_specificty/LEC_specificity.xlsx", sheet = ECs)

# Load ENCODE public ATAC-seq data for comparison (standardized, large-scale dataset)
ENCODE <- read_excel("/Users/tamaravujic/Desktop/BHF 1+3/Enhancer_specificty/LEC_specificity.xlsx", sheet = ENCODE)

# ==============================================================================
# SECTION 3: TRANSFORM ENCODE DATA TO BINARY (PEAK PRESENT/ABSENT)
# ==============================================================================
# Convert ENCODE continuous accessibility scores to binary (0/1) classification.
# a threshold of >= 0.5 is applied to convert to binary "peak present" (1) or 
# "peak absent" (0) 
#
# Operation: Apply the threshold to all columns EXCEPT the first (which contains row labels)
ENCODE[, -1] <- apply(ENCODE[, -1], 2, function(x) ifelse(x >= 0.5, 1, 0))

is.data.frame(ENCODE)

# ==============================================================================
# SECTION 4: RESHAPE DATA FOR VISUALIZATION
# ==============================================================================
# The ggplot2 heatmap requires data in "long" format (one row per observation)
# rather than "wide" format (multiple columns per cell type). 

# Reshape in-house EC data to long format
ECs_data_melted <- melt(enhancers)

# Reshape ENCODE data to long format
ENCODE_data_melted <- melt(ENCODE)

# ==============================================================================
# SECTION 5: DEFINE LEGEND PARAMETERS
# ==============================================================================
# Define the legend and labels for the heatmap color scale.
breaks <- c(0, 1)  
labels <- c("NO", "YES")

# ==============================================================================
# SECTION 6: CREATE AND PLOT IN-HOUSE EC HEATMAP
# ==============================================================================
# Plot structure:
#   - X-axis: cell types (different EC types)
#   - Y-axis: enhancer regions (putative CREs)
#   - Fill color: accessibility (white = no peak, pink = peak present)

ECs_heatmap <- ggplot(ECs_data_melted, aes(x = variable, y = Enhancer_region, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "pink", name = "ATAC peak present", 
                      breaks = breaks, labels = labels) +
  labs(title = "Presence of ATAC peaks at putative CRE regions in different endothelial cell types",
       x = "Cell Type",
       y = "Enhnacer Region") +
    theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 40)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, vjust = 3, hjust = 0.8),
    legend.box.spacing = unit(1, "cm"),
    panel.background = element_rect(color = "black", size = 2.5),
    plot.title = element_text(hjust = 0.6)
  ) +
    guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 0.5)) +
  coord_fixed()

ECs_heatmap

# Save the high-resolution heatmap as PDF for publication
ggsave("ECs_enhancer_specificity_plot.pdf", ECs_heatmap, width = 10, height = 10)

# ==============================================================================
# SECTION 7: CREATE AND PLOT ENCODE HEATMAP
# ==============================================================================
# Generate a heatmap using the ENCODE public dataset for comparison with
# in-house data. Uses identical styling for consistent visualization.

ENCODE_heatmap <- ggplot(ENCODE_data_melted, aes(x = variable, y = Enhancer_region, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "pink", name = "ATAC peak present", 
                      breaks = breaks, labels = labels) +
  labs(title = "Presence of ATAC peaks at putative CRE regions in different endothelial cell types",
       x = "Cell Type",
       y = "Enhnacer Region") +
    theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 40)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, vjust = 3, hjust = 0.8),
    legend.box.spacing = unit(1, "cm"),
    panel.background = element_rect(color = "black", size = 2.5),
    plot.title = element_text(hjust = 0.6)
  ) +
    guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 0.5)) +
    coord_fixed()

# Display the plot in the R graphics window
ENCODE_heatmap

# Save the high-resolution heatmap as PDF for publication
ggsave("ENCODE_enhancer_specificity_plot_averages.pdf", ENCODE_heatmap, width = 25, height = 25)

# ==============================================================================
# SECTION 8: SAVE SESSION INFORMATION FOR REPRODUCIBILITY
# ==============================================================================
# Export complete R session information (version, attached packages, loaded packages,
# package versions, system information, etc.) to a text file.

writeLines(capture.output(sessionInfo()), "Enhancer_specificity_heatmaps_session_info")
