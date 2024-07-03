setwd("/Users/tamaravujic/Desktop/BHF 1+3/Enhancer_specificty")

## Load packages (Install if needed)
library(readxl)
library(ggplot2)
library(tidyr)
library(gplots)


###### IN HOUSE ENDOTHELIAL CELL ATAC SEQEUNCING ###########
#Load data
ECs <- "Sheet2"
ENCODE <- "Sheet3"
enhancers <- read_excel("/Users/tamaravujic/Desktop/BHF 1+3/Enhancer_specificty/LEC_specificity.xlsx", sheet = ECs)
ENCODE <- read_excel("/Users/tamaravujic/Desktop/BHF 1+3/Enhancer_specificty/LEC_specificity.xlsx", sheet = ENCODE)

# Replace values in the dataframe
ENCODE[, -1] <- apply(ENCODE[, -1], 2, function(x) ifelse(x >= 0.5, 1, 0))
is.data.frame(ENCODE)

#put data into correct format for a heatmap
ECs_data_melted <- melt(enhancers)
ENCODE_data_melted <- melt(ENCODE)


# Labels for the legend
breaks <- c(0, 1)  
labels <- c("NO", "YES") 

#ECs Plot heatmap
ECs_heatmap <- ggplot(ECs_data_melted, aes(x = variable, y = Enhancer_region, fill = value)) +
  geom_tile(color = "black") +  # Add black border lines
  scale_fill_gradient(low = "white", high = "pink", name = "ATAC peak present", breaks = breaks, labels = labels) +
  labs(title = "Presence of ATAC peaks at putative CRE regions in different endothelial cell types",
       x = "Cell Type",
       y = "Enhnacer Region") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_text(margin = margin(r = 40)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, vjust = 3, hjust=0.8),
        legend.box.spacing = unit(1, "cm"),
        panel.background = element_rect(color = "black", size = 2.5),
        plot.title = element_text(hjust = 0.6))+
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 0.5)) +  # Add border around the gradient legend
  coord_fixed()

ECs_heatmap


# Save the plot
ggsave("ECs_enhancer_specificity_plot.pdf", ECs_heatmap, width = 10, height = 10)

### ENCODE heatmap
ENCODE_heatmap <- ggplot(ENCODE_data_melted, aes(x = variable, y = Enhancer_region, fill = value)) +
  geom_tile(color = "black") +  # Add black border lines
  scale_fill_gradient(low = "white", high = "pink", name = "ATAC peak present", breaks = breaks, labels = labels) +
  labs(title = "Presence of ATAC peaks at putative CRE regions in different endothelial cell types",
       x = "Cell Type",
       y = "Enhnacer Region") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_text(margin = margin(r = 40)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, vjust = 3, hjust=0.8),
        legend.box.spacing = unit(1, "cm"),
        panel.background = element_rect(color = "black", size = 2.5),
        plot.title = element_text(hjust = 0.6))+
  guides(fill = guide_colorbar(frame.colour = "black", frame.linewidth = 0.5)) # Add border around the gradient legend


ENCODE_heatmap


# Save the plot
ggsave("ENCODE_enhancer_specificity_plot_averages.pdf", ENCODE_heatmap, width = 25, height = 25)

# A txt file with the session info and package versions - will be needed when time to write 
writeLines(capture.output(sessionInfo()), "Enhancer_specificity_heatmaps_session_info")
