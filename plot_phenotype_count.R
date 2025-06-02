rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

library(ggplot2)
library(reshape2)
library(readxl)
library(RColorBrewer)
library(dplyr)

# Transformation functions (50x expansion for 0-0.1%)
stretch_trans <- function(x) {
  ifelse(x <= 0.1, 
         x * 50,               # 0-0.1% → 0-5 (50x visual height)
         5 + (x-0.1)*(95/99.9) # 0.1-100% → 5-100 (linear)
  )
}

inverse_trans <- function(x) {
  ifelse(x <= 5, 
         x / 50,               # Reverse 0-0.1% scaling
         0.1 + (x-5)*(99.9/95) # Reverse 0.1-100% scaling
  )
}

if (!dir.exists("results/phenotype_distribution")) {
  dir.create("results/phenotype_distribution", recursive = TRUE)
}

df <- read_excel("Phenotype_count.xlsx")

df_long <- melt(df, id.vars = c("Blood_Group", "Phenotype", "Type"), 
            variable.name = "Population", value.name = "Value")
df_long$Value <- as.numeric(df_long$Value)
       
# Create a new x-axis label combining Population and Type
df_long$Pop_Type <- paste0(df_long$Population, "_", df_long$Type)

# Define the custom order
custom_order <- c("AFR_Ery", "AFR_Inf", "AMR_Ery", "AMR_Inf", 
                 "EAS_Ery", "EAS_Inf", "EUR_Ery", "EUR_Inf", 
                 "SAS_Ery", "SAS_Inf", "ALL_Ery", "ALL_Inf")

# Convert Pop_Type into a factor with the specified order
df_long$Pop_Type <- factor(df_long$Pop_Type, levels = custom_order)

for (group in unique(df_long$Blood_Group)) {
    df_group <- df_long[df_long$Blood_Group == group, ]

    # Reorder phenotype appearance in the plot
    phenotype_order <- unique(df_group$Phenotype)
    df_group$Phenotype <- factor(df_group$Phenotype, levels = rev(phenotype_order))

    # Define the colours for each phenotype
    colours <- colorRampPalette(c("#004c6d", "#7cb5ec", "#b9ffff"))(length(levels(df_group$Phenotype)))

    plot <- ggplot(df_group, aes(x = Pop_Type, y = Value, fill = Phenotype)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8,
                color = "white", linewidth = 0.5
        ) +
    
        # Bold group separators
        geom_vline(xintercept = seq(2.5, by = 2, length.out = length(custom_order)/2 - 1),
            color = "gray60", linewidth = 0.6, linetype = "dotdash"
        ) +
    
        # Custom y-axis scaling
        scale_y_continuous(
            trans = trans_new("stretch", stretch_trans, inverse_trans),
            breaks = c(0, 0.1, 25, 50, 75, 100), labels = c(0, 0.1, 25, 50, 75, 100),
            expand = c(0, 0), name = "Percentage (%)"
        ) +
    
        # Color and faceting
        scale_fill_manual(values = colours) +
        facet_wrap(~ Blood_Group, scales = "free_x") +
    
        # Labels
        labs(x = "Population Group and Data Source", fill = "Phenotype") +
    

        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
            axis.text.y = element_text(size = 15),
            axis.title = element_text(size = 15, face = "bold"),
            legend.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 15),
            strip.text = element_text(size = 15, face = "bold"),
            plot.margin = unit(c(1, 1, 1, 1), "cm"),
            legend.position = "bottom",
        
            panel.grid.major.y = element_line(
                color = "gray70",
                linewidth = 0.4,
                linetype = "dotted"
            ),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        guides(fill = guide_legend(ncol = 4, keywidth = 0.8, keyheight = 0.8))

    group <- gsub("/", "_", group) # Handle CH/RG Blood Group
    ggsave(paste0("results/phenotype_distribution/bar_chart_", group, ".png"), plot = plot, width = 40, height = 20, units = "cm", dpi = 300)
}