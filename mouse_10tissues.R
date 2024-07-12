library(ggplot2)
library(tidyverse)
#install.packages("readxl")
library(readxl)
library(ggpubr)
setwd("/Users/liweisha/Documents/lipid lengthening with age/human/plasma")  
# Sample data (replace this with your own dataset)
# Example dataset structure:
# lipid_class: Categorized lipid class
# chain_length: Numeric chain length values
# fold_change: Numeric fold change values
data<-Copy_of_BMP_sup_table_1

library(dplyr)
grouped_data <- data %>%
  group_by(Tissues, data$chain_length) %>%
  summarise(Count = n()) %>%
  ungroup()
ggplot(grouped_data, aes(x = grouped_data$`data$chain_length`, y = Count, fill = as.factor(`data$chain_length`))) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Tissues, scales = "free_y", ncol = 1) +
  labs(x = "Chain Length", y = "Count", title = "Frequency of Chain Length in Different Tissues") +
  theme_minimal() +
  theme(legend.position = "none") # Remove the legend since chain length is clear from the x-axis

data$AverageAbundance_young <- rowMeans(data[, 6:13])
data$AverageAbundance_aged <- rowMeans(data[, c(14,16:21)], na.rm = TRUE)


ggplot(data, aes(x = chain_length, y =AverageAbundance_young)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Tissues, scales = "free_y", ncol = 1) +
  labs(x = "Chain Length", y = "Count", title = "abundance of Chain Length in Different Tissues_young") +
  theme_minimal() +
  theme(legend.position = "none") # Remove the legend since chain length is clear from the x-axis

lipids <- data$Lipid_species
# Extract chain length and double bonds using regular expressions
matches <- regmatches(lipids, regexec("(\\d+):(\\d+)", lipids))

# Create a data frame
lipid_data <- data.frame(
  Lipid_species = lipids,
  ChainLength = as.integer(sapply(matches, function(x) x[2])),
  DoubleBonds = as.integer(sapply(matches, function(x) x[3]))
)
data<- full_join(data, lipid_data, by = "Lipid_species")
data$lipid_class <- gsub("[0-9()]", "", data$Lipid_species)

# Clean up any leading or trailing whitespace
data$lipid_class <- trimws(data$lipid_class)
data$lipid_class <- gsub(":", "", data$lipid_class)


# Function to create the plots
create_plots <- function(data) {
  unique_lipid_classes <- unique(data$lipid_class)
  unique_tissues <- unique(data$Tissues)
  
  plot_list <- list()
  
  for (lipid_class in unique_lipid_classes) {
    for (tissue in unique_tissues) {
      plot_data <- data %>%
        filter(lipid_class == !!lipid_class, tissue == !!tissue)
      
      p <- ggplot(plot_data, aes(x = ChainLength, y = DoubleBonds)) +
        geom_point(aes(color = `logFC (old_vs_young)`, size = -log10(`adj.P.Val (old_vs_young)`))) +
        scale_color_gradientn(values = seq(0,1,0.1), colors= c("blue", "white", "red"))+
        labs(title = paste("Lipid Class:", lipid_class, "in", tissue),
             x = "Chain Length", y = "Double Bonds",
             color = "LogFC", size = "Adjusted P-value (-log10)") +
        theme_minimal()
      
      plot_list[[paste(lipid_class, tissue, sep = "_")]] <- p
    }
  }
  
  return(plot_list)
}

# Create the plots
plots <- create_plots(data)

# Display the plots (as an example, we display the first one)
print(plots[[1]])

# Save plots as files (optional)
for (name in names(plots)) {
  ggsave(filename = paste0(name, ".png"), plot = plots[[name]], width = 8, height = 6)
}


plots <- list()
Tissues<-data$Tissues
lipid_class<-data$lipid_class
# Iterate over each tissue
for (tissue in Tissues) {
  # Filter data for the current tissue
  tissue_data <- data %>% filter(Tissues == tissue)
  
  # Create a list to store ggplot objects for each lipid class in the current tissue
  lipid_plots <- list()
  
  # Iterate over each lipid class
  for (lipid_class in lipid_class) {
    # Filter data for the current lipid class
    lipid_class_data <- tissue_data %>% filter(lipid_class == lipid_class)
    
    # Create a ggplot object for the current lipid class and tissue
    lipid_plot <- ggplot(lipid_class_data, aes(x = ChainLength, y = DoubleBonds, color = `logFC (old_vs_young)`, size = log10(`adj.P.Val (old_vs_young)`))) +
      geom_point() +
      scale_color_gradient(name = "logFC (old_vs_young)") +
      scale_size_continuous(name = "adj.P.Val (old_vs_young)") +
      labs(title = paste("Tissue:", tissue, ", Lipid Class:", lipid_class),
           x = "Chain Length", y = "Double Bonds") +
      theme_minimal()
    
    # Add the ggplot object to the list
    lipid_plots[[lipid_class]] <- lipid_plot
  }
  
  # Combine ggplot objects for each lipid class into a single plot for the current tissue
  combined_plot <- patchwork::wrap_plots(lipid_plots, ncol = length(lipid_classes))
  
  # Add the combined plot to the list
  plots[[tissue]] <- combined_plot
}

# Combine ggplot objects for each tissue into a single plot
final_plot <- patchwork::wrap_plots(plots, nrow = length(tissues))

# Display the final plot
final_plot
print(plots[[1]])



data1<-data[data$Tissues=="BAT",]
ggplot(data1, aes(data1$chain_length ,data1$double_bonds)) +
  geom_point( aes(color = `logFC (old_vs_young)`, 
                  size = -log10(`adj.P.Val (old_vs_young)`))) +
  scale_color_gradientn(values = seq(0,1,0.1), colors= c("blue", "white", "red"))+
  facet_wrap(~data1$LipidType,scales = "free")

data1<-data[data$Tissues=="Brain",]
data1<-data



data1$ddatadata1$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
data1$diffexpressed[data1$`logFC (old_vs_young)` > 1 & data1$`adj.P.Val (old_vs_young)` < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data1$diffexpressed[data1$`logFC (old_vs_young)` < -1 & data1$`adj.P.Val (old_vs_young)`< 0.05] <- "DOWN"
# Convert directly in the aes()
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
data1$delabel <- NA
data1$delabel[data1$diffexpressed != "NO"] <- data1$Lipid_species[data1$diffexpressed != "NO"]

ggplot(data1, aes(x=`logFC (old_vs_young)`, y=-log10(`adj.P.Val (old_vs_young)`), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data1, aes(x=`logFC (old_vs_young)`, y=-log10(`adj.P.Val (old_vs_young)`), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
 
library(smplot2)
library(sjPlot)
sp <- ggplot(data1, aes(ChainLength, data1$`logFC (old_vs_young)`)) +
  geom_point()+
  sm_statCorr(corr_method = 'spearman')+
  facet_wrap(~lipid_class,scales = "free")

library(sjPlot)

# Calculate correlation statistics using sm_statCorr
correlation_info <- sm_statCorr(data1, corr_method = 'spearman')

# Extract correlation coefficient and p-value for each lipid class
correlation_values <- correlation_info[, c("cor", "p")]
correlation_info$correlation
# Print the correlation values
print(correlation_values)

sp <- sp + geom_smooth(method = "lm", color = ifelse(data$pvalue < 0.05, "red", "black"))

# Display the plot
print(sp)

?sm_statCorr


# Load required libraries
library(ggplot2)

# Define the plot
sp <- ggplot(data1, aes(ChainLength, `logFC (old_vs_young)`)) +
  geom_point() +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  facet_wrap(~lipid_class, scales = "free") +
  theme_minimal()

# Calculate correlation for each lipid class
for (lipid_class in unique(data$lipid_class)) {
  # Filter data for the current lipid class and remove points with ChainLength below 3
  correlation_data <- subset(data, lipid_class == lipid_class & ChainLength >= 3 & !is.na(ChainLength) & !is.na(`logFC (old_vs_young)`))
  
  # Calculate correlation between ChainLength and logFC (old_vs_young)
  correlation_value <- cor(correlation_data$ChainLength, correlation_data$`logFC (old_vs_young)`, method = "spearman")
  p_value <- cor.test(correlation_data$ChainLength, correlation_data$`logFC (old_vs_young)`, method = "spearman")$p.value
  
  # Add a red line if correlation is significant
  if (p_value < 0.05) {
    sp <- sp + geom_smooth(data = correlation_data, method = "lm", aes(group = NULL), color = "red", se = FALSE) +
      ggtitle(paste("Lipid Class:", lipid_class))
  }
  
  # Print correlation coefficient and p-value
  cat("Lipid Class:", lipid_class, "\n")
  cat("Correlation coefficient:", correlation_value, "\n")
  cat("p-value:", p_value, "\n\n")
}

# Display the plot
print(sp)
