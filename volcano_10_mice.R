library(dplyr)
library(ggplot2)
library(ggrepel)
data<-Copy_of_BMP_sup_table_1
data1<-data[data$Tissues=="Brain",]

###add chain length
split_data <- strsplit(as.character(Copy_of_BMP_sup_table_1$Lipid_species), "[(:)]")

Copy_of_BMP_sup_table_1$chain_length <-sapply(split_data, function(x) x[2]) 
Copy_of_BMP_sup_table_1$double_bonds <-  sapply(split_data, function(x) x[3])


Copy_of_BMP_sup_table_1$double_bonds  <-  as.numeric(gsub("[^0-9]", "", Copy_of_BMP_sup_table_1$double_bonds))
summary(Copy_of_BMP_sup_table_1$double_bonds)

Copy_of_BMP_sup_table_1$chain_length<- as.numeric(gsub("[^0-9]","", Copy_of_BMP_sup_table_1$chain_length))   
summary(Copy_of_BMP_sup_table_1$chain_length)

unique(Copy_of_BMP_sup_table_1$Tissues)


#### OR from here import lipidatamouse set
Copy_of_BMP_sup_table_1$log2FC <- sign (Copy_of_BMP_sup_table_1$`logFC (old_vs_young)`) *log2(abs(Copy_of_BMP_sup_table_1$`logFC (old_vs_young)`) +1)
#first add lipid type column to dataset
# Remove all digits from the Lipid species column
Copy_of_BMP_sup_table_1$LipidType <- gsub("[0-9]", "", Copy_of_BMP_sup_table_1$Lipid_species)
unique(Copy_of_BMP_sup_table_1$LipidType)

# earase () for only specific columns because typ PLE(O) and PC(O) have it in name
# Define the pattern to remove
pattern_to_remove <- "\\(:\\)"

# Define categories to change
categories_to_change <- c("BMP(:)", "CL(:)", "LPG(:)", "LPI(:)", "LPS(:)", "MLCL(:)", "PA(:)", "PG(:)",
                          "PI(:)", "PS(:)", "LPC(:)", "LPE(:)", "PC(:)", "PE(:)", "SLBPA(:)")

# Apply changes only to specified categories
Copy_of_BMP_sup_table_1$LipidType <- ifelse(Copy_of_BMP_sup_table_1$LipidType %in% categories_to_change, 
                                            gsub(pattern_to_remove, "", Copy_of_BMP_sup_table_1$LipidType), 
                                            Copy_of_BMP_sup_table_1$LipidType)
unique(Copy_of_BMP_sup_table_1$LipidType)  #alleen : nog verwijderen bij sommige
Copy_of_BMP_sup_table_1$LipidType <- gsub("[:\\-]", "", Copy_of_BMP_sup_table_1$LipidType)



for (tissue in unique(Tissues)) {
  # Filter data for the current lipid class
  lipid_class_data <- data %>% filter(Tissues == tissue)
  
  lipid_class_data$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  lipid_class_data$diffexpressed[lipid_class_data$`logFC (old_vs_young)` > 1 & lipid_class_data$`adj.P.Val (old_vs_young)` < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  lipid_class_data$diffexpressed[lipid_class_data$`logFC (old_vs_young)` < -1 & lipid_class_data$`adj.P.Val (old_vs_young)`< 0.05] <- "DOWN"
  lipid_class_data$diffexpressed[lipid_class_data$chain_length %% 2 != 0] <- "ODD_CHAIN"
  
  # Highlight lipid chains longer than 50
  #lipid_class_data$diffexpressed[lipid_class_data$chain_length > 40] <- "LONG_CHAIN"
  
  # Convert directly in the aes()
  mycolors <- c("steelblue3", "gray", "hotpink2", "springgreen3")
  names(mycolors) <- c("DOWN", "NO", "UP", "odd chain length")
  
  lipid_class_data$delabel <- NA
  lipid_class_data$delabel[lipid_class_data$diffexpressed != "NO"] <- lipid_class_data$Lipid_species[lipid_class_data$diffexpressed != "NO"]
  
  lipid_plot <- ggplot(lipid_class_data, aes(x=`logFC (old_vs_young)`, y=-log10(`adj.P.Val (old_vs_young)`), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=mycolors) +
    geom_vline(xintercept=c(-1, 1), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  
  # Save the plot to a file
  ggsave(filename=paste0("lipid_plot_", tissue, ".png"), plot=lipid_plot, width=8, height=6)
}
