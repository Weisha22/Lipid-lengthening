
library(tidyverse)
library(dplyr)

#add double bonds and carbon chain length as a saparate column
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



###### check if i have missing values 
sum(complete.cases(Copy_of_BMP_sup_table_1$chain_length, Copy_of_BMP_sup_table_1$log2FC))  ##i have missing values
summary(Copy_of_BMP_sup_table_1$log2FC) 
which(is.na(Copy_of_BMP_sup_table_1$log2FC))

##remove these values
MOUSEdata_clean <- Copy_of_BMP_sup_table_1[-c(183, 1832, 2389, 2825, 2894),] #data is not normally distributed 


# create forloop based on category in column
library(ggplot2)
library(patchwork)


#ALL TOGETHER IN ONE PLOT
#create a list to combine all the plots created in forloop
lipidslist_clean <- split(MOUSEdata_clean, MOUSEdata_clean$LipidType)

output_dir <- "plots_doublebonds"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
plots_list <- list()
for (name in names(lipidslist_clean)) {
  lip_data <- lipidslist_clean[[name]]

    for (cat in unique(lip_data$Tissues)) {
  cat_data <- lip_data[lip_data$Tissues == cat, ]#define cat and create cat_data for each category
  
  correlation_result <- cor.test(cat_data$chain_length, cat_data$log2FC, method = "spearman")
  
  p <- ggplot(cat_data, aes(x = double_bonds, y = log2FC)) + 
    geom_point() +
    ggtitle(cat) +
    geom_smooth(method = "lm", se = FALSE) 
  plots_list[[cat]] <- p + annotate("text", x = Inf, y = -Inf, 
                                  label = paste("p-value:", format(correlation_result$p.value, digits = 4)),
                                  hjust = 1, vjust = -1)

  
    }

combined_plot <- wrap_plots(plots_list) 
combined_plot <- combined_plot + plot_annotation(title = name)
print(combined_plot)
ggsave(filename = file.path(output_dir, paste0(name, ".png")), plot = combined_plot, width = 15, height = 10)
}



## calculate the significante
p_valuesframe <- data.frame(Lipid = character() ,Tissues = character(), Correlation = numeric(),
                            P_value = numeric(), stringsAsFactors = FALSE)


for (name in names(lipidslist_clean)) {
  lip_data <- lipidslist_clean[[name]]
  
  
  for (cat in unique(lip_data$Tissues)) {
    cat_data <- lip_data[lip_data$Tissues == cat, ]#define cat and create cat_data for each category
    
    correlation_result <- cor.test(cat_data$double_bonds, cat_data$log2FC, method = "spearman")
    p_valuesframe <- rbind(p_valuesframe, data.frame(Lipid = name, Tissues = cat,
                                               Correlation = correlation_result$estimate,
                                               P_value = correlation_result$p.value))
  }
 
}
View(p_valuesframe)
which(is.na(p_valuesframe$Correlation))

summary(p_valuesframe$Correlation)

#Correct for multiple testing
#USE multiple testing correction on my calculated P-values
p_valuesframe$adj.P_value <-p.adjust(p_valuesframe$P_value, method = "bonferroni")
p_valuesframe$Hoch.P_value <- p.adjust(p_valuesframe$P_value, method = "hochberg")

p_valuesframe <- p_valuesframe[-175,] #is N.A. value
##add stars for significance


p_valuesframe$Significance <- ifelse(p_valuesframe$adj.P_value < 0.001, "***",
                                      ifelse(p_valuesframe$adj.P_value < 0.01, "**",
                                             ifelse(p_valuesframe$adj.P_value < 0.05, "*", " ")))


p_valuesframe$Significance.Hoch <- ifelse(p_valuesframe$Hoch.P_value < 0.001, "***",
                                     ifelse(p_valuesframe$Hoch.P_value < 0.01, "**",
                                            ifelse(p_valuesframe$Hoch.P_value < 0.05, "*", " ")))




#heatmap

ggplot(p_valuesframe, aes(x = Tissues, y = Lipid, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = Significance.Hoch), show.legend = FALSE, color = "black") + # Removing legend for significance
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Age related lipid lengthening in mice_", x = "Tissues", y = "Lipids")
  

#FOR not corrected
p_valuesframe$SignificanceNOcor <- ifelse(p_valuesframe$P_value < 0.001, "***",
                                                  ifelse(p_valuesframe$P_value < 0.01, "**",
                                                         ifelse(p_valuesframe$P_value < 0.05, "*", " ")))
p_valuesframe$starsNOcor <- ifelse(p_valuesframe$P_value < 0.001, "p < 0.001",
                                           ifelse(p_valuesframe$P_value < 0.01, "0.001 ≤ p < 0.01",
                                                  ifelse(p_valuesframe$P_value < 0.05, "0.01 ≤ p < 0.05", "p ≥ 0.05 ")))



ggplot(p_valuesframe, aes(x = Tissues, y = Lipid, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = SignificanceNOcor), show.legend = FALSE, color = "black") + # Removing legend for significance
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Age related lipid lengthening in mice NO correction", x = "Tissues", y = "Lipid classes")




 