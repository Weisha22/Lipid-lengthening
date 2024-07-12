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

## Georges plot
names(MOUSEdata_clean)[names(MOUSEdata_clean) == "logFC (old_vs_young)"] <- "logFC"
names(MOUSEdata_clean)[names(MOUSEdata_clean) == "P.Value (old_vs_young)"] <- "P.Value"

lipidlist<- split(MOUSEdata_clean, MOUSEdata_clean$LipidType)
plotslist <- list()


for (name in names(lipidlist)) {
  lip_data <- lipidlist[[name]]
  
  for (cat in unique(lip_data$Tissues)) {
    cat_data <- lip_data[lip_data$Tissues == cat, ]
    
    pl <- ggplot(cat_data, aes(chain_length ,double_bonds)) +
      geom_point( aes(color = logFC, 
                      size = -log10(P.Value))) +
      ggtitle(cat) +
      scale_color_gradientn(values = seq(0,1,0.1), colors= c("blue", "white", "red"))
    plotslist[[cat]] <- pl 
  }
  combined_plot <- wrap_plots(plotslist) 
  combined_plot <- combined_plot + plot_annotation(title = name)
  print(combined_plot)
  ggsave(filename = paste0(name, ".png"), plot = combined_plot,width = 20, height = 15)
}

