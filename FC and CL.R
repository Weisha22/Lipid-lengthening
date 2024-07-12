library(ggplot2)
library(tidyverse)
#install.packages("readxl")
library(readxl)
library(ggpubr)
library(smplot2)
setwd("/Users/liweisha/Documents/lipid lengthening with age/human/rental")  
# Sample data (replace this with your own dataset)
# Example dataset structure:
# lipid_class: Categorized lipid class
# chain_length: Numeric chain length values
# fold_change: Numeric fold change values
data<-read.csv("Supplemental_Data_Table_5_MenAgingLipidomics.csv")

lipids <- data$Lipid_Species

# Extract chain length and double bonds using regular expressions
matches <- regmatches(lipids, regexec("(\\d+):(\\d+)", lipids))

# Create a data frame
lipid_data <- data.frame(
  Lipid_Species = lipids,
  ChainLength = as.integer(sapply(matches, function(x) x[2])),
  DoubleBonds = as.integer(sapply(matches, function(x) x[3])),
  lipidclass =  gsub("\\(.*\\)", "", lipids)
)

# Print the resulting data frame
print(lipid_data)



# add the chain length and double bonds to the lipid profile data
data<- full_join(data, lipid_data, by = "Lipid_Species")
write.csv(data, "men lipids with DB and CL.csv")

data<-read.csv("wen lipids with DB and CL.csv")
data<-data[c(-1086,-886,-1085,-887),]
data$logFC..Old_vs_Young. <- as.numeric(data$logFC..Old_vs_Young., na.rm = TRUE)

# Drop any rows with NA values that resulted from the coercion


# Create a scatter plot


sp <- ggplot(data, aes(ChainLength, data$logFC..Old_vs_Young.)) +
  geom_point()+
  sm_statCorr(corr_method = 'spearman')+
  facet_wrap(~lipidclass,scales = "free")

sp

pdf(file="cor_men_muscle_CL_FC.pdf",width=15,height=20)
print(sp)
dev.off()


##kidney
kidney<-Combilist
old_columns <- grep("Old$", colnames(kidney), value = TRUE)
young_columns <- grep("Young$", colnames(kidney), value = TRUE)
mid_columns<- grep("MiddleAged$", colnames(kidney), value = TRUE)
# Calculate the log2 fold change for each pair of columns
oldmean<-rowMeans(kidney[, old_columns])
youngmean<-rowMeans(kidney[, young_columns])
kidney$old_vs_young_log2_fold_change <- log2(oldmean / youngmean)

lipids <- kidney$IsolabMidas

# Extract chain length and double bonds using regular expressions
matches <- regmatches(lipids, regexec("(\\d+):(\\d+)", lipids))

# Create a data frame
lipid_data <- data.frame(
  IsolabMidas = lipids,
  ChainLength = as.integer(sapply(matches, function(x) x[2])),
  DoubleBonds = as.integer(sapply(matches, function(x) x[3]))
)
kidney<- full_join(kidney, lipid_data, by = "IsolabMidas")
write.csv(kidney, "kidney lipids with DB and CL.csv")
sp <- ggplot(kidney, aes(ChainLength, old_vs_young_log2_fold_change)) +
  geom_point()+
  sm_statCorr(corr_method = 'pearson')+
  facet_wrap(~idclasses,scales = "free")

sp

pdf(file="cor_kidney_CL_FC.pdf",width=15,height=20)
print(sp)
dev.off()


data<-read.csv("kidney lipids with DB and CL.csv")
correlation_results <- data.frame(lipidclass = character(), p_value = numeric(),log2_fold_change = numeric(),correlation_r_value = numeric())
lipidclass<-unique(data$idclasses)
lipidclass<-lipidclass[-22]
for (idclass in lipidclass){
  groupsubset<-subset(data,data$idclasses == idclass )
  corr_kidney<-cor.test(groupsubset$ChainLength,as.numeric(groupsubset$old_vs_young_log2_fold_change),method = "spearman")
  R_value <- corr_kidney$estimate
  p_value <- corr_kidney$p.value
  correlation_results <- bind_rows(correlation_results, data.frame(lipidclass = idclass, p_value= p_value, correlation_r_value = R_value))
  
}
correlation_results$group<-"kidney"
write.csv(correlation_results,file = "correlation results kidney.csv")
aa1<-read.csv("correlation results women_muscle.csv")
aa2<-read.csv("correlation results men_muscle.csv")
correlation_results<-rbind(aa1,aa2)

alpha <- 0.05
significance_levels <- ifelse(correlation_results$p_value >= alpha, "",
                              ifelse(correlation_results$p_value < 0.001, "***",
                                     ifelse(correlation_results$p_value < 0.01, "**", "*")))
correlation_results$significance<-Significance
p=ggplot(correlation_results, aes(x = group, y = reorder(lipidclass,correlation_r_value), fill = correlation_r_value)) +
  geom_tile()+
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red") +
  geom_text(aes(label = Significance), color = "black", size = 3) +
  coord_fixed()   
pdf(file="cor_kidney_heatmap.pdf",width=10,height=20)
print(p)
dev.off()

adj.P_value <-p.adjust(correlation_result$p.value, method = "bonferroni")
Hoch.P_value <- p.adjust(correlation_result$p.value, method = "hochberg")




Significance <- ifelse(adj.P_value < 0.001, "***",
                       ifelse(p_valuesframe$adj.P_value < 0.01, "**",
                              ifelse(p_valuesframe$adj.P_value < 0.05, "*", " ")))


Significance.Hoch <- ifelse(Hoch.P_value < 0.001, "***",
                                          ifelse(Hoch.P_value < 0.01, "**",
                                                 ifelse(Hoch.P_value < 0.05, "*", " ")))
