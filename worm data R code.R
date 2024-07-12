#other dataset with all wormstrains and rmean lifespan
splitstrain <- strsplit(as.character(mean_lifespan_all$`;rmean;Strain;;;;`), ",")
mean_lifespan_all$strain<-sapply(splitstrain, function(x) x[1]) 
mean_lifespan_all$rmean <-sapply(splitstrain, function(x) x[2]) 
mean_lifespan_all <- mean_lifespan_all[-1]

splitstrain2 <- strsplit(as.character(mean_lifespan_all$strain),"=")
mean_lifespan_all$strain <- sapply(splitstrain2, function(x) x[2])
splittreatment <- strsplit(as.character(mean_lifespan_all$rmean),";")
mean_lifespan_all$rmean <- sapply(splittreatment, function(x) x[2])

##
rdered_column <- Pheno_Data_2024_04_23_11_50_10[order(Pheno_Data_2024_04_23_11_50_10$rmean), ] #QX589 is reference worm
which(lipid_metadata_CTRL$Strain_code == "QX589") #WR000131  is QX589
#calculate for each lipid WR000131 mean value, by comparing to mean worm for that lipid
##FOR Alkanyl-TG, calculate logFC
lipodomicsFC <- lipid_data_CTRL[c(6,8,9)]
FCdata <- lipid_data_CTRL[11:97]

conditionA <- FCdata[, 2]  # Condition A is column 2 - control worm
conditionB <- FCdata[, c(1:87)]  

log2FC_df <- log2(conditionB / conditionA)
FCdata<- cbind(lipodomicsFC, log2FC_df)

library(ggplot2)
library(patchwork)
#normality check
splitdata <- split(lipid_data_CTRL, lipid_data_CTRL$Lipid.Class)
Alkanyl.TG.norcheck <- splitdata[["Alkanyl-TG"]]
scatter.smooth(Alkanyl.TG.norcheck$WR000137)  #data is normal

wormsplit <- split(FCdata, FCdata$Lipid.Class)
Alkanyl.TG <- wormsplit[["Alkanyl-TG"]]

#calculate r value
r_values <- numeric(length = ncol(Alkanyl.TG) - 3)  # 87 columns from 4 to 90

# Loop through columns 4 to 90 and calculate Spearman correlation
for (i in 4:ncol(Alkanyl.TG)) {
  correlation_result <- cor.test(Alkanyl.TG$cLength, Alkanyl.TG[, i], method = "pearson")
  r_values[i - 3] <- correlation_result$estimate  # Store the r value
}

r_Alkanyl.TG <- data.frame(Barcode = colnames(FCdata[c(4:90)]), R = r_values)

#plot r against lifespan
#get strain and average lifespan from Pheno_data...
barcodestrain <- lipid_metadata_CTRL[c(9,10)]
colnames(barcodestrain)[2] <- "strain"
lifespanbarcode <- merge(barcodestrain, mean_lifespan_all, by = "strain") 
#add r to corresponding worms
lifespanbarcode <- merge(lifespanbarcode, r_Alkanyl.TG, by = "Barcode")

lifespanbarcode$rmean <- as.numeric(as.character(lifespanbarcode$rmean))
summary(lifespanbarcode$rmean)

correlation_test <- cor.test(lifespanbarcode$rmean, lifespanbarcode$R, method = "spearman")
estimate <- correlation_test$estimate
p_value <- correlation_test$p.value
#plot plot
ggplot(lifespanbarcode, aes(x=rmean, y=R)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  scale_x_continuous(breaks = seq(floor(12.20), ceiling(21.97), by = 0.5))  +
  labs(title = "Alkanyl.TG", x = "rmean", y = "R") +
  annotate("text", x = Inf, y = Inf, label = paste("Correlation:", format(estimate, digits = 4), "P-value:", format(p_value, digits = 4)), 
           hjust = 1.1, vjust = 2, size = 4, color = "black")


