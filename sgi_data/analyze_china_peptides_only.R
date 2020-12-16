#Sgi analysis code for the peptides synthesized in china
library(tidyr)
library(ggplot2)
library(dplyr)

full_sgi <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/sgi_data/china_peptides_sgi.csv")
full_sgi <- gather(full_sgi, "plate_col", "mg", A:H, na.rm = T)
full_sgi$mg <- full_sgi$mg/10
full_sgi$Plate <- as.factor(full_sgi$Plate)

#have to fix all the different spellings
full_sgi$Peptide[full_sgi$Peptide == "Mock"] <- "mock"
full_sgi$Peptide[full_sgi$Peptide == "Mock "] <- "mock"
full_sgi$Peptide[full_sgi$Peptide == "mock "] <- "mock"
full_sgi$Peptide[full_sgi$Peptide == "flg22 "] <- "Flg22"
full_sgi$Peptide[full_sgi$Peptide == "flg22"] <- "Flg22"
full_sgi$Peptide[full_sgi$Peptide == "pa"] <- "Pa"
full_sgi$Peptide[full_sgi$Peptide == "pa20"] <- "Pa20"
full_sgi$Peptide[full_sgi$Peptide == "pta"] <- "Pta"
full_sgi$Peptide[full_sgi$Peptide == "PtaDa"] <- "PtaDA"

#change 2609 to 2005
full_sgi$Peptide <- as.character(full_sgi$Peptide)
full_sgi$Peptide[full_sgi$Peptide == "2609"] <- "2005"
full_sgi$Peptide <- factor(full_sgi$Peptide)

#relevel for modelling
full_sgi$Peptide <- relevel(full_sgi$Peptide, ref = "PtaDA")
full_sgi$Plant <- relevel(full_sgi$Plant, ref="Fls2efr")
full_sgi$Experiment <- as.factor(full_sgi$Experiment)
full_sgi <- full_sgi[full_sgi$Experiment != 4,]

ggplot(na.omit(full_sgi[full_sgi$Peptide == "2609",]), aes(Plate, mg, color=Experiment))+
  geom_boxplot(outlier.shape = NA)

library(lme4)
library(lmerTest)
less_rand <- lmer(mg ~ Peptide + Experiment + (1|Plate), data = full_sgi)
summary(less_rand)

ggplot(full_sgi, aes(Peptide, mg)) +
  geom_boxplot(outlier.shape = NA)

t <- summary(less_rand)
s <- as.data.frame( t$coefficients )
colnames(s)[2] <- "StdE"
s_inter <- s[1,]
s <- s[2:(nrow(s)-2),]
row.names(s) <- gsub("Peptide","", row.names(s))
#s$Estimate <- s$Estimate + s_inter$Estimate[1]

#get the sgi data into a a format to add to itol
get_seq <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/youssef_data/old_arrays/finalized_useq_8_22_2018.csv", stringsAsFactors = F)
#The chinese peptides had numbering that I placed uniquely. I would like to get the sequenced based on this file
get_seq_2 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Kate_Vienna_data/peptides_to_order_china.csv")

wheel_data <- c()
for ( i in row.names(s) ) {
  if ( as.character(i) %in% as.character( get_seq$UseqId ) ) {
    stri <- get_seq$Peptides[get_seq$UseqId == i]
  }
  else if ( as.character(i) %in% as.character( get_seq_2$ID ) ) {
    stri <- as.character(get_seq_2$Peptides.to.order[get_seq_2$ID == i])
  }
  else {
    stri <- "A"
  }
  print(stri)
  wheel_data <- append(wheel_data, c(stri, i, as.numeric(s$Estimate[row.names(s) == i]), as.numeric(s$`Pr(>|t|)`[row.names(s) == i]) ))
}
wheel_data <- matrix(wheel_data, ncol = 4, byrow = T)
write.csv(wheel_data,file = "/Users/nicholascolaianni/Documents/dangl_lab/882_data/updated_work_10_23_18/new_clusters/china_variants_sgieffect_wheel_5_23_19.csv" , quote = F, row.names = F)





