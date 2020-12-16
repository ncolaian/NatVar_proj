#this script will analyze the 10_nm data

library(ggplot2)
library(ggpubr)

pos_control <- "Pta"
neg_control <- "PtaDA"

path_ros <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/10nm_data_4_4_19/"

#get file names along with location


file_list <- list.files(path = path_ros, pattern = "*auc_data.csv")

#read the files
csv_list <- list()
for ( i in 1:length(file_list)) {
  csv_list[[i]] <- read.csv(paste(path_ros, file_list[i], sep = ""), stringsAsFactors = F)
  csv_list[[i]]$plate <- i
  if ( i == 1 ) {
    next
  }
  csv_list[[1]] <- rbind(csv_list[[1]], csv_list[[i]])
}

csv_list[[1]]$Name[csv_list[[1]]$Name == "pta"] <- "Pta"
csv_list[[1]]$Name[csv_list[[1]]$Name == "ptada"] <- "PtaDA"
csv_list[[1]]$Name[csv_list[[1]]$Name == "pa"] <- "Pa"
csv_list[[1]]$Name[csv_list[[1]]$Name == "pa20"] <- "Pa20"
csv_list[[1]] <- csv_list[[1]][csv_list[[1]]$Name != "10pa1001186",]
csv_list[[1]]$plate <- factor(as.character(csv_list[[1]]$plate))
csv_list[[1]]$Name <- factor(as.character(csv_list[[1]]$Name))
csv_list[[1]]$Name <- relevel(csv_list[[1]]$Name, ref = "PtaDA")

data_ros <- csv_list[[1]]
data_ros$AUC <- log10(data_ros$AUC)

data_ros[data_ros$AUC == 0,]
library(lme4)
library(lmerTest)
trial <- lm(AUC ~ Name + plate, data = data_ros )
trial2 <- lmer(AUC ~ Name + (1|plate), data = data_ros )
anova(trial2, trial)

s <- summary(trial2)
effect_table <- s$coefficients
effect_table <- effect_table[2:nrow(effect_table),]
effect_table <- as.data.frame(effect_table)
effect_table$names <- gsub("Name", "",row.names(effect_table))

ros_levels <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/882_data/updated_work_10_23_18/new_clusters/all_nat_variants_ptanorm_wheel_11_12_19.csv")
effect_table$names[!(effect_table$names %in% ros_levels$V2)]

effect_table <- effect_table[effect_table$names != "1787",]

effect_table$names <- factor(effect_table$names, levels=unique(ros_levels$V2[order(ros_levels$V3, decreasing = T)]))

effect_table$Estimate <- as.numeric( as.character(effect_table$Estimate) )
effect_table$`Std. Error` <- as.numeric( as.character(effect_table$`Std. Error`) )
effect_table$`Pr(>|t|)` <- as.numeric( as.character(effect_table$`Pr(>|t|)`) )

colnames(effect_table)<- c("Estimate", "se", 'df', "t-val", "pval", "names")
effect_table$pval <- p.adjust(effect_table$pval, method = "fdr")
effect_table$plot_p[effect_table$pval < 0.01] <- "*"

ggplot(effect_table[!(effect_table$names %in% seq(5000,5008)),], aes(names, Estimate))+
  geom_errorbar(aes(ymin=Estimate-se, ymax=Estimate+se))+
  geom_text(inherit.aes = T, aes( y=2, label=plot_p), size=8)+
  xlab("Peptide ID")+
  ylab("log10(AUC) Difference to PtaDA")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

###########
# END of this. Using the AUC did not work because of how low PtaDA is #
#################################
ptada_means_std <- c() 
for ( i in 1:nlevels(data_ros$plate) ) {
  data_ros$AUC[data_ros$plate == i] <- data_ros$AUC[data_ros$plate == i] - min(data_ros$AUC[data_ros$plate == i])+ 1
  ptada_means_std <- append(ptada_means_std, c(i, mean(data_ros$AUC[data_ros$plate == i & data_ros$Name == neg_control]), sd(data_ros$AUC[data_ros$plate == i & data_ros$Name == neg_control])))
}
ptada_means_std <- as.data.frame(matrix(ptada_means_std, ncol = 3, byrow = T))

#**Talk to Corbin about this**
z_scores_per_peptide <- c()
for ( i in ptada_means_std[,1] ) {
  i <- as.character(i)
  for ( j in unique(data_ros$Name[csv_list[[1]]$plate == i]) ) {
    for (f in data_ros$AUC[data_ros$Name == j & data_ros$plate == i ]) {
      z_scores_per_peptide <- append(z_scores_per_peptide, c(j, (f - ptada_means_std$V2[ptada_means_std$V1 == i]) ) )
    }
  }
}
z_scores_per_peptide <- as.data.frame(matrix(z_scores_per_peptide, ncol = 2, byrow = T), stringsAsFactors = F)
z_scores_per_peptide$V2 <- as.numeric(z_scores_per_peptide$V2)


#csv_list[[1]]$AUC <- log10(csv_list[[1]]$AUC)
ggplot(data_ros, aes(Name, AUC, color=plate))+
  geom_boxplot(outlier.colour  = NA)+
  xlab("Peptide ID")+
  ylab("log10(AUC)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##########################
##Get the ROS burst data##
##########################
ros_levels <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/882_data/updated_work_10_23_18/new_clusters/all_nat_variants_ptanorm_wheel_11_12_19.csv")
z_scores_per_peptide$V1[!(z_scores_per_peptide$V1 %in% ros_levels$V2)]

z_scores_per_peptide <- z_scores_per_peptide[z_scores_per_peptide$V1 != "1787",]
ros_levels$V2[ros_levels$V3 > 0 & !(ros_levels$V2 %in% z_scores_per_peptide$V1)]

z_scores_per_peptide$V1 <- factor(z_scores_per_peptide$V1, levels=unique(ros_levels$V2[order(ros_levels$V3, decreasing = T)]))
table(z_scores_per_peptide$V1)

#264 has a really bad outlier
z_scores_per_peptide <- z_scores_per_peptide[!(z_scores_per_peptide$V1 == 264 & z_scores_per_peptide$V2 < 0),]

#Dunnets test
library(multcomp)
z_scores_per_peptide$V1 <- relevel(z_scores_per_peptide$V1, ref = "PtaDA")
anova_a_scores <- aov(V2 ~ V1, z_scores_per_peptide)
dunnet_model <- glht(anova_a_scores, linfct=mcp(V1="Dunnett"))
s <- summary(dunnet_model)
p_vals <- s$test$pvalues
p_vals <- c(1, p_vals)
p_vals[p_vals < 0.05] <- "*"
p_vals[p_vals > 0.05] <- ""
dunnet_frame <- data.frame(levels(z_scores_per_peptide$V1)[levels(z_scores_per_peptide$V1) %in% z_scores_per_peptide$V1], p_vals)
colnames(dunnet_frame) <- c("pep", "pvalue")
z_scores_per_peptide$V1 <- factor(z_scores_per_peptide$V1, levels=unique(ros_levels$V2[order(ros_levels$V3, decreasing = T)]))


ggplot(z_scores_per_peptide, aes(V1, V2))+
  geom_boxplot(outlier.colour  = NA)+
  geom_text(inherit.aes = F,data = dunnet_frame, aes(x = pep, y=16, label=pvalue), size=8)+
  xlab("Peptide ID")+
  ylab("Log10(AUC) PtaDA Z−Score")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

z_scores_per_peptide$V1 <- factor(z_scores_per_peptide$V1, levels = levels(order_names))
ggplot(z_scores_per_peptide, aes(V1, V2))+
  geom_boxplot(outlier.shape = NA)
