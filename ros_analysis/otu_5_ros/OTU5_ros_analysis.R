#This will analyze the otu5 dataset

library(ggplot2)
library(tidyr)
library(Bolstad2)
library(ggpubr)

#early timepoint
late <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/124933.summary_info.csv",stringsAsFactors = F)
#later timepoint
early <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/184838.summary_info.csv" ,stringsAsFactors = F)

plot <- ggplot(late, aes(x=Time, y=measurement, color=Name))+
  geom_line(size=1.5)+
  #  geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd))+
  ylab("Photocount")+
  xlab("Seconds")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') 
plot

#late$Time <- late$Time + max(early$Time) + 120

combined <- rbind(early, late)
combined <- late
combined$measurement <- combined$measurement-min(combined$measurement)

combined$Name <- as.character(combined$Name)

combined <- combined[combined$Name != "1949+Pta22",]

ggplot(combined[combined$Name != "PvB+Pta22",], aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  #geom_vline(aes(xintercept=14))+
  #scale_y_log10()+
  scale_color_brewer(palette = "Dark2")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,4,10,20,30,40,50,60))


auc_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/124933.auc_data.csv")
#auc_data <- auc_data[auc_data$Name != "Pta_100nM",]
#auc_data$Name <- as.character(auc_data$Name)
#auc_data$Name[auc_data$Name == "1949_1uM_fls2"] <- "fls2_1949_1uM_Pta5nM"
#auc_data$Name <- factor(auc_data$Name, levels = c("Pta_5nM", '1949_5nM_Pta5nM', "1949_10nM_Pta5nM", "1949_100nM_Pta5nM", "1949_1uM_Pta5nM", "1949_1uM" ,"fls2_1949_1uM_Pta5nM" ))
auc_data$AUC <- log10(auc_data$AUC)

#get the z-score
auc_data$AUC <- ( auc_data$AUC - mean(auc_data$AUC[auc_data$Name == "mock+Pta22"]) ) / sd(auc_data$AUC[auc_data$Name == "mock+Pta22"])

auc_data <- auc_data[auc_data$Name !="1949+Pta22" & auc_data$Name !="PvB+Pta22",]

auc_sum <- auc_data %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise(total = n(),
            se_d = sd(AUC)/sqrt(total), 
            mean_d = mean(AUC))

set.seed(1994)
library(agricolae)

simple_anova <- aov(AUC ~ Name, auc_data)
summary(simple_anova)

anova_test <- HSD.test(simple_anova, trt = 'Name')
plot(anova_test)

box_auc <- ggplot(auc_data, aes(x=Name, y=AUC))+
  geom_errorbar(data=auc_sum, aes(x=Name, y=mean_d, ymin=mean_d-se_d, ymax=mean_d+se_d), width=.3, size=1.1)+
  geom_point(size=1.5)+
  geom_text(data = anova_test$groups, aes(x=row.names(anova_test$groups), label=groups, y=3), size=5)+
  labs(x="", y="Log10(RLU AUC) Pta Z-Score")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
box_auc


auc_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/190613.auc_data.csv")
auc_data$AUC <- log10(auc_data$AUC)

#get the z-score
auc_data$AUC <- ( auc_data$AUC - mean(auc_data$AUC[auc_data$Name == "mock+Pta22"]) ) / sd(auc_data$AUC[auc_data$Name == "mock+Pta22"])

auc_sum <- auc_data %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise(total = n(),
                   se_d = sd(AUC)/sqrt(total), 
                   mean_d = mean(AUC))

set.seed(1994)
library(agricolae)

simple_anova <- aov(AUC ~ Name, auc_data)
summary(simple_anova)

anova_test <- HSD.test(simple_anova, trt = 'Name')
plot(anova_test)

box_auc <- ggplot(auc_data, aes(x=Name, y=AUC))+
  geom_errorbar(data=auc_sum, aes(x=Name, y=mean_d, ymin=mean_d-se_d, ymax=mean_d+se_d), width=.3, size=1.1)+
  geom_point(size=1.5)+
  geom_text(data = anova_test$groups, aes(x=row.names(anova_test$groups), label=groups, y=8.5), size=5)+
  labs(x="", y="Log10(RLU AUC) Pta Z-Score")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
box_auc
