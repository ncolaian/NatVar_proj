#this script will be to analyze the primer data

library(ggplot2)
library(tidyr)
library(Bolstad2)
library(ggpubr)

#early timepoint
late <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/150539.summary_info.csv",stringsAsFactors = F)
#later timepoint
early <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/141256.summary_info.csv" ,stringsAsFactors = F)

plot <- ggplot(late, aes(x=Time, y=measurement, color=Name))+
  geom_line(size=1.5)+
  #  geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd))+
  ylab("Photocount")+
  xlab("Seconds")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') 

late$Time <- late$Time + max(early$Time) + 120

combined <- rbind(early, late)

pallete <- scale_color_brewer()
combined <- combined[combined$Name != "Pta_100nM",]
combined$Name[combined$Name == "1949_1uM_fls2"] <- "fls2_1949_1uM_Pta5nM"
combined$Name <- factor(combined$Name, levels = c("1949_1uM", "fls2_1949_1uM_Pta5nM", "Pta_5nM", "1949_5nM_Pta5nM", "1949_10nM_Pta5nM", "1949_100nM_Pta5nM", "1949_1uM_Pta5nM", "1949_5uM_Pta5nM"))

#This is figure S3c
ggplot(combined, aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  geom_vline(aes(xintercept=14))+
  #scale_y_log10()+
  scale_color_brewer(palette = "Dark2")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))

#Get the quantifications for the data
auc_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/150539.auc_data.csv")
auc_data <- auc_data[auc_data$Name != "Pta_100nM",]
auc_data$Name <- as.character(auc_data$Name)
auc_data$Name[auc_data$Name == "1949_1uM_fls2"] <- "fls2_1949_1uM_Pta5nM"
auc_data$Name <- factor(auc_data$Name, levels = c("Pta_5nM", '1949_5nM_Pta5nM', "1949_10nM_Pta5nM", "1949_100nM_Pta5nM", "1949_1uM_Pta5nM", "1949_1uM" ,"fls2_1949_1uM_Pta5nM" ))
auc_data$AUC <- log10(auc_data$AUC)
library(dplyr)

#get the z-score
auc_data$AUC <- ( auc_data$AUC - mean(auc_data$AUC[auc_data$Name == "Pta_5nM"]) ) / sd(auc_data$AUC[auc_data$Name == "Pta_5nM"])

auc_sum <- auc_data %>%
  group_by(Name) %>%
  summarise(total = n(),
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

###dat2
late <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/144240.summary_info.csv",stringsAsFactors = F)
#later timepoint
early <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/134942.summary_info.csv" ,stringsAsFactors = F)

late$Time <- late$Time + max(early$Time) + 120

combined <- rbind(early, late)

pallete <- scale_color_brewer()

ggplot(combined, aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  #geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  geom_vline(aes(xintercept=14))+
  #scale_y_log10()+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))

###dat3
# This is figure S3d
late <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/201502.summary_info.csv",stringsAsFactors = F)
#later timepoint
early <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/191923.summary_info.csv" ,stringsAsFactors = F)

late$Time <- late$Time + max(early$Time) + 120

combined <- rbind(early, late)

pallete <- scale_color_brewer()

ggplot(combined, aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  geom_vline(aes(xintercept=14))+
  #scale_y_log10()+
  scale_color_brewer(palette = "Dark2")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))

#Get the quantifications for the data
auc_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/201502.auc_data.csv")
auc_data$AUC <- log10(auc_data$AUC)
library(dplyr)

#get the z-score
auc_data$AUC <- ( auc_data$AUC - mean(auc_data$AUC[auc_data$Name == "Pta_5nM"]) ) / sd(auc_data$AUC[auc_data$Name == "Pta_5nM"])

auc_sum <- auc_data %>%
  group_by(Name) %>%
  summarise(total = n(),
            se_d = sd(AUC)/sqrt(total), 
            mean_d = mean(AUC))

simple_anova <- aov(AUC ~ Name, auc_data)
summary(simple_anova)

anova_test <- HSD.test(simple_anova, trt = 'Name')
plot(anova_test)

box_auc <- ggplot(auc_data, aes(x=Name, y=AUC))+
  geom_errorbar(data=auc_sum, aes(x=Name, y=mean_d, ymin=mean_d-se_d, ymax=mean_d+se_d), width=.3, size=1.1)+
  geom_point(size=1.5)+
  #geom_text(data = anova_test$groups, aes(x=row.names(anova_test$groups), label=groups, y=8.5), size=5)+
  stat_compare_means(comparisons = list(c("Pta_100nM", 'Pta_100nM+1949_1uM'),
                                        c('Pta_10nM', "Pta_10nM+1949_1uM"),
                                        c("Pta_5nM", "Pta_5nM+1949_1uM")),
                     method = "t.test", label = "p.format"
                     )+
  labs(x="", y="Log10(RLU AUC) Pta Z-Score")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))

box_auc
