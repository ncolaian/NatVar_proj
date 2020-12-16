#This script will analyze the 1949 priming data with Pta and elf18

#I have the ROS data before the addition of flg22, but there is not interesting outcomes.

library(ggplot2)
library(tidyr)
library(ggpubr)
library(agricolae)
library(dplyr)

exp3 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/elf18_prim/200523.auc_data.csv")
exp2 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/elf18_prim/182936.auc_data.csv")
exp1 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/elf18_prim/171457.auc_data.csv")

exp3$Experiment <- "3"
exp2$Experiment <- "2"
exp1$Experiment <- "1"

combined <- rbind(exp3, exp2, exp1)

combined$logauc <- log10(combined$AUC)
combined$Experiment <- factor(combined$Experiment)

ggplot(combined, aes(Name, AUC, color=Experiment))+
  geom_boxplot()

exclude <- c("DW+1uM_elf18")
#Try an anova first
combined$Experiment <- factor(combined$Experiment)

anova_simple <- aov(logauc ~ Name + Experiment, data = combined[!combined$Name %in% exclude,])
interact_anova <- aov(logauc ~ Name * Experiment, data = combined[!combined$Name %in% exclude,])
anova(anova_simple, interact_anova)
#simple is better


tukey_test <- HSD.test(anova_simple, trt = "Name")
test <- tukey_test$groups
test$names <- row.names(test)

plot(tukey_test)




##########
exp3 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/elf18_prim/200523.summary_info.csv")
exp2 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/elf18_prim/182936.summary_info.csv")
exp1 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/elf18_prim/171457.summary_info.csv")

exp1$measurement <- exp1$measurement - min(exp1$measurement)
exp2$measurement <- exp2$measurement - min(exp2$measurement)
exp3$measurement <- exp3$measurement - min(exp3$measurement)

ggplot(exp3, aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  #scale_y_log10()+
  scale_color_brewer(palette = "Dark2")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))+
  scale_y_continuous(breaks=c(0,5000,10000,15000,20000,25000), limits = c(0,26000))
  theme(legend.position = "none")
  

ggplot(exp2, aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  #scale_y_log10()+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))+
  scale_y_continuous(breaks=c(0,5000,10000,15000,20000,25000), limits = c(0,26000))
  #theme(legend.position = "none")

#Supplemental figure 3E
ggplot(exp1[exp1$Name != "DW+1uM_Pta",], aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20),
        legend.position = "none")+
  labs(color='Peptide') +
  #scale_y_log10()+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))
  #scale_y_continuous(breaks=c(0,5000,10000,15000,20000,25000), limits = c(0,26000))

####
exp3_auc <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/elf18_prim/200523.auc_data.csv")
exp3_auc$AUC <- log10(exp3_auc$AUC)
exp3_auc <- exp3_auc[exp3_auc$Name != "DW+1uM_Pta",]
exp3_auc$Name <- factor(exp3_auc$Name, levels = c("DW+5nM_Pta", "1uM_1949+5nM_Pta",
                                                  "DW+10nM_elf18", "1uM_1949+10nM_elf18",
                                                  "DW+1uM_elf18"))

#z-score to 5uM Pta
exp3_auc$AUC <- (exp3_auc$AUC - mean(exp3_auc$AUC[exp3_auc$Name == 'DW+5nM_Pta']) )/sd(exp3_auc$AUC[exp3_auc$Name == 'DW+5nM_Pta'])

auc_sum <- exp3_auc %>%
  group_by(Name) %>%
  summarise(total = n(),
            se_d = sd(AUC)/sqrt(total), 
            mean_d = mean(AUC))

box_auc <- ggplot(exp3_auc, aes(x=Name, y=AUC))+
  geom_errorbar(data=auc_sum, aes(x=Name, y=mean_d, ymin=mean_d-se_d, ymax=mean_d+se_d), width=.3, size=1.1)+
  geom_point(size=1.5)+
  labs(x="", y="Log10(RLU AUC) Pta Z-Score")+
  theme_bw(base_size = 16)+
  stat_compare_means(data = exp3_auc, aes(groups=Name), comparisons = list(c('1uM_1949+10nM_elf18', 'DW+10nM_elf18'), c("1uM_1949+10nM_elf18", "DW+1uM_elf18"), c("DW+5nM_Pta", "1uM_1949+5nM_Pta")), 
                     method = "t.test", label = "p.format", label.y = c(10,20,11))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))

box_auc



#I wanna plot the same thing but without the 1uM data

exclude <- c("DW+1uM_elf18", "DW+1uM_Pta")

ggplot(exp3[!(exp3$Name %in% exclude),], aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  #scale_y_log10()+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))+
  scale_y_continuous(breaks=c(0,5000,10000), limits = c(0,10000))+
  theme(legend.position = "none")
  
  
  ggplot(exp2[!(exp2$Name %in% exclude),], aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  #scale_y_log10()+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))+
  scale_y_continuous(breaks=c(0,5000,10000), limits = c(0,10000))+
  theme(legend.position = "none")
  
  ggplot(exp1[!(exp1$Name %in% exclude),], aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20),
        legend.position = "none")+
  labs(color='Peptide') +
  #scale_y_log10()+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,6,13,16,26,30,40,50,60))+
  scale_y_continuous(breaks=c(0,5000,10000), limits = c(0,10000))
