#this script will be to analyze the antagonist ROS data

library(ggplot2)
library(tidyr)
library(patchwork)

data_1471 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak_ros_12_10_19_outliers/150756.summary_info.csv",stringsAsFactors = F)

data_1391 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak_ros_12_10_19_outliers/130527_final.summary_info.csv",stringsAsFactors = F)

data_1857 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak_ros_12_10_19_outliers/134034.summary_info.csv",stringsAsFactors = F)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#662D91", "#D55E00", "#CC79A7", "Black", "grey")

data_1471$Time <- data_1471$Time/60
data_1857$Time <- data_1857$Time/60
data_1391$Time <- data_1391$Time/60

data_1471 <- data_1471[data_1471$Time <= 30,]
data_1857 <- data_1857[data_1857$Time <= 30,]
data_1391 <- data_1391[data_1391$Time <= 30,]

ggplot(data_1471, aes(x=Time, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  theme_classic()+
  scale_color_manual(values = cbPalette)+
  scale_x_continuous(breaks=c(0,4,10,20,30))

ggplot(data_1391, aes(x=Time, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  scale_color_manual(values = cbPalette)+
  theme_classic()+
  scale_x_continuous(breaks=c(0,4,10,20,30))

ggplot(data_1857, aes(x=Time, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  scale_color_manual(values = cbPalette)+
  theme_classic()+
  scale_x_continuous(breaks=c(0,4,10,20,30))

data_1471 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak_ros_12_10_19_outliers/150756.auc_data.csv",stringsAsFactors = F)
data_1391 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak_ros_12_10_19_outliers/130527.auc_data.csv",stringsAsFactors = F)
data_1857 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak_ros_12_10_19_outliers/134034.auc_data.csv",stringsAsFactors = F)

data_1391$log_auc <- log10(data_1391$AUC)
data_1471$log_auc <- log10(data_1471$AUC)
data_1857$log_auc <- log10(data_1857$AUC)

auc_sum <- data_1391 %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise(total = n(),
                   se_d = sd(log_auc)/sqrt(total), 
                   mean_d = mean(log_auc))

set.seed(1994)
library(agricolae)

simple_anova <- aov(log_auc ~ Name, data_1391)
summary(simple_anova)

anova_test <- HSD.test(simple_anova, trt = 'Name')
plot(anova_test)

box_auc <- ggplot(data_1391, aes(x=Name, y=log_auc))+
  geom_errorbar(data=auc_sum, aes(x=Name, y=mean_d, ymin=mean_d-se_d, ymax=mean_d+se_d), width=.3, size=1.1)+
  geom_point(size=1.5)+
  geom_text(data = anova_test$groups, aes(x=row.names(anova_test$groups), label=groups, y=7.7), size=5)+
  labs(x="", y="Log10(RLU AUC)")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
box_auc

auc_sum <- data_1471 %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise(total = n(),
                   se_d = sd(log_auc)/sqrt(total), 
                   mean_d = mean(log_auc))

simple_anova <- aov(log_auc ~ Name, data_1471)
summary(simple_anova)

anova_test <- HSD.test(simple_anova, trt = 'Name')
plot(anova_test)

box_auc <- ggplot(data_1471, aes(x=Name, y=log_auc))+
  geom_errorbar(data=auc_sum, aes(x=Name, y=mean_d, ymin=mean_d-se_d, ymax=mean_d+se_d), width=.3, size=1.1)+
  geom_point(size=1.5)+
  geom_text(data = anova_test$groups, aes(x=row.names(anova_test$groups), label=groups, y=7.7), size=5)+
  labs(x="", y="Log10(RLU AUC)")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
box_auc

auc_sum <- data_1857 %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise(total = n(),
                   se_d = sd(log_auc)/sqrt(total), 
                   mean_d = mean(log_auc))

simple_anova <- aov(log_auc ~ Name, data_1857)
summary(simple_anova)

anova_test <- HSD.test(simple_anova, trt = 'Name')
plot(anova_test)

box_auc <- ggplot(data_1857, aes(x=Name, y=log_auc))+
  geom_errorbar(data=auc_sum, aes(x=Name, y=mean_d, ymin=mean_d-se_d, ymax=mean_d+se_d), width=.3, size=1.1)+
  geom_point(size=1.5)+
  geom_text(data = anova_test$groups, aes(x=row.names(anova_test$groups), label=groups, y=7.7), size=5)+
  labs(x="", y="Log10(RLU AUC)")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
box_auc

#breaks for y-axis
y_breaks <- c(0,5000,10000,15000,20000,25000)

#Alternative black plots
data_1471$plant[data_1471$Name %in% c("bak1-4-1471","bak1-4-Pa")] <- "bak1-4"
data_1471$plant[data_1471$Name %in% c("bak1-5-1471","bak1-5-Pa")] <- "bak1-5"
data_1471$plant[data_1471$Name %in% c("Col-1471","Col-Pa")] <- "Col-0"
data_1471$plant[data_1471$Name %in% c("fls2/efr-1471","fls2/efr-Pa")] <- "fls2/efr"

#Peptide
data_1471$peptide[data_1471$Name %in% c("bak1-4-1471","bak1-5-1471","Col-1471", "fls2/efr-1471" )] <- "1471"
data_1471$peptide[data_1471$Name %in% c("bak1-4-Pa","bak1-5-Pa","Col-Pa", "fls2/efr-Pa" )] <- "Pa"
data_1471$measurement <- data_1471$measurement - min(data_1471$measurement)

p1471 <- ggplot(data_1471[data_1471$Time/60 < 30,], aes(x=Time/60, y=measurement, color=plant, shape=peptide))+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se), size=.3, alpha=.75)+
  geom_line(size=.3)+
  geom_point(size=4, alpha = .75)+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Plant', shape="Peptide") +
  #scale_color_grey()+
  scale_color_manual(values = cbPalette)+
  theme_classic()+
  scale_x_continuous(breaks=c(0,4,10,12,20,30,40,50))+
  scale_y_continuous(breaks = y_breaks, limits = c(0,27500))
####
data_1391$plant[data_1391$Name %in% c("bak1-4-1391","bak1-4-Pa")] <- "bak1-4"
data_1391$plant[data_1391$Name %in% c("bak1-5-1391","bak1-5-Pa")] <- "bak1-5"
data_1391$plant[data_1391$Name %in% c("Col-1391","Col-Pa")] <- "Col-0"
data_1391$plant[data_1391$Name %in% c("fls2/efr-1391","fls2/efr-Pa")] <- "fls2/efr"

#Peptide
data_1391$peptide[data_1391$Name %in% c("bak1-4-1391","bak1-5-1391","Col-1391", "fls2/efr-1391" )] <- "1391"
data_1391$peptide[data_1391$Name %in% c("bak1-4-Pa","bak1-5-Pa","Col-Pa", "fls2/efr-Pa" )] <- "Pa"
data_1391$measurement <- data_1391$measurement - min(data_1391$measurement)

p1391 <- ggplot(data_1391[data_1391$Time/60 < 30,], aes(x=Time/60, y=measurement, color=plant, shape=peptide))+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se), size=.3, alpha=.75)+
  geom_line(size=.3)+
  geom_point(size=4, alpha = .75)+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Plant', shape="Peptide") +
  #scale_color_grey()+
  scale_color_manual(values = cbPalette)+
  theme_classic()+
  scale_x_continuous(breaks=c(0,4,10,12,20,30,40,50))+
  scale_y_continuous(breaks = y_breaks, limits = c(0,27500))

####
data_1857$plant[data_1857$Name %in% c("bak1-4-1857","bak1-4-Pa")] <- "bak1-4"
data_1857$plant[data_1857$Name %in% c("bak1-5-1857","bak1-5-Pa")] <- "bak1-5"
data_1857$plant[data_1857$Name %in% c("Col-1857","Col-Pa")] <- "Col-0"
data_1857$plant[data_1857$Name %in% c("fls2/efr-1857","fls2/efr-Pa")] <- "fls2/efr"

#Peptide
data_1857$peptide[data_1857$Name %in% c("bak1-4-1857","bak1-5-1857","Col-1857", "fls2/efr-1857" )] <- "1857"
data_1857$peptide[data_1857$Name %in% c("bak1-4-Pa","bak1-5-Pa","Col-Pa", "fls2/efr-Pa" )] <- "Pa"
data_1857$measurement <- data_1857$measurement - min(data_1857$measurement)

p1857 <- ggplot(data_1857[data_1857$Time/60 < 30,], aes(x=Time/60, y=measurement, color=plant, shape=peptide))+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se), size=.3, alpha=.75)+
  geom_line(size=.3)+
  geom_point(size=4, alpha = .75)+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Plant', shape= "Peptide") +
  #scale_color_grey()+
  scale_color_manual(values = cbPalette)+
  theme_classic()+
  scale_x_continuous(breaks=c(0,4,10,12,20,30,40,50))+
  scale_y_continuous(breaks = y_breaks, limits = c(0,27500))


(p1391 | p1857 | p1471)

