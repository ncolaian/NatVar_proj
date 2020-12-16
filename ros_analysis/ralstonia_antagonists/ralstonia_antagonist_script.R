#Ralstonia antagonist data

library(ggplot2)
library(tidyr)
library(Bolstad2)

#early timepoint
rep1 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/ralstonia_antagonists/135029.auc_data.csv",stringsAsFactors = F)
#rep2 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/ralstonia_antagonists/172403.auc_data.csv",stringsAsFactors = F)
#rep3 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/ralstonia_antagonists/201659.auc_data.csv",stringsAsFactors = F)
#rep4 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/ralstonia_antagonists/190325.auc_data.csv",stringsAsFactors = F)

rep1$rep <- "1"
rep2$rep <- "2"
#rep3$rep <- "3"
rep4$rep <- "4"

combined <- rbind(rep1, rep2, rep4)
combined <- rep1
combined$log10_val <- log10(combined$AUC)

#Z-score normalize
combined$log10_val <- ( combined$log10_val - mean(combined$log10_val[combined$Name == "mock"]) ) / sd(combined$log10_val[combined$Name == "mock"])

combined <- combined[combined$Name != "5007",]
combined$Name[combined$Name == "mock"] <- "Pta22"
combined$Name <- factor(combined$Name, levels = c("Pta22", "pa20", "5001", "5002", "5003", "5004", "5005", "5007"))

auc_sum <- combined %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise(total = n(),
                   se_d = sd(log10_val)/sqrt(total), 
                   mean_d = mean(log10_val))

set.seed(1994)
library(agricolae)

simple_anova <- aov(log10_val ~ Name, combined)
summary(simple_anova)

anova_test <- HSD.test(simple_anova, trt = 'Name')
plot(anova_test)

box_auc <- ggplot(combined, aes(x=Name, y=log10_val))+
  geom_errorbar(data=auc_sum, aes(x=Name, y=mean_d, ymin=mean_d-se_d, ymax=mean_d+se_d), width=.3, size=1.1)+
  geom_point(size=1.5)+
  geom_text(data = anova_test$groups, aes(x=row.names(anova_test$groups), label=groups, y=3), size=5)+
  labs(x="", y="Log10(RLU AUC) Pta Z-Score")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
box_auc


#Plot the ROS burst curve
late <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/ralstonia_antagonists/135029.summary_info.csv",stringsAsFactors = F)

late$Name <- as.character(late$Name)
late <- late[late$Name != "5007",]
late$Name[late$Name == "mock"] <- "Pta22"
late$Name <- factor(late$Name, levels = c("Pta22", "pa20", "5001", "5002", "5003", "5004", "5005"))

late$measurement <- late$measurement - min(late$measurement)

ggplot(late, aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide', y="ROS Burst (RLU)", x= "Time (Minutes)") +
  #geom_vline(aes(xintercept=14))+
  #scale_y_log10()+
  scale_color_brewer(palette = "Dark2")+
  theme_classic()+
  scale_x_continuous(breaks=c(0,4,10,20,30,40,50,60))
