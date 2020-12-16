#this script will be to analyze the antagonist ROS data

library(ggplot2)
library(tidyr)
library(Bolstad2)

#early timepoint
early <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/161257.summary_info.csv",stringsAsFactors = F)
#later timepoint
late <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/134034.summary_info.csv" ,stringsAsFactors = F)

plot <- ggplot(early, aes(x=Time, y=measurement, color=Name))+
  geom_line(size=1.5)+
  #  geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd))+
  ylab("Photocount")+
  xlab("Seconds")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') 

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


#Trial 2
early <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/123102.summary_info.csv",stringsAsFactors = F)
#later timepoint
late <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/132351.summary_info.csv" ,stringsAsFactors = F)

plot <- ggplot(late, aes(x=Time, y=measurement, color=Name))+
  geom_line(size=1.5)+
  #  geom_errorbar(aes(ymin=measurement-sd, ymax=measurement+sd))+
  ylab("Photocount")+
  xlab("Seconds")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') 

late$Time <- late$Time + max(early$Time) + 120


combined <- rbind(early, late)

cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(ggsci)

#final figure 3k
combined$measurement <- combined$measurement - min(combined$measurement)
ggplot(combined, aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
  ylab("Photocount")+
  xlab("Minutes")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') +
  geom_vline(aes(xintercept=14))+
  #scale_y_log10()+
  #scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values=cbp2)+
  #scale_color_npg()+
  theme_classic(base_size = 14)+
  labs(y="ROS Burst (RLU)")+
  scale_x_continuous(breaks=c(0,6,13,16,26,36,46,56))

combined <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/132351.auc_data.csv" ,stringsAsFactors = F)
#should i get all three results?

set.seed(1994)
library(multcomp)
#first adjust values with a z-score
combined$AUC <- log10(combined$AUC)
combined$Name <- factor(combined$Name, levels = c("Pta", "Pta+Pa20", "Pta+2004", "Pta+1949", "Pta+2009", "Pta+5014", "Pta+1186", "Pta+1410"))

#z-score
combined$AUC <- ( combined$AUC - mean(combined$AUC[combined$Name == "Pta"]) )/sd(combined$AUC[combined$Name == "Pta"])

simple_anova <- aov(AUC ~ Name, combined)
summary(simple_anova)

t <- glht(simple_anova, linfct=mcp(Name="Dunnett"))
summary(t)

names_t <- summary(t)
signif_dunnet <- c()
signif_dunnet <- append(signif_dunnet, "")
for ( i in names_t$test$pvalues ) {
  if (i < 0.05) {
    signif_dunnet <- append(signif_dunnet, "*")
  }
  else {
    signif_dunnet <- append(signif_dunnet, "")
  }
}

dunn_signif <- data.frame(levels(combined$Name), signif_dunnet)


error_bar_data <- combined %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise( numb=n(),
                    se=sd(AUC)/sqrt(numb),
                    mean=mean(AUC))
#final figure 3j
box_auc <- ggplot(combined, aes(x=Name, y=AUC))+
  geom_errorbar(data=error_bar_data, aes(x=Name, y=mean, ymin=mean-se, ymax=mean+se), width=.3, size=1.1)+
  geom_point(size=1.5)+
  geom_text(data = dunn_signif, aes(x=levels.combined.Name., label=signif_dunnet, y=5), size=8)+
  labs(x="", y="Log10(RLU AUC) Pta Z-Score")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))

box_auc

#THIS IS LOOKING AT THE TIME THE MAX ROS BURST IS OBSERVED

combined <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/132351.txt" ,stringsAsFactors = F)
#combined <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/134034.txt" ,stringsAsFactors = F)


combined$Time <- combined$Time/60
combined  <-combined[combined$Time < 31,]

combined <- pivot_longer(combined,A1:A12, names_to = 'well', values_to = 'RLU')

max_data <- c()
for ( i in unique(combined$Name) ) {
  for (j in unique(combined$well) ) {
    if (j == "A12") {
      next
    }
    max_data <- append(max_data, c(i, j, combined$Time[combined$Name == i & combined$well == j & combined$RLU == max(combined$RLU[combined$Name == i & combined$well == j])]) )
  }
}

max_data <- matrix(max_data, ncol=3, byrow = T)
max_data <- as.data.frame(max_data)
max_data$V3 <- as.numeric(as.character(max_data$V3))

ggplot(max_data, aes(x=V1, y=V3))+
  geom_boxplot()+
  geom_point(position = position_jitter())


simple_anova <- aov(V3 ~ V1, max_data)
summary(simple_anova)

t <- glht(simple_anova, linfct=mcp(V1="Dunnett"))
summary(t)

t <- summary(t)
#t$test$pvalues <- p.adjust(t$test$pvalues, method = 'fdr')
signif_dunnet <- c()
signif_dunnet <- append(signif_dunnet, "")
for ( i in t$test$pvalues ) {
  if (i < 0.05) {
    signif_dunnet <- append(signif_dunnet, "*")
  }
  else {
    signif_dunnet <- append(signif_dunnet, "")
  }
}

dunn_signif <- data.frame(levels(max_data$V1), signif_dunnet)


max_data_summary <- max_data %>%
  dplyr::group_by(V1) %>%
  dplyr::summarise( numb=n(),
                    se=sd(V3)/sqrt(numb),
                    mean=mean(V3))

box_auc <- ggplot(max_data, aes(x=V1, y=V3))+
  geom_errorbar(data=max_data_summary, aes(x=V1, y=mean, ymin=mean-se, ymax=mean+se), width=.3, size=1.1)+
  geom_point(size=1.5, position = position_jitter(width = .1))+
  geom_text(data = dunn_signif, aes(x=levels.max_data.V1., label=signif_dunnet, y=30), size=8)+
  labs(x="", y="Time to maximum ROS burst")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
box_auc


library(ggpubr)
ggplot(max_data, aes(x=V1, y=V3))+
  geom_errorbar(data=max_data_summary, aes(x=V1, y=mean, ymin=mean-se, ymax=mean+se), width=.3, size=1.1)+
  geom_point(size=1.5, position = position_jitter(width = .1))+
  #geom_text(data = dunn_signif, aes(x=levels.max_data.V1., label=signif_dunnet, y=30), size=8)+
  stat_compare_means(ref.group = "Pta", method = "t.test")+
  labs(x="", y="Time to maximum ROS burst")+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))


library(agricolae)

hsd_res <- HSD.test(simple_anova,trt="V1")



