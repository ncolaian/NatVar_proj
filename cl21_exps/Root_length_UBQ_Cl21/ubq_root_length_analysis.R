root_length <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/CL21_colonization_experiments/root_length/cl21_rep1_root_length_ubq.csv")
root_length <- rbind(root_length, read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/CL21_colonization_experiments/root_length/cl21_rep2_root_length_ubq.csv"))
root_length <- rbind(root_length, read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/CL21_colonization_experiments/root_length/cl21_rep3_root_length_ubq.csv"))

#we will not be presenting the fliC mutant data
root_length <- root_length[root_length$Strain != "deltafliC",]

root_length$Strain <- as.character(root_length$Strain)
root_length$Strain[root_length$Strain == "p176"] <- "5005"
root_length$Strain[root_length$Strain == "p182"] <- "5003"

root_length$Strain <- factor(root_length$Strain, levels=c("NB","cl21", "5003", "5005"))
#root_length$Strain <- factor(root_length$Strain, levels=c("cl21","NB", "deltafliC","5003", "5005"))


root_length$Rep <- factor(root_length$Rep, levels = c(1,2,3))
root_length$Root_growth..mm. <- as.numeric(as.character(root_length$Root_growth..mm.))
root_length$pid <- paste(root_length$Plate, root_length$Rep, sep = "_")
root_length$pid <- factor(root_length$pid)

ggplot(root_length, aes(Strain, Root_growth..mm.,color=Condition))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitterdodge())


library(ggpubr)
sum_data <- root_length %>%
  dplyr::group_by(Strain, Condition)%>%
  dplyr::summarise(mean_data = median(Root_growth..mm.),
                   se_data = sd(Root_growth..mm.)/(sqrt(n())),
                   count_data = n())

ggplot(root_length, aes(Strain, Root_growth..mm.))+
  geom_boxplot(outlier.colour = NA)+
  geom_point()+
  facet_wrap(facets = 'Condition')+
  geom_text(data=sum_data, aes(x=Strain, y=82, label=count_data))+
  #geom_errorbar(data=sum_data, aes(x=Strain, y=mean_data, ymin=mean_data-se_data, ymax=mean_data+se_data))+
  #stat_compare_means(ref.group = "cl21", method = "wilcox", label = "p.signif", size=5, label.y = 78)+
  labs(y="Main Root Growth (mm)")+
  theme_bw()

trial <- aov(Root_growth..mm. ~ Strain * Condition * Rep, data = root_length)
summary(trial)

library(MuMIn)
dredge(trial)

library(lme4)
library(lmerTest)
root_length$Strain <- factor(root_length$Strain, levels=c("cl21","NB", "5003", "5005"))
plus_flg_model <- lmer(Root_growth..mm. ~ Strain * Rep + (1|pid), data = root_length[root_length$Condition == "flg+",])
minus_flg_model <- lmer(Root_growth..mm. ~ Strain * Rep + (1|pid), data = root_length[root_length$Condition == "flg-",])

#significant random plate effect
summary(plus_flg_model)  
summary(minus_flg_model)

anova(plus_flg_model)

#THIS WILL STAY HERE FOR AFTER THE 3rd REP
#Lets identify groupings within our model using a Tukey post-hoc comparisons and adjusting p-values with the holm method
library(multcomp)
summary(glht(plus_flg_model, linfct = mcp(Strain = "Tukey")), test = adjusted("holm"))
summary(glht(minus_flg_model, linfct = mcp(Strain = "Tukey")), test = adjusted("holm"))
#i manually added the results to the figure
# All a single grouping in the flg- model
# a) Cl21 and 5003 b) NB and 5005 - flg+ groupings

