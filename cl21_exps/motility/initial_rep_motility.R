#motility

motility_r1 <- read.csv('/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/CL21_colonization_experiments/motility/4reps_CL21_motility_assay_Day_4.csv')

motility_r1 <- motility_r1[motility_r1$bacteria %in% c('p176 D8', "p182 H5-1", "delta flic", "CL21 wt"),]
motility_r1$bacteria <- as.character(motility_r1$bacteria)
motility_r1$bacteria[motility_r1$bacteria == "p176 D8"] <- "5005"
motility_r1$bacteria[motility_r1$bacteria == "p182 H5-1"] <- "5003"
motility_r1$bacteria[motility_r1$bacteria == "delta flic"] <- "delta fliC"
motility_r1$bacteria[motility_r1$bacteria == "CL21 wt"] <- "WT"

motility_r1$bacteria <- factor(motility_r1$bacteria, levels = c("delta fliC", "WT", '5003', "5005"))
#motility_r1 <- motility_r1[motility_r1$Strain != "5006" & motility_r1$Strain != "Pa22",]
colnames(motility_r1) <- c("Rep", "Strain", "max_swim")
motility_r1$Rep <- factor(motility_r1$Rep)
 
#percent of motility lost
(mean(motility_r1$max_swim[motility_r1$Strain == "WT"]) - mean(motility_r1$max_swim[motility_r1$Strain == "delta fliC"]))/mean(motility_r1$max_swim[motility_r1$Strain == "WT"]) 

(mean(motility_r1$max_swim[motility_r1$Strain == "WT"]) - mean(motility_r1$max_swim[motility_r1$Strain == "5005"]))/mean(motility_r1$max_swim[motility_r1$Strain == "WT"]) 
(mean(motility_r1$max_swim[motility_r1$Strain == "WT"]) - mean(motility_r1$max_swim[motility_r1$Strain == "5003"]))/mean(motility_r1$max_swim[motility_r1$Strain == "WT"]) 

library(dplyr)

summary_data <- motility_r1 %>%
  dplyr::group_by(Strain) %>%
  dplyr::summarise(median_diam = median(max_swim),
                   total_count = n())

ggplot(motility_r1, aes(Strain, max_swim))+
  geom_point(size=3, position=position_jitter(width = .1))+
  labs(y="Swarming Max Diameter (cm)")+
  theme_bw(base_size = 16)

stats <- aov(max_swim ~ Strain * Rep, data = motility_r1)
summary(stats)
#don't need rep effect
stats <- aov(max_swim ~ Strain, data = motility_r1)

library(agricolae)

tukey <- HSD.test(stats, trt="Strain")
