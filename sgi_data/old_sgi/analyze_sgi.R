#analyze the SGI data

library(tidyr)
library(ggplot2)

#THIS IS KATE'S DATA
fls2_sgi <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/sgi_data/kate_2run_sgi.csv", header = T, stringsAsFactors = F)

#this is Kate's second data
kate2_sgi <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/sgi_data/kate_second_run_sgi_1_19.csv", header = T, stringsAsFactors = F)

fls2_sgi <- fls2_sgi[,1:10]
fls2_sgi$experiment <- 1
#added for kates data
kate2_sgi$plate <- factor(kate2_sgi$plate, levels = 11:27)
fls2_sgi$plate <- factor(fls2_sgi$plate, levels = 1:10)

#go to long format
fls2_sgi <- gather(fls2_sgi, replicate, weight, rep1:rep8)
kate2_sgi <- gather(kate2_sgi, replicate, weight, rep1:rep8)

#get rid of NA's
fls2_sgi <- fls2_sgi[!is.na(fls2_sgi$weight),]
kate2_sgi <- kate2_sgi[!is.na(kate2_sgi$weight),]

#make the replicate a factor
fls2_sgi$replicate <- factor(fls2_sgi$replicate, levels = unique(fls2_sgi$replicate) )
kate2_sgi$replicate <- factor(kate2_sgi$replicate, levels = unique(kate2_sgi$replicate) )

#this handles a value that is obviously bad
fls2_sgi <- fls2_sgi[fls2_sgi$weight < .08,]

kate2_sgi$weight <- as.numeric(kate2_sgi$weight)

  
hist(fls2_sgi$weight, breaks = 30)
hist(kate2_sgi$weight[kate2_sgi$peptide=="Pa20"], breaks = 30)

#relevel the peptides so the comparisons are to mock
kate2_sgi$peptide[kate2_sgi$peptide == "Flg22 Pta "] = "Pta"
kate2_sgi$peptide[kate2_sgi$peptide == "Pta DA"] = "PtaDA"
kate2_sgi$peptide[kate2_sgi$peptide == "FLG20"] = "Pa20"
kate2_sgi$peptide[kate2_sgi$peptide == "FLG22 Pa"] = "Pa"

fls2_sgi$peptide <- factor(fls2_sgi$peptide)
kate2_sgi$peptide <- factor(kate2_sgi$peptide)

fls2_sgi$peptide <- relevel(fls2_sgi$peptide, "mock")
kate2_sgi$peptide <- relevel(kate2_sgi$peptide, "Mock ")

unique(kate2_sgi$peptide)

#compare the models
library(MuMIn)

reg_lm <- lm(weight ~ peptide * plate, fls2_sgi, na.action = na.omit)
reg_lm <- lm(weight ~ peptide + plate, fls2_sgi, na.action = na.omit)
drop1(reg_lm, test = "F")
reg_lm_2 <- lm(weight ~ peptide, fls2_sgi, na.action = na.omit)
model.sel(reg_lm_2, reg_lm) #this shows that the addition of plate does not really
# increase the performance of the lm

dat2_lm <- lm(weight ~ peptide * plate, kate2_sgi, na.action = na.omit)
drop1(dat2_lm, test = "F")
dat2_lm <- lm(weight ~ peptide + plate, kate2_sgi, na.action = na.omit)
drop1(dat2_lm, test = "F")
dat2_lm_2 <- lm(weight ~ peptide, kate2_sgi, na.action = na.omit)
model.sel(dat2_lm_2, dat2_lm)
#Here the simpler model actually performs better

####
# In both cases the simpler model performed better! Meaning plate effect is not significant
####

#make sure the data is normally distributed
require(car)
require(MASS)
require(stats)
#fit each type of distribution as the red bars and look for the one where the most data is within the dashed confidence intervals
#Normal distribution
qqp(fls2_sgi$weight, "norm")
qqp(kate2_sgi$weight, "norm")

#load the linear modeling package
require(lme4)
#by requiring lmerTest we overload lme4 to give additional features, 
#summary() will now include approximate degrees of freedom and p-values using the Satterthwaite approximation.
require(lmerTest)

#REML is set to falso if the data is normal
sgi.lmm <- lmer(weight ~ peptide + (1|plate), data = fls2_sgi, REML = F)
summary(sgi.lmm)
anova(sgi.lmm)

sgi.lmm.2 <- lmer(weight ~ peptide + (1|plate), data = kate2_sgi, REML = F)
summary(sgi.lmm.2)

require(LMERConvenienceFunctions)
#run a general check on the model
#the density plot should look roughly normal, centered on 0. In the quantile plot the black dots should roughly mirror the red line. And in the fitted residual plot most black dots should be between the red lines, with no obvious skew from left to right.
mcp.fnc(sgi.lmm)

#mixed linear model summary
t <- summary(dat2_lm)
t <- summary(reg_lm)

s <- as.data.frame( t$coefficients )
colnames(s)[2] <- "StdE"
s <- s[2:nrow(s),]
s_high <- s[s$`Pr(>|t|)` < .001,]
s_med <- s[s$`Pr(>|t|)` < .01 & s$`Pr(>|t|)` > .001,]
s_low <- s[s$`Pr(>|t|)` < .05 & s$`Pr(>|t|)` > .01,]
p <- ggplot(s, aes(x=row.names(s), y=Estimate))+
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=Estimate-StdE, ymax=Estimate+StdE))+
  geom_text(data = s_med, aes(row.names(s_med)), label = "**", vjust=-1.5)+
  geom_text(data = s_high, aes(row.names(s_high)), label = "***", vjust=-1.5)+
  geom_text(data = s_low, aes(row.names(s_low)), label = "*", vjust =-1.5)
p
#############################
####  order based on ROS ####
############################# 

pos_control <- "Pta"
neg_control <- "PtaDA"

#This path is to all the directories containing the results from the first experiments
first_path_ros <- "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/ros_8_27-30/"

#this path is to the ROS results from the second wave of ROS bursts performed
new_path_ros <- "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/inhib_10_20s_2018/"

#this path is to the two new peptides that I performed ROS on
alan_path <- "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/alan_lab/ros_12_12_2018/"


#here I assume that the first path separates the files by dates, while second path contains all the files regardless of data
get_files_from_directories <- function( first_path_ext, new_path_ext, alan ) {
  #this portion of the function retrieves the old dataset
  file_list <- c()
  num <- c()
  for ( i in 27:31) {
    mm_files <- list.files(path=paste(first_path_ext, "8_", i, "_18", sep = ""), pattern = "*auc_data.csv")
    file_list <- append(file_list, mm_files)
    num <- append(num, rep(i, length(mm_files)))
  }
  csv_list <- list()
  for ( i in 1:length(file_list)) {
    csv_list[[i]] <- read.csv(paste(first_path_ext, "8_", num[i], "_18/", file_list[i], sep = ""), stringsAsFactors = F)
  }
  
  #this portion of code retrieves the new dataset
  file_list <- list.files(path = new_path_ext, pattern = "*auc_data.csv")
  
  prev_point <- length(csv_list) #keeps track of where the old data stops
  
  for ( i in 1:length(file_list) ) {
    csv_list[[i+prev_point]] <- read.csv(paste(new_path_ext, file_list[i], sep = ""), stringsAsFactors = F)
  }
  
  file_list <- list.files(path = alan, pattern = "*auc_data.csv")
  new_point <- length(csv_list) #keeps track of where the old data stops
  
  for ( i in 1:length(file_list) ) {
    csv_list[[i+new_point]] <- read.csv(paste(alan, file_list[i], sep = ""), stringsAsFactors = F)
  }
  return(list(csv_list, prev_point) )
}

#this puts all the data into one dataframe
get_the_combined_dataframe <- function( csv, sep_point ) {
  for ( i in 1:length(csv) ) {
    csv[[i]]$plate <- i
    
    if ( i > sep_point ) {
      csv[[i]]$batch <- 2
    }
    else {
      csv[[i]]$batch <- 1
    }
    
    if ( i > 1 ) {
      csv[[1]] <- rbind(csv[[1]], csv[[i]])
    }
  }
  return(csv[[1]])
}

get_data_frame_rdy_4_plting <- function(df) {
  #before plotting I need to combine the 236 samples, which were run on the same plate
  df$Name[df$Name == "236_a"] <- "236"
  #get rid of D13
  df <- df[df$Name != "D13",]
  #get rid of distilled water samples
  df <- df[df$Name != "DW" & df$Name != "dw",]
  #get rid of ladder samples
  df <- df[!grepl("nm", df$Name),]
  
  df$AUC <- log10(df$AUC)
  df$plate <- factor(as.character(df$plate), levels = as.character(1:length(unique(df$plate)))) 
  df$batch <- factor(as.character(df$batch), levels = as.character(1:2), labels = c("Batch_1", "Batch_2") )
  #make it so pta is first
  df$Name <- factor(as.character(df$Name))
  df$Name <- relevel(df$Name, "Pta")
  
  return(df)
}


### MAIN ###
csv_part <- get_files_from_directories(first_path_ros, new_path_ros, alan_path)
data_ros <- get_the_combined_dataframe(csv_part[[1]], csv_part[[2]])
data_ros <- get_data_frame_rdy_4_plting(data_ros)

ptada_means_std <- c() #THIS IS ACTAULLY MEDIAN******
for ( i in 1:nlevels(data_ros$plate) ) {
  ptada_means_std <- append(ptada_means_std, c(i, mean(data_ros$AUC[data_ros$plate == i & data_ros$Name == neg_control]), sd(data_ros$AUC[data_ros$plate == i & data_ros$Name == neg_control])))
}
ptada_means_std <- as.data.frame(matrix(ptada_means_std, ncol = 3, byrow = T))

z_scores_per_peptide <- c()
for ( i in ptada_means_std[,1] ) {
  for ( j in unique(data_ros$Name[data_ros$plate == i]) ) {
    z_scores_per_peptide <- append(z_scores_per_peptide, c(j, (mean(data_ros$AUC[data_ros$Name == j & data_ros$plate == i ]) - ptada_means_std$V2[i])/ptada_means_std$V3[i] ) )
  }
}
z_scores_per_peptide <- as.data.frame(matrix(z_scores_per_peptide, ncol = 2, byrow = T), stringsAsFactors = F)
z_scores_per_peptide$V2 <- as.numeric(z_scores_per_peptide$V2)

ggplot(z_scores_per_peptide, aes(x=reorder(V1,-V2, median), y=V2))+
  geom_boxplot( outlier.shape = NA)+
  ylab("Log10 AUC PtaDA Z-Score")+
  xlab("Peptide ID")+
  labs(color="Plate Number")+
  ggtitle("ROS Burst of Natural Variants")+
  theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))

####################
#### FINAL PLOT ####
####################

t <- summary(dat2_lm)
t <- summary(reg_lm)

s <- as.data.frame( t$coefficients )
colnames(s)[2] <- "StdE"
s <- s[2:nrow(s),]
row.names(s) <- gsub("peptide","", row.names(s))

order_names <- reorder(z_scores_per_peptide$V1,-z_scores_per_peptide$V2, median)
names_in_sgi <- levels(order_names)[levels(order_names) %in% row.names(s)]

s <- s[row.names(s) %in% names_in_sgi,]
s_high <- s[s$`Pr(>|t|)` < .001,]
s_med <- s[s$`Pr(>|t|)` < .01 & s$`Pr(>|t|)` > .001,]
s_low <- s[s$`Pr(>|t|)` < .05 & s$`Pr(>|t|)` > .01,]
s <- s[row.names(s)[match(names_in_sgi, row.names(s))],]


p <- ggplot(s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate))+
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=Estimate-StdE, ymax=Estimate+StdE))+
  geom_text(data = s_med, aes(row.names(s_med)), label = "**", vjust=-1.5)+
  geom_text(data = s_high, aes(row.names(s_high)), label = "***", vjust=-1.5)+
  geom_text(data = s_low, aes(row.names(s_low)), label = "*", vjust =-1.5)
  

fls2_sgi <- fls2_sgi[fls2_sgi$peptide %in% row.names(s) | fls2_sgi$peptide == "mock",]
fls2_sgi$peptide <- factor(fls2_sgi$peptide, levels=c(row.names(s),"mock"))
p <- ggplot(fls2_sgi, aes(x=peptide, y=weight))+
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=14))+
  ylim(0,0.045)+
  geom_text(data = s_med, aes(row.names(s_med), y=.04), label = "**", vjust=-1.5, inherit.aes = F)+
  geom_text(data = s_high, aes(row.names(s_high), y=.04), label = "***", vjust=-1.5,inherit.aes = F)+
  #geom_text(data = s_low, aes(row.names(s_low),y=.04), label = "*", vjust =-1.5,inherit.aes = F)+
  ylab("Fresh Weight (grams)") +
  xlab("Peptide ID")
p

kate2_sgi <- kate2_sgi[kate2_sgi$peptide %in% row.names(s) | kate2_sgi$peptide == "Mock ",]
kate2_sgi$peptide <- factor(kate2_sgi$peptide, levels=c(row.names(s),"Mock "))
unique(kate2_sgi$peptide)
p <- ggplot(kate2_sgi, aes(x=peptide, y=weight))+
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=14))+
  #ylim(0,0.045)+
  geom_text(data = s_med, aes(row.names(s_med), y=1), label = "**", vjust=-1.5, inherit.aes = F)+
  geom_text(data = s_high, aes(row.names(s_high), y=1), label = "***", vjust=-1.5,inherit.aes = F)+
  #geom_text(data = s_low, aes(row.names(s_low),y=.04), label = "*", vjust =-1.5,inherit.aes = F)+
  ylab("Fresh Weight (mg)") +
  xlab("Peptide ID")
p



ggplot(z_scores_per_peptide[z_scores_per_peptide$V1 %in% kate2_sgi$peptide,], aes(x=reorder(V1,-V2, median), y=V2))+
  geom_boxplot( outlier.shape = NA)+
  ylab("Log10 AUC PtaDA Z-Score")+
  xlab("Peptide ID")+
  labs(color="Plate Number")+
  #ggtitle("ROS Burst of Natural Variants")+
  theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))
