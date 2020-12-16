#this script will be used to create a graph for all the ROS burst assays performed so far. It will also be used
# to explore the plate effects and the possibility of using a glmm to normalize the results
# I will also be looking to the pta results to look at the plate effect and if there's any significance

library(ggplot2)

pos_control <- "Pta"
neg_control <- "PtaDA"

#This path is to all the directories containing the results from the first experiments
first_path_ros <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/ros_8_27-30/"

#this path is to the ROS results from the second wave of ROS bursts performed
new_path_ros <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/inhib_10_20s_2018/"

#this path is to the two new peptides that I performed ROS on
alan_path <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/alan_lab/ros_12_12_2018/"

#This path is to the ROS data from the peptides ordered from China
china_path <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/ch_100_5_9_19/"

#potential antag
pot_antag <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/100nm_ros_potantag_10_16_19/"

#here I assume that the first path separates the files by dates, while second path contains all the files regardless of data
get_files_from_directories <- function( first_path_ext, new_path_ext, alan, china_p, pot_a ) {
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
  file_list <- list.files(path = china_p, pattern = "*auc_data.csv")
  new_point_2 <- length(csv_list) #keeps track of where the old data stops
  
  for ( i in 1:length(file_list) ) {
    csv_list[[i+new_point_2]] <- read.csv(paste(china_p, file_list[i], sep = ""), stringsAsFactors = F)
  }
  
  file_list <- list.files(path = pot_a, pattern = "*auc_data.csv")
  new_point_3 <- length(csv_list) #keeps track of where the old data stops
  
  for ( i in 1:length(file_list) ) {
    csv_list[[i+new_point_3]] <- read.csv(paste(pot_a, file_list[i], sep = ""), stringsAsFactors = F)
  }
  
  return(list(csv_list, c(prev_point, new_point, new_point_2, new_point_3)) )
}

#this puts all the data into one dataframe
get_the_combined_dataframe <- function( csv, sep_point ) {
  for ( i in 1:length(csv) ) {
    csv[[i]]$plate <- i
    if ( i > sep_point[4] ) {
      csv[[i]]$batch <- 5
    }
    else if ( i > sep_point[3] ) {
      csv[[i]]$batch <- 4
    }
    else if ( i > sep_point[2] ) {
      csv[[i]]$batch <- 3
    }
    else if ( i > sep_point[1] ) {
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
  df$Name[df$Name == "1447"] <- "447"
  df <- df[df$Name != "dud",]
  df <- df[df$Name != "dud1",]
  df <- df[df$Name != "dud2",]
  df <- df[df$Name != "dud3",]
  
  df <- na.omit(df)
  
  
  
  df$Name[df$Name == "pta"] <- "Pta"
  df$Name[df$Name == "1787"] <- "Pta"
  df$Name[df$Name == "ptada"] <- "PtaDA"
  df$Name[df$Name == "ptaDA"] <- "PtaDA"
  df$Name[df$Name == "pa"] <- "Pa"
  df$Name[df$Name == "pa20"] <- "Pa20"
  
  df$AUC <- log10(df$AUC)
  df$plate <- factor(as.character(df$plate))
  df$batch <- factor(as.character(df$batch), levels = as.character(1:5), labels = c("Batch_1", "Batch_2", "Batch_3", "Batch_4", "Batch_5") )
  #make it so pta is first
  df$Name <- factor(as.character(df$Name))
  df$Name <- relevel(df$Name, "Pta")
  
  return(df)
}
### MAIN ###

#ROS
csv_part <- get_files_from_directories(first_path_ros, new_path_ros, alan_path, china_path, pot_antag)
data_ros <- get_the_combined_dataframe(csv_part[[1]], csv_part[[2]])
data_ros <- get_data_frame_rdy_4_plting(data_ros)
data_ros <- na.omit(data_ros)

unique(data_ros$batch)
#batch effect plot
ggplot(data_ros, aes(x=reorder(Name,-AUC, median), y=AUC-median(data_ros$AUC[data_ros$Name == neg_control])))+
  geom_boxplot(outlier.shape = NA)+
  ylab("Log10 AUC")+
  xlab("Peptide ID")+
  ggtitle("ROS Values Per Petide")+
  theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))

ggplot(data_ros[data_ros$Name == "Pta",], aes(x=batch, y=AUC))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  ylab("Log10 AUC")+
  xlab("ROS Batch")+
 #labs(color="Plate Number")+
  ggtitle("Pta 100nm Batch Effect")+
  theme(text = element_text(size=20))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = .5))

#simple t-test between batches
t.test(data_ros$AUC[data_ros$batch == "Batch_1"], data_ros$AUC[data_ros$batch == "Batch_2"])

#plate and batch effect
ggplot(data_ros[data_ros$Name == "Pa",], aes(x=batch, y=AUC, color=plate))+
  geom_boxplot(outlier.shape = NA)+
  ylab("Log10 AUC")+
  xlab("ROS Batch")+
  labs(color="Plate Number")+
  ggtitle("PtaDA 100nm Batch and Plate Effect")+
  theme(text = element_text(size=20))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = .5))

#look at specific data
ggplot(data_ros[data_ros$batch == "Batch_4",], aes(x=Name, y=AUC))+
  geom_boxplot(outlier.shape = NA)+
  ylab("Log10 AUC")+
  xlab("ROS Batch")+
  labs(color="Plate Number")+
  ggtitle("PtaDA 100nm Batch and Plate Effect")+
  #theme(text = element_text(size=20))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = .5))

#ussing a glmm like in the r file kate sent me
require(car)
require(MASS)
require(stats)

qqp(data_ros$AUC[data_ros$Name == "Pta"], "norm")
hist(data_ros$AUC[data_ros$Name == "Pta"])
hist(data_ros$AUC[data_ros$Name == "Pta" & data_ros$batch == "Batch_1"])


#try glmm
require(lme4)
#by requiring lmerTest we overload lme4 to give additional features, 
#summary() will now include approximate degrees of freedom and p-values using the Satterthwaite approximation.
require(lmerTest)

#first fit a linear mixed model, with Trial as the random effect
#set REML to FALSE if the data is normal, so we can use the maximum likelihood because the data are normal, we have similar sample sizes between effects, and we only have one random effect
anova(lm(AUC ~ Name + plate, data = data_ros))
ROS.lmm <- lmer(AUC ~ Name + (1 | plate), data = data_ros, REML = FALSE)
summary(ROS.lmm)
anova(ROS.lmm)
#lets check the model assumptions
#load the package
require(LMERConvenienceFunctions)
#run a general check on the model
#the density plot should look roughly normal, centered on 0. In the quantile plot the black dots should roughly mirror the red line. And in the fitted residual plot most black dots should be between the red lines, with no obvious skew from left to right.
mcp.fnc(ROS.lmm)
hist(residuals(ROS.lmm))
qqnorm(residuals(ROS.lmm))
qqline(residuals(ROS.lmm), col="red")
#If the above looks good, then we can get the model characteristics
summary(ROS.lmm)
fixef(ROS.lmm)
ranef(ROS.lmm)
fitted(ROS.lmm)
nrow(data_ros)
length(fitted(ROS.lmm))

data_ros$fit_val <- fitted(ROS.lmm)

rand_ef <- ranef(ROS.lmm)
rand_ef <- rand_ef[[1]]

for ( i in 1:nrow(data_ros) ) {
  data_ros$rremoved[i] <- data_ros$AUC[i] - rand_ef[data_ros$plate[i],1]#,data_ros$batch[i]]
}

ggplot(data_ros[data_ros$Name == "PtaDA",], aes(x=batch, y=rremoved, color=plate))+
  geom_boxplot(outlier.shape = NA)+
  ylab("Log10 AUC")+
  xlab("ROS Batch")+
  labs(color="Plate Number")+
  ggtitle("Pta 100nm Batch and Plate Effect")+
  theme(text = element_text(size=20))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = .5))


ggplot(data_ros[data_ros$Name == "Pta",], aes(x=batch, y=rremoved))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  ylab("Log10 AUC")+
  xlab("ROS Batch")+
  #labs(color="Plate Number")+
  ggtitle("Pta 100nm Batch Effect")+
  theme(text = element_text(size=20))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = .5))

#simple t-test between batches
t.test(data_ros$rremoved[data_ros$batch == "Batch_1"], data_ros$rremoved[data_ros$batch == "Batch_2"])

ggplot(data_ros, aes(x=Name, y=rremoved))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter( aes(color = plate))+
  ylab("Log10 AUC")+
  xlab("Peptide ID")+
  labs(color="Plate Number")+
  theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Get sequences
get_seq <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/youssef_data/old_arrays/finalized_useq_8_22_2018.csv", stringsAsFactors = F)
wheel_data <- c()
for ( i in unique(data_ros$Name) ) {
  if ( as.character(i) %in% as.character( get_seq$UseqId ) ) {
    stri <- get_seq$Peptides[get_seq$UseqId == i]
  }
  else {
    stri <- "A"
  }
  print(stri)
  wheel_data <- append(wheel_data, c(stri, i, median(data_ros$AUC[data_ros$Name == i])))
}
wheel_data <- matrix(wheel_data, ncol = 3, byrow = T)

#write.csv(wheel_data,file = "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/all_nat_variants_analysis/median_vals.csv" , quote = F, row.names = F)

#### Z-score normalization ####
#Using meta-analysis - ie combining the results from multiple tests
#going to calculate the z-score of each peptide against ptada

ptada_means_std <- c() 
for ( i in 1:nlevels(data_ros$plate) ) {
  ptada_means_std <- append(ptada_means_std, c(i, mean(data_ros$AUC[data_ros$plate == i & data_ros$Name == neg_control]), sd(data_ros$AUC[data_ros$plate == i & data_ros$Name == neg_control])))
}
ptada_means_std <- as.data.frame(matrix(ptada_means_std, ncol = 3, byrow = T))

#**Talk to Corbin about this**
z_scores_per_peptide <- c()
for ( i in ptada_means_std[,1] ) {
  for ( j in unique(data_ros$Name[data_ros$plate == i]) ) {
    for (f in data_ros$AUC[data_ros$Name == j & data_ros$plate == i ]) {
      z_scores_per_peptide <- append(z_scores_per_peptide, c(j, (f - ptada_means_std$V2[i])/ptada_means_std$V3[i] ) )
    }
  }
}
z_scores_per_peptide <- as.data.frame(matrix(z_scores_per_peptide, ncol = 2, byrow = T), stringsAsFactors = F)
z_scores_per_peptide$V2 <- as.numeric(z_scores_per_peptide$V2)

#MAIN FIGURE
nv_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/final_figure_antag/supplemental_data/supplemental_table2.csv")

ggplot(z_scores_per_peptide[(z_scores_per_peptide$V1 %in% as.character(nv_data$PeptideID)),], aes(x=reorder(V1,-V2, median), y=V2))+
  geom_boxplot( outlier.shape = NA)+
  ylab("Log10 AUC PtaDA Z-Score")+
  xlab("Peptide ID")+
  labs(color="Plate Number")+
  ggtitle("ROS Burst of Natural Variants")+
  theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))+
  ylim(c(-10, 18))

get_seq <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/youssef_data/old_arrays/finalized_useq_8_22_2018.csv", stringsAsFactors = F)
#The chinese peptides had numbering that I placed uniquely. I would like to get the sequenced based on this file
get_seq_2 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Kate_Vienna_data/peptides_to_order_china.csv")
#Need to also put the peptides that are from the pot_antag
get_seq_3 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/882_data/updated_work_10_23_18/new_clusters/potential_antag_ralstonia.csv")

wheel_data <- c()
for ( i in unique(z_scores_per_peptide$V1) ) {
  if ( as.character(i) %in% as.character( get_seq$UseqId ) ) {
    stri <- get_seq$Peptides[get_seq$UseqId == i]
  }
  else if ( as.character(i) %in% as.character( get_seq_2$ID ) ) {
    stri <- as.character(get_seq_2$Peptides.to.order[get_seq_2$ID == i])
  }
  else if ( as.character(i) %in% as.character( get_seq_3$Peptide.ID ) ) {
    stri <- as.character(get_seq_3$Sequence[get_seq_3$Peptide.ID == i])
  }
  else if (as.character(i) == "Pta") {
    stri <- "TRLSSGLKINSAKDDAAGLQIA"
  }
  else if (as.character(i) == "Pa") {
    stri <- "QRLSTGSRINSAKDDAAGLQIA"
  }
  else if (as.character(i) == "PtaDA") {
    stri <- "TRLSSGLKINSAKADAAGLQIA"
  }
  else if (as.character(i) == "Pa20") {
    stri <- "QRLSTGSRINSAKDDAAGLQ"
  }
  else {
    stri <- "A"
  }
  print(stri)
  wheel_data <- append(wheel_data, c(stri, i, median(as.numeric(z_scores_per_peptide$V2[z_scores_per_peptide$V1 == i]))))
}
wheel_data <- matrix(wheel_data, ncol = 3, byrow = T)

write.csv(wheel_data,file = "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/882_data/updated_work_10_23_18/new_clusters/all_nat_variants_ptanorm_wheel_11_12_19.csv" , quote = F, row.names = F)

###############################
### THE REST IS DEPRICATESD ###
###############################

#stoufferz z-score method
stoufferz_zscores <- c()
for ( i in unique(z_scores_per_peptide$V1) ) {
  k <- length(z_scores_per_peptide$V2[z_scores_per_peptide$V1 == i])
  stoufferz_zscores <- append(stoufferz_zscores, c(i, (sum(z_scores_per_peptide$V2[z_scores_per_peptide$V1==i])/sqrt(k))) )
}
stoufferz_zscores <- as.data.frame(matrix(stoufferz_zscores, ncol = 2, byrow = T), stringsAsFactors = F)
stoufferz_zscores$V2 <- as.numeric(stoufferz_zscores$V2)
ggplot(stoufferz_zscores, aes(x=reorder(V1,-V2), y=V2))+
  geom_bar(stat = "identity")+
  ylab("Z-Score")+
  xlab("Peptide ID")+
  labs(color="Plate Number")+
  ggtitle("ROS Stouffer's Z-scores")+
  theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))


## THIS WAS ADDED TO SUBSAMPLE MY DATA WITH KATE"S GRAPH ##
sgi_peptides <- c(1195,264,1843,1090,"Pa", "Pta", 1477, 166, 17, 1391, 468, 470, 1943, 1213, 1128, 295, 1635, 94, 289, 1184, 1814, 1471, 1042, 1877, 1544, 1294, "PtaDA", 161, 236, 1186, 108, 359, 457, 59, 133, 118, 459, 20, 255, 1952, 321, 99, 146, 102, 1776, 1297, "Pa20")
z_scores_per_peptide <- z_scores_per_peptide[ z_scores_per_peptide$V1 %in% sgi_peptides ,]
ggplot(z_scores_per_peptide, aes(x=reorder(V1,-V2, median), y=V2))+
  geom_boxplot( outlier.shape = NA)+
  ylab("Log10 AUC PtaDA Z-Score")+
  xlab("Peptide ID")+
  labs(color="Plate Number")+
  ggtitle("ROS Burst of Natural Variants")+
  theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))

####THIS WILL GET A NON-PLATE SPECIFIC Z-Score
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

#I quickly would like to look at the correlation between identity to canonical and signal


