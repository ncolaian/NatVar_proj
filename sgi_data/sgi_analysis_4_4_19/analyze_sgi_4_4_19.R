#This script will analyze all of kate's sgi data with a linear mixed model

library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)

#load up the file
full_sgi <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/sgi_data/sgi_analysis_4_4_19/full_kate_sgi_data.csv")
full_sgi$W2 <- as.numeric(as.character(full_sgi$W2))
full_sgi$W6 <- as.numeric(as.character(full_sgi$W6))
full_sgi <- full_sgi[,1:12]

#Get the ROS burst data
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

#I need to change the format from wide to long
library(tidyr)

full_sgi <- gather(full_sgi, "plate_col", "gramsdtenthousand", W1:W8, na.rm = T)
full_sgi$Plate <- factor(full_sgi$Plate)
full_sgi$Exp <- factor(full_sgi$Exp)

colo_sgi <- full_sgi[full_sgi$Plant == "colo",]
colo_sgi$ID <- as.character(colo_sgi$ID)
colo_sgi$ID[colo_sgi$ID == "flg22"] <- "Pta" 
colo_sgi$ID <- factor(colo_sgi$ID)
colo_sgi$ID <- relevel(colo_sgi$ID, ref = "mock")
colo_sgi$gramsdtenthousand <- colo_sgi$gramsdtenthousand/10

ggplot(na.omit(colo_sgi[colo_sgi$ID == "mock",]), aes(Plate, gramsdtenthousand, color=Exp))+
  geom_boxplot(outlier.shape = NA)

#i am using a nested linear mixed model where I am nesting the plates within experiment
#I am only going to use the col-o data
full_model <- lmer(gramsdtenthousand ~ ID + (1|Exp/Plate), data = colo_sgi)

#going to compare this to a fixed exp effect and a random plate effect
less_rand <- lmer(gramsdtenthousand ~ ID + Exp + (1|Plate), data = colo_sgi)

summary(full_model)

anova(full_model, less_rand)

### The less random model fit much better ###

summary(less_rand)

ggplot(colo_sgi, aes(ID, gramsdtenthousand)) +
  geom_boxplot(outlier.shape = NA)

#use the functions to get the ROS
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

order_names <- reorder(z_scores_per_peptide$V1,-z_scores_per_peptide$V2, median)

#plot the results
t <- summary(less_rand)
s <- as.data.frame( t$coefficients )
colnames(s)[2] <- "StdE"
s_inter <- s[1,]
s <- s[2:(nrow(s)-4),]
row.names(s) <- gsub("ID","", row.names(s))
#s$Estimate <- s$Estimate + s_inter$Estimate[1]

names_in_sgi <- levels(order_names)[levels(order_names) %in% row.names(s)]
s <- s[row.names(s)[match(names_in_sgi, row.names(s))],]

s_high <- s[s$`Pr(>|t|)` < .001,]
s_med <- s[s$`Pr(>|t|)` < .01 & s$`Pr(>|t|)` > .001,]
s_low <- s[s$`Pr(>|t|)` < .05 & s$`Pr(>|t|)` > .01,]

#raw data with code
colo_sgi$ID <- as.character(colo_sgi$ID)
colo_sgi$ID[colo_sgi$ID == "flg22"] <- "Pta" 
colo_sgi$ID <- factor(colo_sgi$ID, levels=c(row.names(s),"mock"))
mock_med <- median(na.omit(colo_sgi$gramsdtenthousand[colo_sgi$ID == "mock"]))

#effects
p <- ggplot(s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate))+
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))+
  geom_errorbar(aes(ymin=Estimate-StdE, ymax=Estimate+StdE))+
  geom_text(data = s_med, aes(row.names(s_med), y=9), label = "**", vjust=-1.5)+
  geom_text(data = s_high, aes(row.names(s_high), y=9), label = "***", vjust=-1.5)+
  geom_text(data = s_low, aes(row.names(s_low), y=9), label = "*", vjust =-1.5)
p
p+geom_point(data = na.omit(colo_sgi), aes(x=ID, y=gramsdtenthousand-mock_med))

#Final Plot
p <- ggplot(na.omit(colo_sgi), aes(x=ID, y=gramsdtenthousand-mock_med))+
  geom_point(outlier.shape = NA, color="grey") +
  #geom_jitter(aes(fill=Plate))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=14))+
  ylim(-25,40)+
  geom_text(data = s_med, aes(row.names(s_med), y=37), label = "**", vjust=-1.5, inherit.aes = F)+
  geom_text(data = s_high, aes(row.names(s_high), y=37), label = "***", vjust=-1.5,inherit.aes = F)+
  geom_text(data = s_low, aes(row.names(s_low),y=37), label = "*", vjust =-1.5,inherit.aes = F)+
  ylab("Mock Centered Fresh Weight (mg)") +
  xlab("Peptide ID")
p <- p + geom_errorbar(data=s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate, ymin=Estimate-StdE, ymax=Estimate+StdE))
p

#get the sgi data into a a format to add to itol
get_seq <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/youssef_data/old_arrays/finalized_useq_8_22_2018.csv", stringsAsFactors = F)
#The chinese peptides had numbering that I placed uniquely. I would like to get the sequenced based on this file
get_seq_2 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Kate_Vienna_data/peptides_to_order_china.csv")

wheel_data <- c()
for ( i in row.names(s) ) {
  if ( as.character(i) %in% as.character( get_seq$UseqId ) ) {
    stri <- get_seq$Peptides[get_seq$UseqId == i]
  }
  else if ( as.character(i) %in% as.character( get_seq_2$ID ) ) {
    stri <- as.character(get_seq_2$Peptides.to.order[get_seq_2$ID == i])
  }
  else {
    stri <- "A"
  }
  print(stri)
  wheel_data <- append(wheel_data, c(stri, i, as.numeric(s$Estimate[row.names(s) == i]), as.numeric(s$`Pr(>|t|)`[row.names(s) == i]) ))
}
wheel_data <- matrix(wheel_data, ncol = 4, byrow = T)

write.csv(wheel_data,file = "/Users/nicholascolaianni/Documents/dangl_lab/882_data/updated_work_10_23_18/new_clusters/all_nat_variants_sgieffect_wheel_5_23_19.csv" , quote = F, row.names = F)


###############################
#Analyze the fls2 data 
fls2_sgi <- full_sgi[full_sgi$Plant=="fls2efr",]
fls2_sgi$gramsdtenthousand <- fls2_sgi$gramsdtenthousand/10

ggplot(fls2_sgi, aes(x=ID, y=gramsdtenthousand)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means()

#######################################
#Just get the experiment with fls2/efr and then make the correct comparisons
plus_comparisons <- full_sgi[full_sgi$Exp == "20190221",]
plus_comparisons$gramsdtenthousand <- plus_comparisons$gramsdtenthousand/10

#get the rows with a peptide with a positive significant effect
signif_pos_peptides <- c()
for ( i in row.names(rbind(s_high, s_med)) ) {
  if ( i %in% row.names(s_high)) {
    if ( s_high[row.names(s_high) == i,][1] > 0 ) {
      signif_pos_peptides <- append(signif_pos_peptides, i)
    }
  }
  else{
    if ( s_med[row.names(s_med) == i,][1] > 0 ) {
      signif_pos_peptides <- append(signif_pos_peptides, i)
    }
  }
}

plus_comparisons <- plus_comparisons[plus_comparisons$ID %in% c(signif_pos_peptides,"mock"),]

ggplot(plus_comparisons, aes(Plant, gramsdtenthousand))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(facets="ID")+
  stat_compare_means(method="anova", label = "p.signif", label.x = 1.5)

for ( i in unique(plus_comparisons$Plate )) {
  p_mock_med <- median(plus_comparisons$gramsdtenthousand[plus_comparisons$ID == "mock" & 
                                                            plus_comparisons$Plate == i])
  
  plus_comparisons$gramsdtenthousand[plus_comparisons$Plate==i] <- plus_comparisons$gramsdtenthousand[plus_comparisons$Plate==i] - p_mock_med
}

ggplot(plus_comparisons, aes(Plant, gramsdtenthousand))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(facets="ID")+
  stat_compare_means(method = "anova") +
  ylab("Mock Centered Fresh Weight (mg)") +
  xlab("Plant")
  