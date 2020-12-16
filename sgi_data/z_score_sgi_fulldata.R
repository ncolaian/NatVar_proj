# I am going to analyze all of the sgi data together with z-score
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(tidyr)

#load up the file
full <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/sgi_data/sgi_analysis_4_4_19/full_kate_sgi_data.csv")
full$W2 <- as.numeric(as.character(full$W2))
full$W6 <- as.numeric(as.character(full$W6))
full <- full[,1:12]


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
set.seed(1994)
#ROS
csv_part <- get_files_from_directories(first_path_ros, new_path_ros, alan_path, china_path, pot_antag)
data_ros <- get_the_combined_dataframe(csv_part[[1]], csv_part[[2]])
data_ros <- get_data_frame_rdy_4_plting(data_ros)
data_ros <- na.omit(data_ros)
#SGI pre-china
full <- gather(full, "plate_col", "gramsdtenthousand", W1:W8, na.rm = T)
full$gramsdtenthousand <- full$gramsdtenthousand/10

#SGI CH-apeptide peptides
full_sgi <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj//sgi_data/china_peptides_sgi.csv")
full_sgi <- gather(full_sgi, "plate_col", "mg", A:H, na.rm = T)
full_sgi$mg <- full_sgi$mg/10
colnames(full_sgi) <- c("ID", "Plate", "Exp","Plant","plate_col", "gramsdtenthousand")
full_sgi <- rbind(full_sgi[,c(1,2,3,4,6)], full[,c(1,2,3,4,6)])

#SGI Potential antag peptides
full_sgi2 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj//sgi_data/potential_antag_vienna.csv")
full_sgi2 <- gather(full_sgi2, "plate_col", "mg", A:H, na.rm = T)
full_sgi2$mg[full_sgi2$mg < 1 ] <- full_sgi2$mg[full_sgi2$mg < 1 ] * 10000
full_sgi2$mg <- full_sgi2$mg/10
colnames(full_sgi2) <- c("ID", "Plate", "Exp","Plant","plate_col", "gramsdtenthousand")
full_sgi <- rbind(full_sgi, full_sgi2[,c(1,2,3,4,6)])

full_sgi2$Plant
full_sgi$plate_exp <- paste(full_sgi$Plate, full_sgi$Exp, sep = "_")

#have to fix all the different spellings
full_sgi$ID <- as.character(full_sgi$ID)
full_sgi$ID[full_sgi$ID == "Mock"] <- "mock"
full_sgi$ID[full_sgi$ID == "Mock "] <- "mock"
full_sgi$ID[full_sgi$ID == "mock "] <- "mock"
full_sgi$ID[full_sgi$ID == "flg22 "] <- "Pa"
full_sgi$ID[full_sgi$ID == "flg22"] <- "Pa"
full_sgi$ID[full_sgi$ID == "pa"] <- "Pa"
full_sgi$ID[full_sgi$ID == "pa20"] <- "Pa20"
full_sgi$ID[full_sgi$ID == "pta"] <- "Pta"
full_sgi$ID[full_sgi$ID == "PtaDa"] <- "PtaDA"
full_sgi$ID[full_sgi$ID == "ptada"] <- "PtaDA"
full_sgi$ID[full_sgi$ID == "Pta flg22"] <- "PtaDA"

full_sgi$Plant[full_sgi$Plant == "colo"] <- "col"
full_sgi$Plant[full_sgi$Plant == "Col"] <- "col"
full_sgi$Plant[full_sgi$Plant == "Fls2efr"] <- "fls2efr"


#change 2609 to 2005
full_sgi$ID[full_sgi$ID == "2609"] <- "2005"
full_sgi$ID <- factor(full_sgi$ID)
#full_sgi$plate_exp <- factor(full_sgi$plate_exp)
full_sgi$Exp <- factor(full_sgi$Exp)

full_sgi$Plant <- factor(as.character(full_sgi$Plant))
levels(full_sgi$Plant)
colo_sgi <- full_sgi[full_sgi$Plant == "col",]

#Shows Experiment Bias - along with random plate bias
#ggplot(na.omit(colo_sgi[colo_sgi$ID == "mock",]), aes(Exp, gramsdtenthousand, color=plate_exp))+
#  geom_boxplot(outlier.shape = NA)
###
### I am going to shy away from this ###
###
#turn all the data to z-scores

#Create the histogram of the sd
hist_data <- colo_sgi[colo_sgi$ID == "mock",]
middle<-median(hist_data$gramsdtenthousand)
ggplot(hist_data, aes(x=gramsdtenthousand))+ 
  geom_histogram(binwidth = 1)+
  theme_bw()+
  scale_x_continuous(breaks=c(0, round(middle-(2*sd(colo_sgi$gramsdtenthousand[colo_sgi$ID == "mock"]))), 
         round( middle-(sd(hist_data$gramsdtenthousand))),
         round( middle-(sd(hist_data$gramsdtenthousand)*.5 )),
         middle,
         round( middle+(sd(hist_data$gramsdtenthousand)*.5 )),
         round( middle+(sd(hist_data$gramsdtenthousand))),
         round(middle+(2*sd(hist_data$gramsdtenthousand))),
         round(max(hist_data$gramsdtenthousand)))
       )
sd(hist_data$gramsdtenthousand)*.40711

#just a quick trial on loooking at all the medians and SD
hist_meds <- c()
hist_sds <- c()
for ( i in hist_data$plate_exp ) {
  hist_meds <- append(hist_meds, median(hist_data$gramsdtenthousand[hist_data$plate_exp==i]))
  hist_sds <- append(hist_sds, sd(hist_data$gramsdtenthousand[hist_data$plate_exp==i]))
}

hist(hist_meds) + abline(v=median(hist_data$gramsdtenthousand))
hist(hist_sds) + abline(v=sd(hist_data$gramsdtenthousand))

plot(hist_meds, hist_sds)

for ( i in unique(colo_sgi$plate_exp ) ) {
  colo_sgi$gramsdtenthousand[colo_sgi$plate_exp == i] <- (colo_sgi$gramsdtenthousand[colo_sgi$plate_exp == i] - median(colo_sgi$gramsdtenthousand[colo_sgi$plate_exp == i & colo_sgi$ID == "mock"], na.rm = T))/sd(colo_sgi$gramsdtenthousand[colo_sgi$plate_exp == i & colo_sgi$ID == "mock"], na.rm = T)
}
#CHECK TO SEE IF THE ADJUSTMENT WORKED
#ggplot(na.omit(colo_sgi[colo_sgi$ID == "mock",]), aes(Exp, gramsdtenthousand, color=plate_exp))+
#  geom_boxplot(outlier.shape = NA)
#Z-score lines up the mocks

#relevel
colo_sgi$ID <- relevel(colo_sgi$ID, ref = "mock")


#I AM ACTUALLY MODELING THE Z-SCORES PRODUCED. SHOULD I STILL KEEP THE PLATE AND EXPERIMENTAL EFFECT?
#Model I am using to call significance
less_rand <- lmer(gramsdtenthousand ~ ID + Exp + (1|plate_exp), data = colo_sgi)
summary(less_rand)
trial <- lm(gramsdtenthousand ~ ID + Exp, data = colo_sgi)
summary(trial)
#Compare the models with different data 
anova(less_rand, trial)
#The model with the random effect fits better

#Get the ROS data for visualization

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

order_names <- reorder(z_scores_per_peptide$V1,-z_scores_per_peptide$V2, median)

#plot the results
t <- summary(less_rand)
s <- as.data.frame( t$coefficients )
colnames(s)[2] <- "StdE"
s_inter <- s[1,]
s <- s[2:(nrow(s)-7),]
row.names(s) <- gsub("ID","", row.names(s))
#s$Estimate <- s$Estimate + s_inter$Estimate[1]

names_in_sgi <- levels(order_names)[levels(order_names) %in% row.names(s)]
s <- s[row.names(s)[match(names_in_sgi, row.names(s))],]

#I am putting the p-values through fdr correction
s$`Pr(>|t|)` <- p.adjust(s$`Pr(>|t|)`, method = "fdr")


s_high <- s[s$`Pr(>|t|)` < .001,]
s_med <- s[s$`Pr(>|t|)` < .01 & s$`Pr(>|t|)` > .001,]
s_low <- s[s$`Pr(>|t|)` < .05 & s$`Pr(>|t|)` > .01,]

#raw data with code
colo_sgi$ID <- as.character(colo_sgi$ID)
colo_sgi$ID[colo_sgi$ID == "flg22"] <- "Pta" 
colo_sgi$ID <- factor(colo_sgi$ID, levels=c(row.names(s),"mock"))
#mock_med <- median(na.omit(colo_sgi$gramsdtenthousand[colo_sgi$ID == "mock"]))
min_exps <- 100
for (i in unique(colo_sgi$ID)) {
  if ( nrow(na.omit(colo_sgi[colo_sgi$ID == i,])) < min_exps & nrow(na.omit(colo_sgi[colo_sgi$ID == i,])) != 0 ) {
    min_exps <- nrow(na.omit(colo_sgi[colo_sgi$ID == i,]))
  }
}
length(unique(colo_sgi$ID))
sum(as.character(unique(colo_sgi$ID)) %in% append(row.names(s_med), row.names(s_high)))

#I am going to change the outliers to reasonable values for plotting - Significance tests will not
#be affected because they have already been run.
colo_sgi$gramsdtenthousand[colo_sgi$gramsdtenthousand > 7] <- 7
colo_sgi$gramsdtenthousand[colo_sgi$gramsdtenthousand < -7] <- -7

colnames(s)[ncol(s)] <- "pval"
s$plot_pval <- s$pval
s$plot_pval[s$plot_pval>.05] <- .05

#Final Plot
alternate <- 0
#colo_sgi$plot_grams <- colo_sgi$gramsdtenthousand-median(colo_sgi$gramsdtenthousand[colo_sgi$ID=="mock"], na.rm = T)

#only plot nat variants from NV
nv_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/final_figure_antag/supplemental_data/supplemental_table2.csv")

p <- ggplot(na.omit(colo_sgi[colo_sgi$ID %in% c(as.character(nv_data$PeptideID), "mock"),]), aes(x=ID, y=gramsdtenthousand))+
  geom_point(color="grey", size=.25) +
  geom_errorbar(data=s[row.names(s) %in% c(as.character(nv_data$PeptideID), "mock"),], aes(x=row.names(s)[row.names(s) %in% c(as.character(nv_data$PeptideID), "mock")], y=Estimate, ymin=Estimate-StdE, ymax=Estimate+StdE, color=plot_pval), width=.75)+
  geom_errorbar(data=s[row.names(s) %in% c(as.character(nv_data$PeptideID), "mock"),], aes(x=row.names(s)[row.names(s) %in% c(as.character(nv_data$PeptideID), "mock")], y=Estimate, ymin=Estimate, ymax=Estimate, color=plot_pval))+
  scale_color_gradient(low = "firebrick2", high = "#333333")+
  stat_summary(fun.data = function(y){return(data.frame(y=sample(seq(4,4.5,.1), 1), label=paste(length(y), sep = ""), size=2))},
               geom = "text")+
  #geom_jitter(aes(fill=Plate))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(text = element_text(size=14))+
  #ylim(-9,9)+
  #geom_text(data = s_med, aes(row.names(s_med), y=2.5), label = "*", vjust=-1.5, inherit.aes = F, size=5)+
  #geom_text(data = s_high, aes(row.names(s_high), y=2.5), label = "*", vjust=-1.5,inherit.aes = F, size=5)+
  #geom_text(data = s_low, aes(row.names(s_low),y=12.5), label = "*", vjust =-1.5,inherit.aes = F, size=5)+
  ylab("Mock Adjusted Fresh Weight (z-score)") +
  xlab("Peptide ID") +
  ggtitle("SGI of Natural Variants")+
  scale_y_continuous( breaks=c(-4,-2,-1,-.5, 0, .5, 1,2,4))
#p <- p + geom_errorbar(data=s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate, ymin=Estimate-StdE, ymax=Estimate+StdE))
p

#Make a heatmap of significant sgi effects
important <- rbind(s_high, s_med)

mapping <- c()
for (i in levels(colo_sgi$ID)) {
  if ( i == "mock"){
    next
  }
  print(i)
  if ( i %in% row.names(important)) {
    if ( s$Estimate[row.names(s) == i] < 0 ){
      mapping <- append(mapping, c(i, s$pval[row.names(s) == i], "Significant", s$Estimate[row.names(s) == i]))
    }
    else{
      mapping <- append(mapping, c(i, .1, "NA", s$Estimate[row.names(s) == i]))
    }
  }
  else{
    mapping <- append(mapping, c(i, s$pval[row.names(s) == i], "NA", s$Estimate[row.names(s) == i]))
  }
}
signif_frame <- data.frame(matrix(mapping, ncol=4, byrow = T))
signif_frame$X1 <- factor(as.character(signif_frame$X1), levels = levels(colo_sgi$ID))
signif_frame$X3 <- as.character(signif_frame$X3)
signif_frame$X2 <- as.numeric(as.character(signif_frame$X2))
signif_frame$X2[signif_frame$X2>.1] <- .1
signif_frame$X4 <- as.numeric(as.character(signif_frame$X4))
signif_frame$X4[signif_frame$X4 > 0] <- 0

#219AE0
length(unique(signif_frame$X1[!(signif_frame$X1 %in% seq(5000,5007))]))

sgraph <- ggplot(signif_frame[signif_frame$X1 %in% as.character(nv_data$PeptideID),], aes(X1, "SGI"))+
  geom_tile(aes(fill=X4,color=X3, width=0.8, height=0.8), size=1)+
  scale_color_manual(values = c("White", "Black"))+
  scale_fill_gradient2(low="#662D91", mid="#219AE0", high = "#D9E021", midpoint = -1.6 )+
  scale_y_discrete(expand = c(0, 0))+
  scale_x_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.margin = unit(c(1,1,0,1), "cm"))+
  theme(text = element_text(size=14))+
  labs(x = "Peptide ID", y = "")+
  labs(fill="SGI effect")
 # theme(legend.position = "none")
sgraph

length(unique(z_scores_per_peptide$V1[(z_scores_per_peptide$V1 %in% as.character(nv_data$PeptideID))]))
length(signif_frame$X1[signif_frame$X1 %in% as.character(nv_data$PeptideID)])
#######
# I WILL BE TRYING TO PLOT THE SGI WITH THE ROS DATA #
#######
p1 <- ggplot(z_scores_per_peptide[(z_scores_per_peptide$V1 %in% as.character(nv_data$PeptideID)),], aes(x=reorder(V1,-V2, median), y=V2))+
  geom_boxplot( outlier.shape = NA)+
  ylab("Log10 AUC PtaDA Z-Score")+
 # xlab("Peptide ID")+
  labs(color="Plate Number")+
  ggtitle("Characterization of Natural Variants")+
  theme(text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 0),
        plot.title = element_text(hjust=.5),
        plot.margin = unit(c(1,1,0,1), "cm"))
  #geom_errorbar(data=s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate, ymin=Estimate-StdE, ymax=Estimate+StdE, color="Red"))

library(grid)
library(gridExtra)
#Plot the two on the same plot with height differences
gA <- ggplotGrob(p1)
gB <- ggplotGrob(sgraph)
gA$widths <- gB$widths
grid.newpage()
grid.draw(arrangeGrob(gA,gB, heights = c(1, .3)) )

#######################################
#I will be getting the ROS burst and SGI data for the Ralstonia peptides ##
p <- ggplot(na.omit(colo_sgi[colo_sgi$ID %in% c(1186, 5001:5007, 1410, "Pta", "Pa", "PtaDA", "Pa20", "mock"),]), aes(x=ID, y=gramsdtenthousand))+
  geom_point(color="grey", size=.25) +
  geom_errorbar(data=s[row.names(s) %in% c(1186, 5001:5007, 1410, "Pta", "Pa", "PtaDA", "Pa20", "mock"),], aes(x=row.names(s)[row.names(s) %in% c(1186, 5001:5007, 1410, "Pta", "Pa", "PtaDA", "Pa20", "mock")], y=Estimate, ymin=Estimate-StdE, ymax=Estimate+StdE, color=plot_pval), width=.75)+
  geom_errorbar(data=s[row.names(s) %in% c(1186, 5001:5007, 1410, "Pta", "Pa", "PtaDA", "Pa20", "mock"),], aes(x=row.names(s)[row.names(s) %in% c(1186, 5001:5007, 1410, "Pta", "Pa", "PtaDA", "Pa20", "mock")], y=Estimate, ymin=Estimate, ymax=Estimate, color=plot_pval))+
  scale_color_gradient(low = "firebrick2", high = "#333333")+
  stat_summary(fun.data = function(y){return(data.frame(y=sample(seq(4,4.5,.1), 1), label=paste(length(y), sep = ""), size=2))},
               geom = "text")+
  #geom_jitter(aes(fill=Plate))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(text = element_text(size=14))+
  #ylim(-9,9)+
  #geom_text(data = s_med, aes(row.names(s_med), y=2.5), label = "*", vjust=-1.5, inherit.aes = F, size=5)+
  #geom_text(data = s_high, aes(row.names(s_high), y=2.5), label = "*", vjust=-1.5,inherit.aes = F, size=5)+
  #geom_text(data = s_low, aes(row.names(s_low),y=12.5), label = "*", vjust =-1.5,inherit.aes = F, size=5)+
  ylab("Mock Adjusted Fresh Weight (z-score)") +
  xlab("Peptide ID") +
  ggtitle("SGI of Natural Variants")+
  scale_y_continuous( breaks=c(-4,-2,-1,-.5, 0, .5, 1,2,4))
p

sgraph <- ggplot(signif_frame[signif_frame$X1 %in% c(1186, 5001:5007, 1410, "Pta", "Pa", "PtaDA", "Pa20"),], aes(X1, "SGI"))+
  geom_tile(aes(fill=X4,color=X3, width=0.8, height=0.8), size=1)+
  scale_color_manual(values = c("White", "Black"))+
  scale_fill_gradient2(low="#662D91", mid="#219AE0", high = "#D9E021", midpoint = -1.6 )+
  scale_y_discrete(expand = c(0, 0))+
  scale_x_discrete(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.margin = unit(c(1,1,0,1), "cm"))+
  theme(text = element_text(size=14))+
  labs(x = "Peptide ID", y = "")+
  labs(fill="SGI effect")
# theme(legend.position = "none")
sgraph

length(unique(z_scores_per_peptide$V1[(z_scores_per_peptide$V1 %in% as.character(nv_data$PeptideID))]))
length(signif_frame$X1[signif_frame$X1 %in% as.character(nv_data$PeptideID)])

p1 <- ggplot(z_scores_per_peptide[(z_scores_per_peptide$V1 %in% c(1186, 5001:5007, 1410, "Pta", "Pa", "PtaDA", "Pa20")),], aes(x=reorder(V1,-V2, median), y=V2))+
  geom_boxplot( outlier.shape = NA)+
  ylab("Log10 AUC PtaDA Z-Score")+
  # xlab("Peptide ID")+
  labs(color="Plate Number")+
  ggtitle("Characterization of Natural Variants")+
  theme(text = element_text(size=14))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 0),
        plot.title = element_text(hjust=.5),
        plot.margin = unit(c(1,1,0,1), "cm"))
#geom_errorbar(data=s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate, ymin=Estimate-StdE, ymax=Estimate+StdE, color="Red"))


#Plot the two on the same plot with height differences
gA <- ggplotGrob(p1)
gB <- ggplotGrob(sgraph)
gA$widths <- gB$widths
grid.newpage()
grid.draw(arrangeGrob(gA,gB, heights = c(1, .3)) )

#############################
#####Check the fls2 data#####
#############################

unique(full_sgi$Plant)
fls2_sgi <- full_sgi[full_sgi$Plant == "fls2efr",]
unique(fls2_sgi$ID)
unique(colo_sgi$ID[!( colo_sgi$ID %in% fls2_sgi$ID)])

#turn values to z-scores
for ( i in unique(fls2_sgi$plate_exp ) ) {
  fls2_sgi$gramsdtenthousand[fls2_sgi$plate_exp  == i] <- (fls2_sgi$gramsdtenthousand[fls2_sgi$plate_exp  == i] - median(fls2_sgi$gramsdtenthousand[fls2_sgi$plate_exp  == i & fls2_sgi$ID == "mock"], na.rm = T))/sd(fls2_sgi$gramsdtenthousand[fls2_sgi$plate_exp  == i & fls2_sgi$ID == "mock"], na.rm = T)
}
ggplot(na.omit(fls2_sgi[fls2_sgi$ID == "mock",]), aes(Exp, gramsdtenthousand, color=plate_exp ))+
  geom_boxplot(outlier.shape = NA)
unique(colo_sgi$ID[!( colo_sgi$ID %in% fls2_sgi$ID)])

#model the same as col-0
less_rand <- lmer(gramsdtenthousand ~ ID + Exp + (1|plate_exp), data = fls2_sgi)
summary(less_rand)

#plot the results
t <- summary(less_rand)
s <- as.data.frame( t$coefficients )
colnames(s)[2] <- "StdE"
s_inter <- s[1,]
s <- s[2:(nrow(s)-1),]
row.names(s) <- gsub("ID","", row.names(s))
#s$Estimate <- s$Estimate + s_inter$Estimate[1]

names_in_sgi <- levels(order_names)[levels(order_names) %in% row.names(s)]
s <- s[row.names(s)[match(names_in_sgi, row.names(s))],]

#I am putting the p-values through fdr correction
s$`Pr(>|t|)` <- p.adjust(s$`Pr(>|t|)`, method = "fdr")

s_high <- s[s$`Pr(>|t|)` < .001,]
s_med <- s[s$`Pr(>|t|)` < .01 & s$`Pr(>|t|)` > .001,]
s_low <- s[s$`Pr(>|t|)` < .05 & s$`Pr(>|t|)` > .01,]

#raw data with code
fls2_sgi$ID <- as.character(fls2_sgi$ID)
fls2_sgi$ID[fls2_sgi$ID == "flg22"] <- "Pta" 
fls2_sgi$ID <- factor(fls2_sgi$ID, levels=c(row.names(s),"mock"))
#mock_med <- median(na.omit(colo_sgi$gramsdtenthousand[colo_sgi$ID == "mock"]))
min_exps <- 100
for (i in unique(fls2_sgi$ID)) {
  if ( nrow(na.omit(fls2_sgi[fls2_sgi$ID == i,])) < min_exps & nrow(na.omit(fls2_sgi[fls2_sgi$ID == i,])) != 0 ) {
    min_exps <- nrow(na.omit(fls2_sgi[fls2_sgi$ID == i,]))
  }
}
unique(colo_sgi$ID[!( colo_sgi$ID %in% fls2_sgi$ID)])
fls2_sgi[fls2_sgi$ID == "PtaDA",]


## Adjust plot
fls2_sgi$gramsdtenthousand[fls2_sgi$gramsdtenthousand > 7] <- 7
fls2_sgi$gramsdtenthousand[fls2_sgi$gramsdtenthousand < -7] <- -7

colnames(s)[ncol(s)] <- "pval"
s$plot_pval <- s$pval
s$plot_pval[s$plot_pval>.05] <- .05

#Final Plot
#colo_sgi$plot_grams <- colo_sgi$gramsdtenthousand-median(colo_sgi$gramsdtenthousand[colo_sgi$ID=="mock"], na.rm = T)
p <- ggplot(na.omit(fls2_sgi), aes(x=ID, y=gramsdtenthousand))+
  geom_point(color="grey", size=.5) +
  geom_errorbar(data=s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate, ymin=Estimate-StdE, ymax=Estimate+StdE, color=plot_pval), width=.75)+
  geom_errorbar(data=s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate, ymin=Estimate, ymax=Estimate, color=plot_pval))+
  scale_color_gradient(low = "firebrick2", high = "#333333")+
  stat_summary(fun.data = function(y){return(data.frame(y=-sample(seq(4,5,.1), 1), label=paste(length(y), sep = ""), size=4))},
             geom = "text")+
  #geom_jitter(aes(fill=Plate))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(text = element_text(size=14))+
  #ylim(-9,9)+
  #geom_text(data = s_med, aes(row.names(s_med), y=2.5), label = "*", vjust=-1.5, inherit.aes = F, size=5)+
  #geom_text(data = s_high, aes(row.names(s_high), y=2.5), label = "*", vjust=-1.5,inherit.aes = F, size=5)+
  #geom_text(data = s_low, aes(row.names(s_low),y=12.5), label = "*", vjust =-1.5,inherit.aes = F, size=5)+
  ylab("Mock Adjusted Fresh Weight (z-score)") +
  xlab("Peptide ID") +
  ggtitle("SGI of Natural Variants on fls2/efr")+
  scale_y_continuous( breaks=c(-4,-2,-1,-.5, 0, .5, 1,2,4))
#p <- p + geom_errorbar(data=s, aes(x=factor(row.names(s), levels=row.names(s)), y=Estimate, ymin=Estimate-StdE, ymax=Estimate+StdE))
p

