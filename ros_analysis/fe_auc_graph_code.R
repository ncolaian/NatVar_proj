#this script will get all the files and create dimensionless units to combine auc results
#THIS WAS MODIFIED FOR the fls2/efr mutant

library(ggplot2)
library(ggpubr)

pos_control <- "cPta"
neg_control <- "fePta"

path_ros <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/fls2_efr/"
second_ros <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/fls2efr_10_16_19/"

#get file names along with location


file_list <- list.files(path = path_ros, pattern = "*auc_data.csv")

frst_path <- length(file_list)

file_list <- append(file_list, list.files(path=second_ros, pattern = "*auc_data.csv"))

#read the files
csv_list <- list()
for ( i in 1:(length(file_list))) {
  if ( i <= frst_path ) {
    csv_list[[i]] <- read.csv(paste(path_ros, file_list[i], sep = ""), stringsAsFactors = F)
    #csv_list[[i]]$AUC <- log10(csv_list[[i]]$AUC)-min(log10(csv_list[[i]]$AUC))
  }
  else {
    csv_list[[i]] <- read.csv(paste(second_ros, file_list[i], sep = ""), stringsAsFactors = F)
    #csv_list[[i]]$AUC <- log10(csv_list[[i]]$AUC-min(csv_list[[i]]$AUC)+1)
  }
  csv_list[[i]]$plate = i
  if ( i == 1 ) {
    next
  }
  csv_list[[1]] <- rbind(csv_list[[1]], csv_list[[i]])
}

#remove the things i don't want in my figure
csv_list[[1]] <- csv_list[[1]][csv_list[[1]]$Name != "dud" & csv_list[[1]]$Name != "water" & csv_list[[1]]$Name != "cPpnlp20" & csv_list[[1]]$Name != "efPpnlp20" &
                                 csv_list[[1]]$Name != "cMpnlp20" & csv_list[[1]]$Name != "efMpnlp20" & csv_list[[1]]$Name != "DW" & csv_list[[1]]$Name != "feDW", ]
unique(csv_list[[1]]$Name)
table(csv_list[[1]]$Name)
#add a column for efr vs col-o
for ( i in 1:nrow(csv_list[[1]]) ) {
  if ( grepl("c", csv_list[[1]]$Name[i]) ) {
    csv_list[[1]]$plant[i] <- "Col-O"
    csv_list[[1]]$Name[i] <- strsplit(csv_list[[1]]$Name[i], "c")[[1]][2]
  }
  else{
    csv_list[[1]]$plant[i] <- "Fls2/Efr"
    if ( grepl("ef", csv_list[[1]]$Name[i]) ) {
      csv_list[[1]]$Name[i] <- strsplit(csv_list[[1]]$Name[i], "ef")[[1]][2]
    }
    else{
      csv_list[[1]]$Name[i] <- strsplit(csv_list[[1]]$Name[i], "fe")[[1]][2]
    }
  }
}

csv_list[[1]]$Name[csv_list[[1]]$Name == 1493] <- 1943

csv_list[[1]]$Name[csv_list[[1]]$Name == "pa"] <- "Pa"
csv_list[[1]]$AUC <- log10(csv_list[[1]]$AUC)
csv_list[[1]][csv_list[[1]]$Name == "Pa",]
#On each respective plate i want to get rid of systematic bias
for ( i in unique(csv_list[[1]]$plate) ) {
    csv_list[[1]]$AUC[csv_list[[1]]$plate == i] <- (csv_list[[1]]$AUC[csv_list[[1]]$plate == i] - median(csv_list[[1]]$AUC[csv_list[[1]]$plate == i & csv_list[[1]]$plant == "Fls2/Efr"]))
    #csv_list[[1]]$AUC[csv_list[[1]]$plate == i & csv_list[[1]]$Name == j] <- log10(csv_list[[1]]$AUC[csv_list[[1]]$plate == i & csv_list[[1]]$Name == j]+1)
    #csv_list[[1]]$AUC[csv_list[[1]]$plate == i & csv_list[[1]]$Name == j] <- (csv_list[[1]]$AUC[csv_list[[1]]$plate == i & csv_list[[1]]$Name == j] - min(csv_list[[1]]$AUC[csv_list[[1]]$plate == i& csv_list[[1]]$Name == j ]))

}

#csv_list[[1]]$AUC <- log10(csv_list[[1]]$AUC+1)
#for ( i in unique(csv_list[[1]]$plate) ) {
#  for ( j in unique(csv_list[[1]]$Name[csv_list[[1]]$plate == i]) ) {
#    csv_list[[1]]$AUC[csv_list[[1]]$plate == i & csv_list[[1]]$Name == j] <- ((csv_list[[1]]$AUC[csv_list[[1]]$plate == i & csv_list[[1]]$Name == j] - median(csv_list[[1]]$AUC[csv_list[[1]]$plate == i& csv_list[[1]]$Name == j &csv_list[[1]]$plant=="Fls2/Efr" ]))/sd(csv_list[[1]]$AUC[csv_list[[1]]$plate == i& csv_list[[1]]$Name == j &csv_list[[1]]$plant=="Fls2/Efr" ]))
#  }
#}
#I need to order the peptides based on the ros data
#I am going to read in the plot data and then order the ID's based on that
ros_levels <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/882_data/updated_work_10_23_18/new_clusters/all_nat_variants_ptanorm_wheel_11_12_19.csv")
csv_list[[1]]$Name[!(csv_list[[1]]$Name %in% ros_levels$V2)]
csv_list[[1]]$Name[!(csv_list[[1]]$Name %in% ros_levels$V2)]
ros_levels$V2[(!ros_levels$V2[ros_levels$V3 > 0] %in% csv_list[[1]]$Name)]
csv_list[[1]]$Name[csv_list[[1]]$Name == "498"] <- "468"
csv_list[[1]] <- csv_list[[1]][csv_list[[1]]$Name != "1787",]
ros_levels$V2[ros_levels$V3 > 0 & !(ros_levels$V2 %in% csv_list[[1]]$Name)]
csv_list[[1]]$Name <- factor(csv_list[[1]]$Name, levels=unique(ros_levels$V2[order(ros_levels$V3, decreasing = T)]))

#Get rid of some obvious problems with plates etc. 
csv_list[[1]] <- csv_list[[1]][ !(csv_list[[1]]$Name == "Pa" & csv_list[[1]]$AUC > 2 & csv_list[[1]]$plant == "Fls2/Efr"), ]
csv_list[[1]] <- csv_list[[1]][!(csv_list[[1]]$Name == "1544" & csv_list[[1]]$plate %in% c(2,3) ),]
table(csv_list[[1]]$Name)
csv_list[[1]] <- csv_list[[1]][csv_list[[1]]$AUC > -10,]
csv_list[[1]] <- csv_list[[1]][!(csv_list[[1]]$Name %in% seq(5000,5008)),]
library(plyr)
trial <- csv_list[[1]] %>% ddply(.(Name,plant), "summarise",
                        N    = sum(!is.na(AUC)),
                        mean = mean(AUC, na.rm=TRUE),
                        med = median(AUC, na.rm=TRUE),
                        sd   = sd(AUC, na.rm=TRUE),
                        se   = sd / sqrt(N) )

csv_list[[1]]$plate <- factor(csv_list[[1]]$plate)
ggplot(csv_list[[1]], aes(Name, AUC, color=plant))+
  geom_errorbar(inherit.aes = F, data = trial, aes(x=Name,y=med,color=plant, ymin=med-se, ymax=med+se))+
  geom_point()+
  stat_compare_means(aes(group=plant), method = "t.test", label = "p.signif", size=9, symnum.args= list(cutpoints = c(0, 0.05, 1), symbols = c("*", "")))+
  xlab("Peptide ID")+
  labs(color="Genotype")+
  ylab("Plate Adjusted Log10(AUC)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#######################################################################
#### DEPRICATED ####
####################

# go through and get all the data for each number
inconstant_mat <- matrix(ncol = 4, nrow = 0)

for ( i in 1:length(csv_list) ) {
  pos_mean <- mean(csv_list[[i]]$AUC[csv_list[[i]]$Name == pos_control])
  neg_mean <- mean(csv_list[[i]]$AUC[csv_list[[i]]$Name == neg_control])
  pos_sd <- sd(csv_list[[i]]$AUC[csv_list[[i]]$Name == pos_control])
  for ( j in unique(csv_list[[i]]$Name) ) {
    #I am scaling every value to the variance of the control - this allows me to compare values across samples because
    # The values are unitless
    j_mean <- mean(csv_list[[i]]$AUC[csv_list[[i]]$Name == j])
    val <- (j_mean - neg_mean)/pos_sd
    
    # I am looking to see where the values are
    if ( j_mean >= pos_mean) {
      position <- 1
    }
    else if ( j_mean <= neg_mean ) {
      position <- 2
    }
    else if ( j_mean >= pos_mean-(sqrt(pos_sd)) ) {
      position <- 3
    }
    else {
      position <- 4
    }
    
    inconstant_mat <- rbind(inconstant_mat, matrix(c(j,i,val,position), ncol = 4))
  }
  csv_list[[i]]$plate <- i
  if ( i > 1 ) {
    csv_list[[1]] <- rbind(csv_list[[1]], csv_list[[i]])
  }
}
inconstant_mat <- as.data.frame(inconstant_mat, stringsAsFactors = F)
colnames(inconstant_mat) <- c("Name", "Exp", "Unitless_measure", "Part")

portion_count <- matrix(ncol = 5, nrow = 0)
for ( i in unique(inconstant_mat$Name) ) {
  trimws(i)
  print(i)
  add <- c(i, 0, 0, 0, 0)
  for ( part in inconstant_mat$Part[inconstant_mat$Name == i] ) {
    add[as.numeric(part)+1] <- as.numeric(add[as.numeric(part)+1])+1
  }
  portion_count <- rbind(portion_count, matrix( add, ncol = 5))
}
colnames(portion_count) <- c("Peptide","Greater_pos", "Less_neg", "1sd_pos", "Else")
write.csv(portion_count,file = "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/27-31_analysis.csv" , quote = F, row.names = F)

pta_scores <- inconstant_mat[inconstant_mat$Name == pos_control,]
pta_scores$factor <- max(as.numeric(pta_scores$Unitless_measure))/as.numeric(pta_scores$Unitless_measure)

scaled_score <- c()
for (i in 1:nrow(inconstant_mat) ) {
  scale_fact <- pta_scores$factor[pta_scores$Exp == inconstant_mat[i,2] ]# the column with the experiment is 2
  scaled_score <- append(scaled_score, (as.numeric(inconstant_mat$Unitless_measure[i])*scale_fact) )
}

#get rid of natural varian ids
ids_to_exclude <- c("161","108","359")
inconstant_mat$scaled_score <- scaled_score
inconstant_mat <- inconstant_mat[!(inconstant_mat$Name %in% ids_to_exclude),]
#get rid of bad plate
plot_scaled_scores <- ggplot(inconstant_mat, aes(x=Name, y=scaled_score))+
  geom_boxplot()+
  xlab("UniqueID")+
  ylab("Pta Scaled AUC Score")
plot_scaled_scores

get_seq <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/youssef_data/finalized_useq_8_22_2018.csv", stringsAsFactors = F)

wheel_data <- c()
for ( i in unique(inconstant_mat$Name) ) {
  score <- mean(inconstant_mat$scaled_score[inconstant_mat$Name == i])
  if ( i %in% as.character( get_seq$UseqId ) ) {
    stri <- get_seq$Peptides[get_seq$UseqId == as.numeric(i)]
  }
  else {
    stri <- "A"
  }
  print(stri)
  wheel_data <- append(wheel_data, c(stri, i, score))
}
wheel_data <- matrix(wheel_data, ncol = 3, byrow = T)

write.csv(wheel_data,file = "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/27-31_analysis_wheel.csv" , quote = F, row.names = F)


#CSV LIST 1 CONTAINS ALL THE AUC SCORES IN ONE
#This allows me to create a nice graph
csv_list[[1]]$AUC <- log10(csv_list[[1]]$AUC)
csv_list[[1]]$plate <- as.character(csv_list[[1]]$plate)
ggplot(csv_list[[1]], aes(x=Name, y=AUC))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter( aes(color = plate))+
  ylab("Log10 AUC")+
  xlab("Peptide ID")+
  labs(color="Plate Number")
