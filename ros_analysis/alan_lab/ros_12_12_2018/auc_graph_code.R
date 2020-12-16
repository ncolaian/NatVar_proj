#this script will get all the files and create dimensionless units to combine auc results
#THIS WAS MODIFIED FOR the fls2/efr mutant

library(ggplot2)

pos_control <- "100nm"
neg_control <- "1nm"

path_ros <- "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/alan_lab/ros_12_12_2018/"

#get file names along with location


file_list <- list.files(path = path_ros, pattern = "*auc_data.csv")
#num <- c()
#for ( i in 30:31) {
#  mm_files <- list.files(path=paste(path_ros, "8_", i, "_18", sep = ""), pattern = "*auc_data.csv")
#  file_list <- append(file_list, mm_files)
#  num <- append(num, rep(i, length(mm_files)))
#}

#read the files
csv_list <- list()
for ( i in 1:length(file_list)) {
  csv_list[[i]] <- read.csv(paste(path_ros, file_list[i], sep = ""), stringsAsFactors = F)
}
csv_list[[3]]

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

pta_scores <- inconstant_mat[inconstant_mat$Name == pos_control,]
pta_scores$factor <- max(as.numeric(pta_scores$Unitless_measure))/as.numeric(pta_scores$Unitless_measure)

scaled_score <- c()
for (i in 1:nrow(inconstant_mat) ) {
  scale_fact <- pta_scores$factor[pta_scores$Exp == inconstant_mat[i,2] ]# the column with the experiment is 2
  scaled_score <- append(scaled_score, (as.numeric(inconstant_mat$Unitless_measure[i])*scale_fact) )
}

#get rid of natural varian ids
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
csv_list[[3]]$AUC <- log10(csv_list[[3]]$AUC)
#csv_list[[3]]$plate <- as.character(csv_list[[3]]$plate)
p <- ggplot(csv_list[[3]], aes(x=Name, y=AUC))+
  geom_boxplot(outlier.shape = NA)+
  #geom_jitter( aes(color = plate))+
  ylab("Log10 AUC")+
  xlab("P. Tabaci Concentration")+
  labs(color="Plate Number")+
  theme(text = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

p <- p + stat_compare_means(method = "t.test", label = "p.signif", ref.group = "100nm", label.y = 8)
#p <- p+compare_means(AUC ~ Name, data=csv_list[[1]], method="anova")
p
