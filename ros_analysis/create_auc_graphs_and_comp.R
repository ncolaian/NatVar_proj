#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

library(getopt)
library(tidyr)
library(ggplot2)
library(Bolstad2)
library(ggpubr)

#try to reproduce clarks code in R

params = matrix(c(
  "input_file", "i", 1, "character",
  "out_dir", "o", 1, "character",
  "ref_name", "r",1,"character",
  "base_out_name", "b",1,"character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

##############
#### TEST #### 
##############
#opt$input_file = "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/reelic_1391/135921.txt"
#opt$ref_name = "mock"

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

############################

file <- read.csv(opt$input_file, stringsAsFactors = F)

#first need to go from wide to long
if ( sum(is.na(file$A9)) == 0 ) {
  file <- gather(file, key = "replicate", value = "measurement", A1:A12)
}
if ( sum(is.na(file$A9)) != 0 ) {
  file <- gather(file, key = "replicate", value = "measurement", A1:A8)
}
if ( sum(is.na(file$A6)) != 0 ) {
  file <- gather(file, key = "replicate", value = "measurement", A1:A5)
}
file <- na.omit(file)

#file <- file[!(file$replicate == "A1"),]
sum_file <- summarySE(file, measurevar = "measurement", groupvars = c("Name","Time"))

pd <- position_dodge()
plot <- ggplot(sum_file, aes(x=Time/60, y=measurement, color=Name))+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se))+
    ylab("Photocount")+
  xlab("Seconds")+
  theme(text = element_text(size=20))+
  labs(color='Peptide') 

#integrate the area under the curve
mat_mm <- c()
for ( n in unique(file$Name) ) {
  for ( r in unique(file$replicate) ) {
    integral_calc <- sintegral(file$Time[file$Name == n & file$replicate == r], file$measurement[file$Name == n & file$replicate == r])$int
    mat_mm <- append(mat_mm, c(n,r,integral_calc))
  }
}
int_mat <- matrix(mat_mm, ncol = 3, byrow = T)
int_mat <- data.frame(int_mat, stringsAsFactors = F)
colnames(int_mat) <- c("Name", "Replicate", "AUC")
int_mat$AUC <- as.numeric(int_mat$AUC)
box_auc <- ggplot(int_mat, aes(x=Name, y=AUC))+
  geom_point()+
  geom_boxplot()+
  theme(text = element_text(size=20))

box_auc <- box_auc + stat_compare_means(method = "t.test",label = "p.signif", ref.group = opt$ref_name, label.y = max(int_mat$AUC)+sd(int_mat$AUC))
box_auc <- box_auc+xlab("Peptide")
ggsave(paste(opt$out_dir, opt$base_out_name, ".auc_means.pdf", sep=""),plot=box_auc)
ggsave(paste(opt$out_dir, opt$base_out_name, ".averaged_signal.pdf", sep=""),plot=plot)
write.csv(sum_file,paste(opt$out_dir, opt$base_out_name,".summary_info.csv", sep = ""), row.names = F, quote = F)
write.csv(int_mat, paste(opt$out_dir, opt$base_out_name, ".auc_data.csv", sep = ""), row.names = F, quote = F)
          