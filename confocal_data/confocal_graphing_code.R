# This script is an extension of z_score_sgi_fulldata.r
# it is located at: /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/sgi_data/z_score_sgi_fulldata.R

#Run that above file to the SGI graphing portion and then stop.
#I want 3 figures that correspond to the confocal data

confocal_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/confocal_data/confocal_data_formated.csv")

library(tidyr)

confocal_data <- gather(confocal_data, key="Mes_num", value = "Photon", M1:M12)

confocal_data <- na.omit(confocal_data)
confocal_data$log10_photons <- log10(confocal_data$Photon)

ggplot(confocal_data, aes(Peptide, log10(log10_photons)))+
  geom_boxplot()

#going to change flg22 to Pa
confocal_data$Peptide <- as.character(confocal_data$Peptide)
confocal_data$Peptide[confocal_data$Peptide == "flg22"] <- "Pa"
confocal_data$Peptide[confocal_data$Peptide == "flg20_n"] <- "flg20"
confocal_data$Peptide[confocal_data$Peptide == "flg22_n"] <- "Pa"

#going to go ahead and line up the data to the ROS burst
names_in_confocal <- c(levels(order_names)[levels(order_names) %in% confocal_data$Peptide], as.character(unique(confocal_data$Peptide)[!(unique(confocal_data$Peptide) %in% levels(order_names))])) 

confocal_data$Peptide <- factor(confocal_data$Peptide, levels = names_in_confocal)

ggplot(confocal_data, aes(Peptide, log10_photons))+
  geom_boxplot()

# I will add a one way anova with a tukey test
conf_aov <- aov(log10_photons ~ Peptide, confocal_data)
summary(conf_aov)
#perform the tukey test
library(agricolae)
confocal_posttukey <- HSD.test(conf_aov, trt = "Peptide")
confocal_posttukey <- HSD.test(conf_aov, trt = "Peptide")
TukeyHSD(conf_aov)

conf_res <- as.data.frame(confocal_posttukey$groups)
conf_res$pep <- row.names(confocal_posttukey$groups)

ggplot(confocal_data, aes(Peptide, log10_photons))+
  geom_boxplot()+
  geom_text(data=conf_res,aes(x=pep, label=groups, y=6.5))

# This is a bit cluttered. I am going to try and compare all of the peptides to flg20_n
library(ggpubr)
ggplot(confocal_data, aes(Peptide, log10_photons))+
  geom_boxplot()+
  stat_compare_means(ref.group = "flg20_n", method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     label = "p.signif")


#I color the data based on SGI since it is ordered by ROS
#get SGI
s$peptide <- rownames(s)

for ( i in as.character(unique(confocal_data$Peptide)) ) {
  if ( !(i %in% s$peptide) ) {
    next()
  }
  if ( s$Estimate[s$peptide == i] < 0 & s$pval[s$peptide == i] < 0.05 ) {
    confocal_data$sgi_val[confocal_data$Peptide == i] <- "Positive"
   # print(i)
  }
  else if (s$Estimate[s$peptide == i] > 0 | s$pval[s$peptide == i] > 0.05) {
    confocal_data$sgi_val[confocal_data$Peptide == i] <- "Negative"
   # print(i)
  }
  else {
    print(i)
  }
}

#confocal_data[confocal_data$Peptide == 94,]
s[s$peptide == 1877,]


ggplot(confocal_data, aes(Peptide, log10_photons))+
  geom_boxplot()+
  stat_compare_means(ref.group = "flg20", method = "t.test", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     label = "p.signif")+
  ylab("MYB51 Expression - Log10(photon)")

#Perform a Dunnets test to control for multiple testing against a single control

#####################
#### MYB51 Graph ####
#####################

set.seed(1994)
library(multcomp)
names_in_confocal <- c(levels(order_names)[levels(order_names) %in% confocal_data$Peptide], as.character(unique(confocal_data$Peptide)[!(unique(confocal_data$Peptide) %in% levels(order_names))])) 
confocal_data$Peptide <- factor(confocal_data$Peptide, levels = names_in_confocal)
confocal_data$Peptide <- relevel(confocal_data$Peptide, ref = "flg20")
conf_aov <- aov(log10_photons ~ Peptide, confocal_data)
t <- glht(conf_aov, linfct=mcp(Peptide="Dunnett"))
summary(t)
#signif_dunnet <- c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "", "*", "*", "","","*",
#                   "*","", "*", "", "", "*", "","", "*", "", "*","", "", "")

names_in_confocal <- c(levels(order_names)[levels(order_names) %in% confocal_data$Peptide], as.character(unique(confocal_data$Peptide)[!(unique(confocal_data$Peptide) %in% levels(order_names))])) 
confocal_data$Peptide <- factor(confocal_data$Peptide, levels = names_in_confocal)
names_t <- summary(t)
signif_dunnet <- c()
for ( i in names_t$test$pvalues ) {
  if (i < 0.01) {
    signif_dunnet <- append(signif_dunnet, "*")
  }
  else {
    signif_dunnet <- append(signif_dunnet, "")
  }
}
signif_dunnet <- append(signif_dunnet, "")

dunn_signif <- data.frame(levels(confocal_data$Peptide), signif_dunnet)

ggplot(confocal_data, aes(Peptide, log10_photons))+
  geom_boxplot(outlier.colour = NA)+
  geom_text(data = dunn_signif, aes(x=levels.confocal_data.Peptide., label=signif_dunnet, y=6.29), size=6)+
  stat_summary(fun.data = function(y){return(data.frame(y=6.34, label=paste(length(y), sep = ""), size=3))},
               geom = "text")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))+
  ylab("MYB51 Expression - Log10(photon)")


#This should be run after z_score_sgi_fulldata so that the bottom of the graph will have the same data as previous graphs
MYB51_signif_frame <- signif_frame[signif_frame$X1 %in% c(as.character(confocal_data$Peptide), "Pa20"),]


sgraph <- ggplot(MYB51_signif_frame, aes(X1, "SGI"))+
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
sgraph

######################
#### WRKY11 Graph ####
######################
wrky11_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/confocal_data/wrky11_confocal_graphing.csv")
colnames(wrky11_data) <- c("Peptide", "M1", "M2", "M3", "M4", "M5", "M6")
wrky11_data <- gather(wrky11_data, key="Mes_num", value = "Photon", M1:M6)
wrky11_data <- na.omit(wrky11_data)
wrky11_data$log10_photons <- log10(wrky11_data$Photon)
names_in_confocal <- c(levels(order_names)[levels(order_names) %in% wrky11_data$Peptide], as.character(unique(wrky11_data$Peptide)[!(unique(wrky11_data$Peptide) %in% levels(order_names))])) 
wrky11_data$Peptide <- factor(wrky11_data$Peptide, levels = names_in_confocal)

wrky11_data$Peptide <- relevel(wrky11_data$Peptide, ref = "flg20")
conf_aov <- aov(log10_photons ~ Peptide, wrky11_data)
t <- glht(conf_aov, linfct=mcp(Peptide="Dunnett"))
summary(t)

signif_dunnet <- c("*", "*", "*", "*", "*", "*", "", "", "*", "", "*", "","", "*","", "")
wrky11_data$Peptide <- factor(wrky11_data$Peptide, levels = names_in_confocal)

dunn_signif <- data.frame(levels(wrky11_data$Peptide), signif_dunnet)

ggplot(wrky11_data, aes(Peptide, log10_photons))+
  geom_boxplot(outlier.color = NA)+
  geom_text(data = dunn_signif, aes(x=levels.wrky11_data.Peptide., label=signif_dunnet, y=6), size=6)+
  stat_summary(fun.data = function(y){return(data.frame(y=6.17, label=paste(length(y), sep = ""), size=2.5))},
               geom = "text")+
  theme_classic2()+
  ylab("WRKY11 Expression - Log10(photon)")

WRKY11_signif_frame <- signif_frame[signif_frame$X1 %in% c(as.character(wrky11_data$Peptide), "Pa20"),]
sgraph <- ggplot(WRKY11_signif_frame, aes(X1, "SGI"))+
  geom_tile(aes(fill=X4,color=X3, width=0.8, height=0.8), size=1)+
  scale_color_manual(values = c("White", "Black"))+
  scale_fill_gradient2(low="#662D91", mid="#219AE0", high = "#D9E021", midpoint = -1.6 )+
  scale_y_discrete(expand = c(0, 0))+
  scale_x_discrete(expand = c(0,0))+
  theme_void()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),
        plot.margin = unit(c(1,1,0,1), "cm"))+
  theme(text = element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x = "Peptide ID", y = "")+
  labs(fill="SGI effect")
sgraph
