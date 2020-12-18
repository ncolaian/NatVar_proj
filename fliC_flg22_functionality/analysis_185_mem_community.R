#Lets make a figure

#I am trying to understand the relationship of flg22 and defense in the context of a bacterial community

genome_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/882_data/updated_work_10_23_18/new_clusters/t6ss_comb_meta/genome_metadata.tsv")
isai_data <- readRDS("/Users/nicholascolaianni/Documents/dangl_lab/Isai_data/stress_experiment_data/dat_amplicon_useq97_4stresses.RDS")
#pull in the sequence to id data 
map_seq_name <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Omri_Data/no_dup_clusters.csv")

sgi_ros_data <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/final_figure_antag/supplemental_data/supplemental_table2.csv")


isai_data$RelativeAbundance$Tab
map_seq_name$Useq

#the genomes that have antagonist peptides
ids_antag <- c(2009,2004,1410, 1949, 2003, 1186, 5009, 5013, 5014, 5015)
gene_data <- read.delim("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/882_data/updated_work_10_23_18/new_clusters/t6ss_comb_meta/gene_metadata_file.tsv")
gene_data$full_seq <- sapply(1:nrow(gene_data), function(x){
  paste(lapply(gene_data[x,13:34],as.character),collapse = "" )
})
antag_seqs <- c("QRLSTGLRVEQAADNAAYWSIA", "DRVSSGEKVGKAADNVAYWSIS", "QRLSTGLRVNSAQDDAAAYAVA", 
                "NRISTGLKIGEAKDNAAYWAIS", "NRVSSGYKVSQASDNVAYWSIS", "QRLSTGLRVNSAKDDAAAYAVA",
                "SRINTGLKVASTKDDSASYTIA", "TRLSSGLKINSAKDDAAGMQIA", "ARLSSGARIGGASDDPAGQALA",
                "ARLSSGTSIASAADDPAGQAIT")
genomes_with_antag <- gene_data$taxon_oid[gene_data$full_seq %in% antag_seqs]

genome_data$antag[genome_data$taxon_oid %in% genomes_with_antag] <- 1
genome_data[genome_data$taxon_oid %in% genomes_with_antag,]
root_rel_abund[root_rel_abund$taxon_oid %in% genomes_with_antag,]

isai_data$Rarefied$Map[which(isai_data$Rarefied$Map$condition %in% c("1000Pi", "5.5pH", "21C") & isai_data$Rarefied$Map$Tissue == "Root"),]

root_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$Tissue == "Root")]
shoot_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$Tissue == "Shoot")]
agar_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$Tissue == "Agar")]

isai_data$Rarefied$Map[which(isai_data$Rarefied$Map$ID_Matrix == "P1SalinityC5"),]

###################
#### FliC DATA ####
###################

#This function will get the relative abundance of strains with each type of FliC protein
get_fliC_percentages <- function(rel_abundance_linked_to_bact_information, rel_abund_start = 1, rel_abund_end, site, condition) {
  rt_clust_prop <- c()
  rt_clust_prop <- append(rt_clust_prop,(colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end] * rel_abundance_linked_to_bact_information$clust1, na.rm = T)-colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end] * rel_abundance_linked_to_bact_information$clust1[duplicated(rel_abundance_linked_to_bact_information$Useq)], na.rm = T))/
                            (colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end]) - colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end])) )
  
  rt_clust_prop <- append(rt_clust_prop,(colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end] * rel_abundance_linked_to_bact_information$clust2, na.rm = T)-colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end] * rel_abundance_linked_to_bact_information$clust2[duplicated(rel_abundance_linked_to_bact_information$Useq)], na.rm = T))/
                            (colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end]) - colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end])) )
  
  rt_clust_prop <- append(rt_clust_prop,(colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end] * rel_abundance_linked_to_bact_information$clust3, na.rm = T)-colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end] * rel_abundance_linked_to_bact_information$clust3[duplicated(rel_abundance_linked_to_bact_information$Useq)], na.rm = T))/
                            (colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end]) - colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end])) )
  
  rt_clust_prop <- append(rt_clust_prop,(colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end] * rel_abundance_linked_to_bact_information$antag, na.rm = T)-colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end] * rel_abundance_linked_to_bact_information$antag[duplicated(rel_abundance_linked_to_bact_information$Useq)], na.rm = T))/
                            (colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end]) - colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end])) )
  
  rt_clust_prop <- append(rt_clust_prop,(colSums(rel_abundance_linked_to_bact_information[is.na(rel_abundance_linked_to_bact_information$NumfliC),rel_abund_start:rel_abund_end], na.rm = T)-colSums(rel_abundance_linked_to_bact_information[is.na(rel_abundance_linked_to_bact_information$NumfliC),rel_abund_start:rel_abund_end][duplicated(rel_abundance_linked_to_bact_information$Useq[is.na(rel_abundance_linked_to_bact_information$NumfliC)]),rel_abund_start:rel_abund_end], na.rm = T))/
                            (colSums(rel_abundance_linked_to_bact_information[,rel_abund_start:rel_abund_end]) - colSums(rel_abundance_linked_to_bact_information[duplicated(rel_abundance_linked_to_bact_information$Useq),rel_abund_start:rel_abund_end])) )
  
  rt_clust_prop <- matrix(rt_clust_prop, ncol = 1, byrow = F)
  rt_clust_prop <- as.data.frame(rt_clust_prop)
  rt_clust_prop$clust <- c(rep("1", rel_abund_end-rel_abund_start+1),rep("2", rel_abund_end-rel_abund_start+1), rep("3", rel_abund_end-rel_abund_start+1),rep("A", rel_abund_end-rel_abund_start+1), rep("No_fliC", rel_abund_end-rel_abund_start+1) )
  rt_clust_prop$place <- site
  rt_clust_prop$cond <- condition
  return(rt_clust_prop)
}

library(dplyr)
root_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$Tissue == "Root")]
root_rel_abund <- as.data.frame(root_rel_abund) 
num_condit <- ncol(root_rel_abund)
root_rel_abund$Useq <- row.names(root_rel_abund)
root_rel_abund <- left_join(root_rel_abund, map_seq_name, by="Useq" )
#root
root_rel_abund <- left_join(root_rel_abund, genome_data, by="taxon_oid")
unique(root_rel_abund$Useq)
#colSums(root_rel_abund[is.na(root_rel_abund$NumfliC),][duplicated(root_rel_abund$Useq[is.na(root_rel_abund$NumfliC)]), c(1,2)], na.rm = T)
#Go through conditions that are control and those that were detrimental to the plant
#root temperature
root <- matrix(nrow = 0, ncol = 4)
for ( i in unique(isai_data$Rarefied$Map$condition ) ) {
  if(i == 'Inoculum'){
    next
  }
  root_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$condition == i & isai_data$RelativeAbundance$Map$Tissue == "Root")]
  root_rel_abund <- as.data.frame(root_rel_abund) 
  num_condit <- ncol(root_rel_abund)
  root_rel_abund$Useq <- row.names(root_rel_abund)
  root_rel_abund <- left_join(root_rel_abund, map_seq_name, by="Useq" )
  root_rel_abund <- left_join(root_rel_abund, genome_data, by="taxon_oid")
  root <- rbind(root, get_fliC_percentages(root_rel_abund, rel_abund_end = num_condit, site = "root", condition=i) )
}

root$cond <- factor(root$cond, levels = unique(root$cond))

ggplot(root, aes(place, V1, color=clust))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = .2))+
  facet_wrap(facets = "cond")

#shoot
shoot <- matrix(nrow = 0, ncol = 4)
for ( i in unique(isai_data$Rarefied$Map$condition ) ) {
  if(i == 'Inoculum'){
    next
  }
  shoot_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$condition == i & isai_data$RelativeAbundance$Map$Tissue == "Shoot")]
  shoot_rel_abund <- as.data.frame(shoot_rel_abund) 
  num_condit <- ncol(shoot_rel_abund)
  shoot_rel_abund$Useq <- row.names(shoot_rel_abund)
  shoot_rel_abund <- left_join(shoot_rel_abund, map_seq_name, by="Useq" )
  shoot_rel_abund <- left_join(shoot_rel_abund, genome_data, by="taxon_oid")
  shoot <- rbind(shoot, get_fliC_percentages(shoot_rel_abund, rel_abund_end = num_condit, site = "shoot", condition=i) )
}

shoot$cond <- factor(shoot$cond, levels = unique(shoot$cond))

ggplot(shoot, aes(place, V1, color=clust))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = .2))+
  facet_wrap(facets = "cond")

#agar
agar <- matrix(nrow = 0, ncol = 4)
for ( i in unique(isai_data$Rarefied$Map$condition ) ) {
  if(i == 'Inoculum'){
    next
  }
  agar_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$condition == i & isai_data$RelativeAbundance$Map$Tissue == "Agar")]
  agar_rel_abund <- as.data.frame(agar_rel_abund) 
  num_condit <- ncol(agar_rel_abund)
  agar_rel_abund$Useq <- row.names(agar_rel_abund)
  agar_rel_abund <- left_join(agar_rel_abund, map_seq_name, by="Useq" )
  agar_rel_abund <- left_join(agar_rel_abund, genome_data, by="taxon_oid")
  agar <- rbind(agar, get_fliC_percentages(agar_rel_abund, rel_abund_end = num_condit, site = "agar", condition=i) )
}
agar$cond <- factor(agar$cond, levels = unique(agar$cond))

#First I am going to plot the control conditions 
library(ggpubr)
control_conditions <- c("1000Pi", "21C")
fliC_percent <- rbind(root, shoot, agar)
fliC_percent <- fliC_percent[fliC_percent$cond %in% control_conditions,]

#antagonists are not a clade
fliC_percent <- fliC_percent[fliC_percent$clust != "A",]

#get the data in a format to perform a tukey test
fliC_percent$tukey_groups <- paste(fliC_percent$clust, fliC_percent$place, sep = '_') 

total_anova <- aov(V1 ~ tukey_groups * cond, data=fliC_percent )
summary(total_anova)
library(agricolae)

tukey_res <- (HSD.test(total_anova, trt = "tukey_groups"))


count_sum <- fliC_percent %>% 
  group_by(clust, place) %>%
  dplyr::summarise( count_d = n(), 
                    total_abund = sum(V1))

ggplot(fliC_percent, aes(place, V1, shape=clust))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitterdodge(jitter.width = .2), size=3)+
  #facet_wrap(facets = "cond")+
  ylab("Relative Abundance of Colonized Bacteria With Each FliC")+
  xlab("Fraction")+
  labs(shape="Clade")+
  ggtitle("FliC type across controls")+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5))

summary <- fliC_percent %>%
  group_by(clust, place, cond) %>%
  dplyr::summarise(count=n(),
                   med = median(V1),
                   sd_mes = sd(V1))
print(summary)

#get the fliC data for the sick conditions
sick_fliC_percent <- rbind(root, shoot, agar)
sick_conditions <- c("200NaCl" )
sick_fliC_percent <- sick_fliC_percent[sick_fliC_percent$cond %in% sick_conditions,]

#antagonists are not a clade
sick_fliC_percent <- sick_fliC_percent[sick_fliC_percent$clust != "A",]

#get the data in a format to perform a tukey test
sick_fliC_percent$tukey_groups <- paste(sick_fliC_percent$clust, sick_fliC_percent$place, sep = '_') 

total_anova <- aov(V1 ~ tukey_groups, data=sick_fliC_percent )
summary(total_anova)

tukey_res <- (HSD.test(total_anova, trt = "tukey_groups"))

count_sum <- sick_fliC_percent %>% 
  group_by(clust, place) %>%
  dplyr::summarise( count_d = n())

pdf(file='/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/185_trial/final_figs/flic_type_high_salt.pdf', height = 5, width = 8)
ggplot(sick_fliC_percent, aes(place, V1, shape=clust))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitterdodge(jitter.width = .2), size=2.3)+
  ylab("Relative Abundance of Bacteria With Each FliC")+
  xlab("Fraction")+
  labs(shape="Clade")+
  ggtitle("Community FliC Profiles in 200 mM NaCl")+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5))

dev.off()


root_rel_abund[root_rel_abund$full_seq == "TRLSSGLKINSAKDDAAGLQIA",]

########################
#### Flg22 analysis ####
########################
#### REAL FLG22 ANALYSIS ####

get_flg22_percentages <- function(rel_abundance_linked_to_bact_information, rel_abund_start = 1, rel_abund_end, site, condition, col_clust_start, useq_col, flg22_seq_col=2) {
  library(matrixStats)
  library(dplyr)
  library(tidyr)
  flg22_percentages <- matrix(ncol=4, nrow = 0)
  colnames(flg22_percentages) <- c("full_seq", "clust1", "clust2", 'clust3')
  for ( i in rel_abund_start:rel_abund_end ) {
    col_name <- colnames(rel_abundance_linked_to_bact_information)[i+2]
    #2 is the full flg22 seqeunce
    variants_and_rel_abund <- na.omit(rel_abundance_linked_to_bact_information[,c(flg22_seq_col,i+2,col_clust_start,col_clust_start+1,col_clust_start+2)][rel_abundance_linked_to_bact_information$NumfliC > 0 & !duplicated(rel_abundance_linked_to_bact_information[c(flg22_seq_col,useq_col)]),] %>%
                                        dplyr::group_by(full_seq, clust1, clust2, clust3) %>%
                                        dplyr::summarise(
                                          sum_vals = sum(!!rlang::sym(col_name))
                                        ))
    variants_and_rel_abund <- as.data.frame(variants_and_rel_abund)
    colnames(variants_and_rel_abund)[5] <- col_name
    flg22_percentages <- merge(flg22_percentages, variants_and_rel_abund, by=c("full_seq", "clust1",   "clust2",   "clust3"), all = T, no.dups = T)
  }
  flg22_percentages <- pivot_longer(flg22_percentages, 5:ncol(flg22_percentages), names_to = "exp_id", values_to = "flg22_abundance")
  
  
  #flg22_percentages$median <- rowMedians(as.matrix(flg22_percentages[5:ncol(flg22_percentages)]), na.rm = T)
  #sum(flg22_percentages$med, na.rm = T)
  
  flg22_percentages$cond <- condition
  flg22_percentages$place <- site
  # print(flg22_percentages)
  return(flg22_percentages)
}


root_flg22 <- matrix(nrow = 0, ncol = 7)
for ( i in unique(isai_data$Rarefied$Map$condition ) ) {
  if(i == 'Inoculum'){
    next
  }
  root_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$condition == i & isai_data$RelativeAbundance$Map$Tissue == "Root")]
  root_rel_abund <- as.data.frame(root_rel_abund) 
  num_condit <- ncol(root_rel_abund)
  root_rel_abund$Useq <- row.names(root_rel_abund)
  root_rel_abund <- left_join(root_rel_abund, map_seq_name, by="Useq" )
  root_rel_abund <- left_join(root_rel_abund, genome_data, by="taxon_oid")
  root_rel_abund <- left_join(gene_data[,c(1,ncol(gene_data))],root_rel_abund, by="taxon_oid" )
  root_flg22 <- rbind(root_flg22, get_flg22_percentages(root_rel_abund, rel_abund_end = num_condit, site = "root", condition=i, col_clust_start = which(colnames(root_rel_abund) == "clust1"), useq_col = which(colnames(root_rel_abund) == "Useq")) )
}

root_flg22$cond <- factor(root_flg22$cond, levels = unique(root$cond))

#get the shoot flg22 seq data
shoot_flg22 <- matrix(nrow = 0, ncol = 7)
for ( i in unique(isai_data$Rarefied$Map$condition ) ) {
  if(i == 'Inoculum'){
    next
  }
  shoot_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$condition == i & isai_data$RelativeAbundance$Map$Tissue == "Shoot")]
  shoot_rel_abund <- as.data.frame(shoot_rel_abund) 
  num_condit <- ncol(shoot_rel_abund)
  shoot_rel_abund$Useq <- row.names(shoot_rel_abund)
  shoot_rel_abund <- left_join(shoot_rel_abund, map_seq_name, by="Useq" )
  shoot_rel_abund <- left_join(shoot_rel_abund, genome_data, by="taxon_oid")
  shoot_rel_abund <- left_join(gene_data[,c(1,ncol(gene_data))],shoot_rel_abund, by="taxon_oid" )
  shoot_flg22 <- rbind(shoot_flg22, get_flg22_percentages(shoot_rel_abund, rel_abund_end = num_condit, site = "shoot", condition=i, col_clust_start = which(colnames(shoot_rel_abund) == "clust1"), useq_col = which(colnames(shoot_rel_abund) == "Useq")) )
}

shoot_flg22$cond <- factor(shoot_flg22$cond, levels = unique(shoot$cond))

#get the agar flg22 seq data
agar_flg22 <- matrix(nrow = 0, ncol = 7)
for ( i in unique(isai_data$Rarefied$Map$condition ) ) {
  if(i == 'Inoculum'){
    next
  }
  agar_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$condition == i & isai_data$RelativeAbundance$Map$Tissue == "Agar")]
  agar_rel_abund <- as.data.frame(agar_rel_abund) 
  num_condit <- ncol(agar_rel_abund)
  agar_rel_abund$Useq <- row.names(agar_rel_abund)
  agar_rel_abund <- left_join(agar_rel_abund, map_seq_name, by="Useq" )
  agar_rel_abund <- left_join(agar_rel_abund, genome_data, by="taxon_oid")
  agar_rel_abund <- left_join(gene_data[,c(1,ncol(gene_data))],agar_rel_abund, by="taxon_oid" )
  agar_flg22 <- rbind(agar_flg22, get_flg22_percentages(agar_rel_abund, rel_abund_end = num_condit, site = "agar", condition=i, col_clust_start = which(colnames(agar_rel_abund) == "clust1"), useq_col = which(colnames(agar_rel_abund) == "Useq")) )
}

agar_flg22$cond <- factor(agar_flg22$cond, levels = unique(agar$cond))

total_flg22_data <- rbind(root_flg22, shoot_flg22, agar_flg22)
control_total <- total_flg22_data[total_flg22_data$cond %in% control_conditions,]


control_total <- dplyr::left_join(control_total, sgi_ros_data, by=c('full_seq'='PeptideSequence'))

#This will get the percentage of flg22 abundance I can functionally annotate
functional_data <- c()
for( i in control_conditions ) {
  for ( j in unique(control_total$place) ) {
    for ( l in unique(control_total$exp_id[control_total$cond == i & control_total$place == j])) {
      functional_data <- append(functional_data, c(i,j,
                                                 sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & is.na(control_total$ROS.Zscore)],na.rm = T)/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j], na.rm = T),
                                                 sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & !is.na(control_total$ROS.Zscore)], na.rm = T)/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j], na.rm = T),
                                                 sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & is.na(control_total$ROS.Zscore) & control_total$clust1 == 1], na.rm = T)/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & is.na(control_total$ROS.Zscore)], na.rm = T),
                                                 sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & is.na(control_total$ROS.Zscore) & control_total$clust2 == 1], na.rm = T)/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & is.na(control_total$ROS.Zscore)], na.rm = T),
                                                 sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & is.na(control_total$ROS.Zscore) & control_total$clust3 == 1], na.rm = T)/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & is.na(control_total$ROS.Zscore)], na.rm = T)
    ))
    }
  }
}

functional_data <- as.data.frame(matrix(functional_data, ncol=7, byrow = T))
colnames(functional_data) <- c("Cond", 'Place', "unknown", "known", "clade1_u", "clade2_u", 'clade3_u')

#plot known phenotypes
functional_data$known <- as.numeric(as.character(functional_data$known))
ggplot(functional_data, aes(Place,known))+
  geom_boxplot()+
  geom_point(position = position_jitter(width = .2))+
  labs(y="Percentage of flg22 Abundance Tested", x="Fraction")+
  #ylim(c(0,1))+
  theme_bw(base_size = 14)

#What clades are the unknown percentage from?
long_function <-functional_data %>%
  pivot_longer(cols = clade1_u:clade3_u, names_to = "clade_u", values_to = 'percentage')

long_function$percentage <- as.numeric(as.character(long_function$percentage))
pdf("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/185_trial/What_clade_are_we_missing.pdf", height = 4, width = 4)
ggplot(long_function, aes(clade_u,percentage, color=Place))+
  geom_point(position = position_jitter(width = .2))+
  labs(y="Percent of Untested flg22 variants", x="Fraction")+
  ylim(c(0,1))+
  theme_bw(base_size = 14)
dev.off()

immunogenicity_data <- c()
for( i in control_conditions ) {
  for ( j in unique(control_total$place) ) {
    for ( l in unique(control_total$exp_id[control_total$cond == i & control_total$place == j]))
      immunogenicity_data <- append(immunogenicity_data, c(i,j,l,
                                                         sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & control_total$ROS.Zscore > 0 & control_total$SGI.Zscore < 0 & control_total$fdr_pvalue < 0.05], na.rm = T)/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j], na.rm = T), #This is for the active
                                                         sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & control_total$ROS.Zscore > 0 & (control_total$SGI.Zscore > 0 | control_total$fdr_pvalue > 0.05 )],na.rm = T)/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j], na.rm = T), #deviant peptides
                                                         sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & control_total$ROS.Zscore < 0 & (control_total$SGI.Zscore > 0 | control_total$fdr_pvalue > 0.05 ) & !(control_total$full_seq %in% antag_seqs) ], na.rm = T)/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j], na.rm = T), #evading peptides
                                                         sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & control_total$full_seq %in% antag_seqs ])/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j], na.rm = T),
                                                         sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j & is.na(control_total$ROS.Zscore) ])/sum(control_total$flg22_abundance[control_total$exp_id == l & control_total$cond == i & control_total$place == j], na.rm = T)#antag flg22 sequences
    ))
  }
}

immunogenicity_data <- as.data.frame(matrix(immunogenicity_data, ncol=8, byrow = T))
colnames(immunogenicity_data) <- c("Cond", 'Place',"exp_id", "Immuno", "Deviant", "Evading", "Antag", "Not_tested")

immunogenicity_data <- immunogenicity_data[!is.na(immunogenicity_data$Immuno),]

long_immune <-immunogenicity_data %>%
  pivot_longer(cols = Immuno:Not_tested, names_to = "Fxn", values_to = 'percentage')

long_immune$percentage <- as.numeric(as.character(long_immune$percentage))

#pdf("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/185_trial/flg22_activity.pdf", height = 5, width = 5)
ggplot(long_immune, aes(Place,percentage, color=Fxn))+
  geom_boxplot()+
  geom_point(position = position_dodge(width = .75))+
  labs(y="Median Percentage of flg22 functional Type In Samples", x="Fraction")+
  theme_bw(base_size = 14)
#dev.off()

#THIS PART OF THE CODE BUILDS THE LOGO
logo_data <- c()
for ( i in unique(control_total$full_seq) ) {
  for ( j in c("root", "shoot", "agar") ) {
    logo_data <- append(logo_data, c(i, j, median(control_total$flg22_abundance[control_total$full_seq == i & control_total$place == j])))
  }
}
logo_data <- matrix(logo_data, ncol = 3, byrow = T)
logo_data <- as.data.frame(logo_data)
colnames(logo_data) <- c("seq", 'place', "med_abund")

logo_data$med_abund <- as.numeric(as.character( logo_data$med_abund ))
logo_data <- logo_data[logo_data$med_abund > 0,]

root <- apply(logo_data[logo_data$place == "root",], 1, function(x){
  rep(x[1], round(as.numeric(x[3])*1000) )
})
root_print <- c()
for ( i in root) {
  root_print <- append(root_print, i)
}
length(root_print)

shoot <- apply(logo_data[logo_data$place == "shoot",], 1, function(x){
  rep(x[1], round(as.numeric(x[3])*1000) )
})
shoot_print <- c()
for ( i in shoot) {
  shoot_print <- append(shoot_print, i)
}
length(shoot_print)

agar <- apply(logo_data[logo_data$place == "agar",], 1, function(x){
  rep(x[1], round(as.numeric(x[3])*1000) )
})
agar_print <- c()
for ( i in agar) {
  agar_print <- append(agar_print, i)
}
length(agar_print)

#I'm using the scripts from simple-pssm-builder
seq_data_all <- c()
root_print <- as.matrix(root_print, ncol=1)
seq_data_all <- sapply(root_print, function(x){ seq_data_all <- append(seq_data_all, strsplit(x, "") )})
seq_data_all <- as.character(unlist(seq_data_all))

seq_data <- cbind(root_print, matrix(seq_data_all, ncol = 22, byrow = T))
seq_data <- as.data.frame(seq_data)
seq_data$name <- "root"

table(seq_data[,16])

base_out <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/185_trial/root"
#regular_pssm(seq_data, "name", 2, 23, base_out)

seq_data_all <- c()
shoot_print <- as.matrix(shoot_print, ncol=1)
seq_data_all <- sapply(shoot_print, function(x){ seq_data_all <- append(seq_data_all, strsplit(x, "") )})
seq_data_all <- as.character(unlist(seq_data_all))

seq_data <- cbind(shoot_print, matrix(seq_data_all, ncol = 22, byrow = T))
seq_data <- as.data.frame(seq_data)
seq_data$name <- "root"

table(seq_data[,16])

base_out <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/185_trial/shoot"
#regular_pssm(seq_data, "name", 2, 23, base_out)

seq_data_all <- c()
agar_print <- as.matrix(agar_print, ncol=1)
seq_data_all <- sapply(agar_print, function(x){ seq_data_all <- append(seq_data_all, strsplit(x, "") )})
seq_data_all <- as.character(unlist(seq_data_all))

seq_data <- cbind(agar_print, matrix(seq_data_all, ncol = 22, byrow = T))
seq_data <- as.data.frame(seq_data)
seq_data$name <- "agar"

table(seq_data[,16])

base_out <- "/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/185_trial/agar"
#regular_pssm(seq_data, "name", 2, 23, base_out)




####Sick
sick_total <- total_flg22_data[total_flg22_data$cond %in% sick_conditions,]

sick_total <- dplyr::left_join(sick_total, sgi_ros_data, by=c('full_seq'='PeptideSequence'))

immunogenicity_data_bad <- c()
for( i in sick_conditions ) {
  for ( j in unique(sick_total$place) ) {
    for ( l in unique(sick_total$exp_id[sick_total$cond == i & sick_total$place == j] ) ) {
      immunogenicity_data_bad <- append(immunogenicity_data_bad, c(i,j,l,
                                                                 sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & sick_total$ROS.Zscore > 0 & sick_total$SGI.Zscore < 0 & sick_total$fdr_pvalue < 0.05], na.rm = T)/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T), #This is for the active
                                                                 sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & sick_total$ROS.Zscore > 0 & (sick_total$SGI.Zscore > 0 | sick_total$fdr_pvalue > 0.05 )],na.rm = T)/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T), #deviant peptides
                                                                 sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & sick_total$ROS.Zscore < 0 & (sick_total$SGI.Zscore > 0 | sick_total$fdr_pvalue > 0.05 ) & !(sick_total$full_seq %in% antag_seqs) ], na.rm = T)/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T), #evading peptides
                                                                 sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & sick_total$full_seq %in% antag_seqs ])/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T), #antag flg22 sequences
                                                                 sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & is.na(sick_total$ROS.Zscore) ])/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T)
    ))
    }
  }
}

immunogenicity_data_bad <- as.data.frame(matrix(immunogenicity_data_bad, ncol=8, byrow = T))
colnames(immunogenicity_data_bad) <- c("Cond", 'Place',"exp_id" , "Immuno", "Deviant", "Evading", "Antag", 'Not_tested')

long_immune_bad <-immunogenicity_data_bad %>%
  pivot_longer(cols = Immuno:Not_tested, names_to = "Fxn", values_to = 'percentage')

long_immune_bad$percentage <- as.numeric(as.character(long_immune_bad$percentage))

ggplot(long_immune_bad, aes(Place,percentage, color=Fxn))+
  geom_boxplot(outlier.color = NA)+
  geom_point(position = position_jitterdodge(jitter.width = .2))+
  labs(y="Percentage of flg22 Abundance Tested", x="Fraction")+
  #ylim(c(0,1))+
  theme_bw(base_size = 14)

#plot the healthy and sick data side by side
long_immune_bad$plant <- "200mM_NaCl"
long_immune$plant <- "Control"
total <- rbind(long_immune_bad, long_immune)


#pdf("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/185_trial/flg22_activity_env_health.pdf", height = 5, width = 7)
ggplot(total[total$plant == 'Control',], aes(Place,percentage, shape=Fxn))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitterdodge(jitter.width = .05), size=2.5)+
  labs(y="Relative Abundace of flg22 Variant Class", x="Fraction")+
  ggtitle("Flg22 Profile in 21 C and 1000 Pi Communities")+
  #ylim(c(0,1))+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = .5))

sum_control <- total[total$plant == 'Control',] %>%
  group_by(Place, Fxn) %>%
  dplyr::summarise(total_n = n())


#dev.off()

ggplot(total[total$plant == '200mM_NaCl',], aes(Place,percentage, shape=Fxn))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitterdodge(jitter.width = .1), size=2)+
  labs(y="Median Percentage of flg22 functional Type In Samples", x="Fraction")+
  ggtitle("Flg22 Profile in 200 mM NaCl Communities")+
  #ylim(c(0,1))+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = .5))

sum_control <- total[total$plant == '200mM_NaCl',] %>%
  group_by(Place, Fxn) %>%
  dplyr::summarise(total_n = n())

total$tukey_group <- paste(total$plant, total$Fxn, total$Place, sep = "_")


anova_model_control <- aov(percentage ~ tukey_group * Cond, data = total[total$plant == 'Control',])
summary(anova_model_control)

total_tukey <- HSD.test(anova_model_control, trt="tukey_group")
total_tukey

anova_model_sick <- aov(percentage ~ tukey_group, data = total[total$plant == '200mM_NaCl',])
summary(anova_model_sick)

total_tukey_sick <- HSD.test(anova_model_sick, trt="tukey_group")
total_tukey_sick


#Analyze the immunogenic sequences
control_total$plant <- "Control"
sick_total$plant <- "Salt"
total_flg22_data <- rbind(sick_total, control_total)
for ( i in unique(total_flg22_data$exp_id) ) {
  print(i)
  print(sum(total_flg22_data$flg22_abundance[total_flg22_data$exp_id == i]))
  flg22_abund_sample <- sum(total_flg22_data$flg22_abundance[total_flg22_data$exp_id == i])
  
  total_flg22_data$flg22_abundance[total_flg22_data$exp_id == i] <- total_flg22_data$flg22_abundance[total_flg22_data$exp_id == i]/flg22_abund_sample
}

total_flg22_data$plant <- factor(total_flg22_data$plant, levels = c("Control", 'Salt'))

ros_order <- sgi_ros_data$PeptideSequence[order(sgi_ros_data$ROS.Zscore, decreasing = T)]
total_flg22_data <- total_flg22_data[total_flg22_data$full_seq %in% ros_order,]
total_flg22_data$full_seq <- factor(as.character(total_flg22_data$full_seq),levels = ros_order[!duplicated(ros_order)])
total_flg22_data$deviant[ total_flg22_data$ROS.Zscore > 0 & (total_flg22_data$SGI.Zscore > 0 | total_flg22_data$fdr_pvalue > 0.05 )] <- 'Deviant'

#I need to normalize the relative abundance as is displayed in the previous figures

total_flg22_data <- total_flg22_data[total_flg22_data$ROS.Zscore > 0,]


library(grid)
library(gridExtra)
library(MuMIn)

trial <- lm( flg22_abundance ~ 0 + full_seq, data = total_flg22_data, na.action = "na.fail")
dredge(trial)
#complex model is the best model
summary(trial)
ag_trial <- lm( flg22_abundance ~ 0 + full_seq, data = total_flg22_data[total_flg22_data$place =="agar",], na.action = "na.fail")
sh_trial <- lm( flg22_abundance ~ 0 + full_seq, data = total_flg22_data[total_flg22_data$place =="shoot",], na.action = "na.fail")
rt_trial <- lm( flg22_abundance ~ 0 + full_seq, data = total_flg22_data[total_flg22_data$place =="root",], na.action = "na.fail")
summary(ag_trial)
summary(sh_trial)
summary(rt_trial)


#I am only going to look at the prevelant and abundant sequences from at least one fraction
#I am defining them as peptides that were significantly different than 0 in above model
prev_abund_flg22_seqs <- c("QRLSSGLRINSAKDDAAGLAIS","TRLSSGLKINSAKDDAAGLQIA", "ERLSTGKRINSAKDDAAGLAIA", "EKLSSGLRINRAADDAAGLSIS",
                           "QRLSSGLRINSAKDDAAGQAIA", "TRLSSGKRINSAADDAAGLAIS", "ERLSSGMRINSAKDDAAGQAIA", "EKLSSGFRINRAADDAAGLAIS",
                           "EKLSSGLRINRAGDDAAGLAIS", "ERLSSGMRINSAKDDAAGQAIA", "ERLSTGKRINSAKDDAAGLAIA")

prev_abund_data <- total_flg22_data[total_flg22_data$full_seq %in% prev_abund_flg22_seqs,]
trial <- lm( flg22_abundance ~ 0 + full_seq * place * plant, data = prev_abund_data, na.action = "na.fail")
prev_abund_data$label <- paste(prev_abund_data$place, prev_abund_data$plant, prev_abund_data$full_seq, sep = "_")
hsd_anova <- aov(flg22_abundance ~ label, data = prev_abund_data[prev_abund_data$plant == "Control",])
t <- HSD.test(hsd_anova, trt="label")
trial <- lm( flg22_abundance ~ 0 + full_seq * place, data = prev_abund_data[prev_abund_data$plant == "Control",], na.action = "na.fail")
summary(trial)

#pdf("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/185_trial/flg22_activity_all_samples-2.pdf", width = 8, height = 5)
p1<- ggplot(prev_abund_data, aes( full_seq, flg22_abundance, color=plant, shape=place))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitterdodge(jitter.width = .05), size=2)+
  #geom_label(aes(label=deviant, y=.5))+
  #facet_wrap(facets = "full_seq", scales = "free_y")+
  #stat_compare_means(comparisons = list(c("agar", "root"), c("agar", 'shoot')), label = "p.signif",method = "t.test", aes(x=place, color=plant))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust =.5))
dev.off()
total_flg22_data$sgi_signif[total_flg22_data$SGI.Zscore < 0 & total_flg22_data$fdr_pvalue < 0.01 ] <- "Signif"
total_flg22_data$sgi_signif[ is.na(total_flg22_data$sgi_signif)] <- "Not"


rgraph <- ggplot(total_flg22_data[total_flg22_data$full_seq %in% prev_abund_flg22_seqs,], aes(full_seq, y="ROS", fill=ROS.Zscore, color=sgi_signif ))+
  geom_tile(size=2)+
  scale_fill_gradient2(high="#662D91", mid="#219AE0", low = "#D9E021", midpoint = 3)+
  scale_color_manual(values = c(NA, "black"))

#sgraph <- ggplot(total_flg22_data[total_flg22_data$full_seq %in% prev_abund_flg22_seqs,], aes(full_seq, y="SGI", fill=SGI.Zscore))+
#  geom_tile()+
# scale_fill_gradient2(high="#D9E021", mid="#219AE0", low = "#662D91", midpoint = -1.5)
library(patchwork)
p1/(rgraph/sgraph)

gA <- ggplotGrob(p1)
gB <- ggplotGrob(rgraph)
#gC <- ggplotGrob(sgraph)
gA$widths <- gB$widths
#gC$widths <- gB$widths
grid.newpage()
grid.draw(arrangeGrob(gA,gB, heights = c(1,.2)))
dev.off()

levels(total_flg22_data$full_seq)
total_flg22_data$place <- factor(total_flg22_data$place, levels=c("agar", "root", "shoot"))
total_flg22_data$sick <- factor(total_flg22_data$sick, levels=c("Control", "Challenged"))

total_flg22_data <- total_flg22_data[total_flg22_data$ROS.Zscore > 0 & total_flg22_data$cond != "0Pi",]

trial <- lm( flg22_abundance ~ 0 + full_seq * sick * place, data = total_flg22_data, na.action = "na.fail")
dredge(trial)
summary(trial)
trial <- aov(flg22_abundance ~ 0 + full_seq * sick * place, data = total_flg22_data)
summary(trial)
dredge_top_model <- lm(flg22_abundance ~ 0 + full_seq * sick + place + full_seq:place, data = total_flg22_data, na.action = "na.fail")
summary(dredge_top_model)


#Look at specific microbes with these sequences
root_rel_abund <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$Tissue == "Root")]
root_rel_abund <- as.data.frame(root_rel_abund) 
num_condit <- ncol(root_rel_abund)
root_rel_abund$Useq <- row.names(root_rel_abund)
root_rel_abund <- left_join(root_rel_abund, map_seq_name, by="Useq" )
#root
root_rel_abund <- left_join(root_rel_abund, genome_data, by="taxon_oid")
root_rel_abund <- left_join(gene_data[,c(1,ncol(gene_data))],root_rel_abund, by="taxon_oid" )
root_rel_abund[root_rel_abund$full_seq=="QRLSSGLRINSAKDDAAGLAIS" & !is.na(root_rel_abund$Useq),]
root_rel_abund[root_rel_abund$full_seq=="TRLSSGLKINSAKDDAAGLQIA" & !is.na(root_rel_abund$Useq),]

isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$condition %in% control_conditions & isai_data$RelativeAbundance$Map$Tissue == "Root")]



### Look at the plants from lower salt concentrations
sick_conditions <- c("150NaCl", "100NaCl", "50NaCl")
sick_total <- total_flg22_data[total_flg22_data$cond %in% sick_conditions,]

sick_total <- dplyr::left_join(sick_total, sgi_ros_data, by=c('full_seq'='PeptideSequence'))

immunogenicity_data_bad <- c()
for( i in sick_conditions ) {
  for ( j in unique(sick_total$place) ) {
    for ( l in unique(sick_total$exp_id[sick_total$cond == i & sick_total$place == j] ) ) {
      immunogenicity_data_bad <- append(immunogenicity_data_bad, c(i,j,l,
                                                                   sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & sick_total$ROS.Zscore > 0 & sick_total$SGI.Zscore < 0 & sick_total$fdr_pvalue < 0.05], na.rm = T)/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T), #This is for the active
                                                                   sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & sick_total$ROS.Zscore > 0 & (sick_total$SGI.Zscore > 0 | sick_total$fdr_pvalue > 0.05 )],na.rm = T)/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T), #deviant peptides
                                                                   sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & sick_total$ROS.Zscore < 0 & (sick_total$SGI.Zscore > 0 | sick_total$fdr_pvalue > 0.05 ) & !(sick_total$full_seq %in% antag_seqs) ], na.rm = T)/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T), #evading peptides
                                                                   sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & sick_total$full_seq %in% antag_seqs ])/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T), #antag flg22 sequences
                                                                   sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j & is.na(sick_total$ROS.Zscore) ])/sum(sick_total$flg22_abundance[sick_total$exp_id == l & sick_total$cond == i & sick_total$place == j], na.rm = T)
      ))
    }
  }
}

immunogenicity_data_bad <- as.data.frame(matrix(immunogenicity_data_bad, ncol=8, byrow = T))
colnames(immunogenicity_data_bad) <- c("Cond", 'Place',"exp_id" , "Immuno", "Deviant", "Evading", "Antag", 'Not_tested')

long_immune_bad <-immunogenicity_data_bad %>%
  pivot_longer(cols = Immuno:Not_tested, names_to = "Fxn", values_to = 'percentage')

long_immune_bad$percentage <- as.numeric(as.character(long_immune_bad$percentage))

ggplot(long_immune_bad, aes(Place,percentage, shape=Fxn))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(position = position_jitterdodge(jitter.width = .1), size=2)+
  labs(y="Median Percentage of flg22 functional Type In Samples", x="Fraction")+
  ggtitle("Flg22 Profile in NaCl Communities")+
  #ylim(c(0,1))+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = .5))

#long_immune_bad <- long_immune_bad[long_immune_bad$Cond == "150NaCl",]
long_immune_bad$tukey_group <- paste(long_immune_bad$Fxn, long_immune_bad$Place, sep = "_")
anova_model_control <- aov(percentage ~ tukey_group*Cond, data = long_immune_bad)
summary(anova_model_control)

total_tukey <- HSD.test(anova_model_control, trt="tukey_group")
total_tukey

