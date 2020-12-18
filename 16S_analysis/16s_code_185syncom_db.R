##### THIS IS THE DATA FROM FINKEL ET AL
#They do a soil experiment and then go through their SynCom.

#I will basically be repeating the right side of figure 4B, but with the soil data
soil_consensus <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/reviewer_comments/finkel_et_al_suppdata/soil_relative_abundance_info.csv")
syncom185 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/reviewer_comments/finkel_et_al_suppdata/185_inoculum_info.csv")

#some 16S abundance is not assigned to an ASV, thus we will normalize it all the Relative abundance values to 1 for convenience
soil_consensus$Cumulative.relative.abuncance.in.soil <- soil_consensus$Cumulative.relative.abuncance.in.soil/sum(soil_consensus$Cumulative.relative.abuncance.in.soil)
soil_consensus$Cumulative.relative.abuncance.in.root <- soil_consensus$Cumulative.relative.abuncance.in.root/sum(soil_consensus$Cumulative.relative.abuncance.in.root)
soil_consensus$Cumulative.relative.abuncance.in.shoot <-soil_consensus$Cumulative.relative.abuncance.in.shoot/sum(soil_consensus$Cumulative.relative.abuncance.in.shoot)
sum(soil_consensus$Cumulative.relative.abuncance.in.soil)
sum(soil_consensus$Cumulative.relative.abuncance.in.shoot)
sum(soil_consensus$Cumulative.relative.abuncance.in.root)
soil_consensus[soil_consensus$Phylum == "Proteobacteria",]
table(soil_consensus$Phylum)
table(soil_consensus$Class[soil_consensus$Phylum == "Proteobacteria"])

#my database data
arab_db <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/final_figure_antag_paper/supplemental_data/supplemental_table1.csv")
table(arab_db$IMG_Class)[names(table(arab_db$IMG_Class) ) %in% names(table(soil_consensus$Class) )]

#get a plot taxa to get the taxa in the db from the taxa found in the soil
soil_consensus$plot_taxa <- soil_consensus$Class
soil_consensus$plot_taxa <- as.character(soil_consensus$plot_taxa)
#classes in the DB not in the soil census
table(arab_db$IMG_Class)[!names(table(arab_db$IMG_Class) ) %in% names(table(soil_consensus$Class) )]

# we will manually go in and identify the classes that are not identical between the datasets to make them comparable
# In the soil census data betaprroteobacteria are an order, which we will adjust
soil_consensus$plot_taxa[soil_consensus$Order == "Betaproteobacteriales"] <- "Betaproteobacteria"

#flavobacteriia is not a class in the soil consensus. we will take all orders in flavobacteriia and manually annotate them
#only one order Flavobacteriales
soil_consensus$plot_taxa[soil_consensus$Order == "Flavobacteriales"] <- "Flavobacteriia"

#Sphingobacteriia is not a class in the soil census - the order is Sphingobacteriales
soil_consensus$plot_taxa[soil_consensus$Order == "Sphingobacteriales"] <- "Sphingobacteriia"

table(soil_consensus$Order[soil_consensus$Class == "Bacteroidia"])

#cytophagales is a cytophagia
soil_consensus$plot_taxa[soil_consensus$Order == "Cytophagales"] <- "Cytophagia"

#Deinococci is not a class in the soil census. The order is not found.

table(soil_consensus$Order)


#I will create the figures from finkel et al. first
#table(soil_consensus$plot_taxa)
#soil_consensus$plot_taxa[!soil_consensus$plot_taxa %in% arab_db$IMG_Class] <- "Other" 
#soil_consensus[soil_consensus$plot_taxa == "Other" & soil_consensus$Cumulative.relative.abuncance.in.root > 0.01,]
#soil_consensus$plot_taxa[soil_consensus$plot_taxa == "Other"] <- as.character(soil_consensus$Phylum[soil_consensus$plot_taxa == "Other"])

soil_consensus$Phylum <- soil_consensus$Phylum %>% factor(levels = c("Acidobacteria","Actinobacteria","Bacteroidetes",
                                               "Chloroflexi","Cyanobacteria","Firmicutes",
                                               "Gemmatimonadetes","Patescibacteria","Proteobacteria",
                                               "Verrucomicrobia","Other"))
soil_consensus$Phylum[is.na(soil_consensus$Phylum)] <- "Other"

soil_consensus$plot_taxa[soil_consensus$Phylum == "Other" & !soil_consensus$plot_taxa %in% arab_db$IMG_Class] <- "Not-In-Database" 

table(soil_consensus$plot_taxa)

library(dplyr)
plot_soil_census <- soil_consensus %>%
  dplyr::group_by(Phylum) %>%
  summarise(total_soil = sum(Cumulative.relative.abuncance.in.soil),
            total_root = sum(Cumulative.relative.abuncance.in.root),
            total_shoot = sum(Cumulative.relative.abuncance.in.shoot))

colnames(plot_soil_census) <- c("Phylum", "Soil", "Root", "Shoot")
plot_soil_census <- pivot_longer(plot_soil_census, Soil:Shoot, names_to = 'location', values_to = "rel_abund")

plot_soil_census$location <- plot_soil_census$location %>% factor(levels = c("Soil", "Root", "Shoot"))

library(paletteer)

t_color <- c(paletteer_d(palette = "ggsci::category20_d3",n = 10), "Black")

ggplot(plot_soil_census, aes(location, y=rel_abund, fill=Phylum) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Phylum")+
  scale_fill_manual(values = t_color)+
  theme_classic()


#Trial of the class dataset
plot_soil_census_class <- soil_consensus %>%
  dplyr::group_by(plot_taxa) %>%
  summarise(total_soil = sum(Cumulative.relative.abuncance.in.soil),
            total_root = sum(Cumulative.relative.abuncance.in.root),
            total_shoot = sum(Cumulative.relative.abuncance.in.shoot))

#Reduce some of the low abundant classes
plot_soil_census_class$plot_taxa <- as.character(plot_soil_census_class$plot_taxa)
plot_soil_census_class$plot_taxa[ !plot_soil_census_class$plot_taxa %in% arab_db$IMG_Class] <- "Not-In-Database"

plot_soil_census_class <- plot_soil_census_class %>%
  dplyr::group_by(plot_taxa) %>%
  summarise(total_soil = sum(total_soil),
            total_root = sum(total_root),
            total_shoot = sum(total_shoot))

colnames(plot_soil_census_class) <- c("Class", "Soil", "Root", "Shoot")
plot_soil_census_class <- pivot_longer(plot_soil_census_class, Soil:Shoot, names_to = 'location', values_to = "rel_abund")
plot_soil_census_class$location <- plot_soil_census_class$location %>% factor(levels = c("Soil", "Root", "Shoot"))

plot_soil_census_class$Class <- plot_soil_census_class$Class %>% factor(levels = c("Actinobacteria", "Alphaproteobacteria", "Bacilli", "Betaproteobacteria", 'Cytophagia', "Deinococci", 'Flavobacteriia', 'Gammaproteobacteria', 'Sphingobacteriia', 'Not-In-Database', "Non-Isolate-RA"))
t_color <- c(paletteer_d(palette = "ggsci::category20_d3",n = 8), "Black")

ggplot(plot_soil_census_class, aes(location, y=rel_abund, fill=Class) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Class")+
  scale_fill_manual(values = t_color)+
  theme_classic()

# I will now get the inoculum data

#I WANT TO VISUALIZE ALL THIS DATA IN THE CONTEXT OF CLASS ISOLATED MICROBES
root_microbe_not_in_db <- plot_soil_census_class$rel_abund[plot_soil_census_class$Class == "Not-In-Database" & plot_soil_census_class$location == 'Root']

isai_data <- readRDS("/Users/nicholascolaianni/Documents/dangl_lab/Isai_data/stress_experiment_data/dat_amplicon_useq97_4stresses.RDS")
inoculum <- isai_data$RelativeAbundance$Tab[,which(isai_data$RelativeAbundance$Map$Tissue == "Inoculum")]

inoculum_mean <- rowMeans(inoculum)
inoculum_mean <- inoculum_mean/sum(inoculum_mean)

#this gets the taxa information for all the Useqs
classes <- c()
for ( i in isai_data$RelativeAbundance$Tax$Taxonomy) {
  class <- strsplit(i, ";")[[1]][4]
  class <- gsub(" c__", "", class)
  classes <- append(classes, class)
}

all.equal.character(names(inoculum_mean), as.character(isai_data$RelativeAbundance$Tax$ID))
#the order is the same

inoculum_data <- as.data.frame( matrix(c(inoculum_mean, classes), ncol = 2, byrow = F) )
inoculum_data$location <- "Inoculum"
inoculum_data$V1 <- as.numeric(as.character(inoculum_data$V1))
inoculum_data$V1 <- inoculum_data$V1 * (1-root_microbe_not_in_db)
inoculum_data <- inoculum_data[,c(2,3,1)]
colnames(inoculum_data) <- colnames(plot_soil_census_class)

tw_data <- as.data.frame(matrix(c("Non-Isolate-RA", "Inoculum", root_microbe_not_in_db), ncol = 3))
colnames(tw_data) <- colnames(plot_soil_census_class)
inoculum_data <- rbind(inoculum_data, tw_data)
inoculum_data$Class <- as.character(inoculum_data$Class)
inoculum_data$Class[!inoculum_data$Class %in% levels(plot_soil_census_class$Class)] <- "Not-In-Database"

plot_soil_census_class <- rbind(plot_soil_census_class, inoculum_data)
plot_soil_census_class$rel_abund <- as.numeric(as.character(plot_soil_census_class$rel_abund))

t_color <- c(paletteer_d(palette = "ggsci::category20_d3",n = 9), "Black", "White")
ggplot(plot_soil_census_class, aes(location, y=rel_abund, fill=Class) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Class")+
  scale_fill_manual(values = t_color)+
  theme_classic()

# i will now get the community data
arab_db_class <- table(arab_db$IMG_Class)
arab_db_class <- arab_db_class/sum(arab_db_class)

arab_db_class <- matrix(c(names(arab_db_class), arab_db_class), ncol = 2, byrow = F)
arab_db_class <- as.data.frame(arab_db_class)
arab_db_class$names <- "Arab_Database"
arab_db_class$V2 <- as.numeric(as.character(arab_db_class$V2))
arab_db_class <- arab_db_class[,c(1,3,2)]
colnames(arab_db_class) <- colnames(plot_soil_census_class)
#scale the RA based on the community composition that is representative of wild soil
arab_db_class$rel_abund <- arab_db_class$rel_abund * (1-root_microbe_not_in_db)
#add the data to make it comparable
tw_data <- as.data.frame(matrix(c("Non-Isolate-RA", "Arab_Database", root_microbe_not_in_db), ncol = 3))
colnames(tw_data) <- colnames(plot_soil_census_class)
arab_db_class <- rbind(arab_db_class, tw_data)

plot_soil_census_class <- rbind(plot_soil_census_class, arab_db_class)
plot_soil_census_class$rel_abund <- as.numeric(as.character(plot_soil_census_class$rel_abund))

plot_soil_census_class$location <- plot_soil_census_class$location %>% factor(levels = c('Soil', "Root", "Shoot", "Inoculum", "Arab_Database") )

t_color <- c(paletteer_d(palette = "ggsci::category20_d3",n = 9), "Black", "White")
t_color_t <- c(paletteer_d(palette = "ggsci::category20_d3",n = 9), "Black", "White")
#I will be keeping the colors the same from previous figures
t_color[2] <- t_color[1]
t_color[1] <- '#9400D3'
t_color[3] <- "#008000"
t_color[4] <- "#DFE120"
t_color[5] <- t_color_t[8]
t_color[7] <- "#FF8000"
t_color[8] <- "#FF0000"
t_color[9] <- t_color_t[7]

ggplot(plot_soil_census_class, aes(location, y=rel_abund, fill=Class) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Class")+
  scale_fill_manual(values = t_color)+
  theme_classic()

### THIS is the phylum analysis

############
#          #
# OUTDATED #
#          #
############

#This is the data directly from the paper's figure
soil_consensus <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/reviewer_comments/finkel_et_al_suppdata/census_data.csv")

for (i in unique(soil_consensus$Fraction)) {
  soil_consensus$Abundance[soil_consensus$Fraction == i] <- soil_consensus$Abundance[soil_consensus$Fraction == i]/sum(soil_consensus$Abundance[soil_consensus$Fraction == i])
}
soil_consensus$Taxon <- soil_consensus$Taxon %>% factor(levels = c("Acidobacteria","Actinobacteria","Bacteroidetes",
                                            "Chloroflexi","Cyanobacteria","Firmicutes",
                                            "Gemmatimonadetes","Patescibacteria","Proteobacteria",
                                            "Verrucomicrobia","Other"))
soil_consensus$Fraction <- soil_consensus$Fraction %>% factor(levels = c("Soil", "Root", "Shoot"))

ggplot(soil_consensus, aes(Fraction, y=Abundance, fill=Taxon) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Phylum")+
  scale_fill_manual(values = t_color)+
  theme_classic()

census_phylum <- sapply(syncom185$Taxonomic.Classification, function(x){
  strsplit(as.character(x), split = ";")[[1]][2]
})

census_phylum <- gsub("p__", "", census_phylum)
census_phylum <- census_phylum[!is.na(census_phylum)]
table(census_phylum)

############
lundberg <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/Paulo_MAE/nature21417-s2/2012_lundberg.csv")
lundberg$Phylum <- as.character(lundberg$Phylum)
lundberg$Phylum[!(lundberg$Phylum %in% soil_consensus$Taxon)] <- "Other"

lundberg <- lundberg %>%
  dplyr::group_by(Phylum, fraction) %>%
  summarise(total = sum(Ra))

lundberg <- as.data.frame(lundberg)
lundberg <- lundberg[,c(1,3,2)]
colnames(lundberg) <- colnames(t)

lundberg$Freq <- lundberg$Freq / 100
lundberg$Var1 <- lundberg$Var1 %>% factor(levels = c("Acidobacteria","Actinobacteria","Bacteroidetes",
                                                                   "Chloroflexi","Cyanobacteria","Firmicutes",
                                                                   "Gemmatimonadetes","Patescibacteria","Proteobacteria",
                                                                   "Verrucomicrobia","Other"))
lundberg$data <- as.character(lundberg$data)
lundberg$data[lundberg$data == "S"] <- "Lundberg_Soil"
lundberg$data[lundberg$data == "R"] <- "RZ"
lundberg$data[lundberg$data == "EC"] <- "Lundberg_EC"
lundberg <- lundberg[lundberg$data %in% c('Lundberg_Soil', "Lundberg_EC"),]
lundberg$data <- factor(lundberg$data, levels = c("Lundberg_Soil", "Lundberg_EC"))

soil_consensus <- soil_consensus[,c(1,3,2)]
colnames(lundberg) <- colnames(soil_consensus)

soil_consensus <- rbind(lundberg,soil_consensus)

ggplot(soil_consensus, aes(Fraction, y=Abundance, fill=Taxon) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Phylum")+
  scale_fill_manual(values = t_color)+
  theme_classic()

#Now we will add the 185-member inoculum
syncom_185 <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/reviewer_comments/finkel_et_al_suppdata/data_Fig4B.csv")
syncom_185 <- syncom_185[syncom_185$typebyTissue == "Inoculum",]
syncom_185$typebyTissue <- as.character(syncom_185$typebyTissue)
syncom_185$typebyTissue[syncom_185$typebyTissue == "Inoculum"] <- '185-SynCom_Inoculum'

syncom_185 <- syncom_185[,c(1,3,2)]
colnames(syncom_185) <- colnames(soil_consensus)

soil_consensus <- rbind(soil_consensus, syncom_185)

ggplot(soil_consensus, aes(Fraction, y=Abundance, fill=Taxon) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Phylum")+
  scale_fill_manual(values = t_color)+
  theme_classic()

#Look into the Arabidopsis Community
arab_db <- table(arab_db$IMG_Phylum)
arab_db <- arab_db/sum(arab_db)

arab_db <- matrix(c(names(arab_db), arab_db), ncol = 2, byrow = F)
arab_db <- as.data.frame(arab_db)
arab_db$names <- "Arab_Database"

colnames(arab_db) <- colnames(soil_consensus)
soil_consensus <- rbind(soil_consensus, arab_db)
soil_consensus$Taxon[!(soil_consensus$Taxon %in% c("Acidobacteria","Actinobacteria","Bacteroidetes",
                                                   "Chloroflexi","Cyanobacteria","Firmicutes",
                                                   "Gemmatimonadetes","Patescibacteria","Proteobacteria",
                                                   "Verrucomicrobia","Other"))] <- "Other"
soil_consensus$Abundance <- as.numeric(as.character(soil_consensus$Abundance))

ggplot(soil_consensus, aes(Fraction, y=Abundance, fill=Taxon) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Phylum")+
  scale_fill_manual(values = t_color)+
  theme_classic()


#This is the final figure for now


ggplot(soil_consensus[!soil_consensus$Fraction %in% c('Lundberg_Soil', "Lundberg_EC"),], aes(Fraction, y=Abundance, fill=Taxon) )+
  geom_bar(stat="identity", position = "stack", width = .5)+
  labs(y="Relative Abundance", x="Dataset", fill="Phylum")+
  scale_fill_manual(values = t_color)+
  theme_classic()


