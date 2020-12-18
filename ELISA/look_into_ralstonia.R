#Ralstonia ELISA data

elis_data <- read.csv('/Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/elisa_res/kate_december_final_antag_data/ralstoonia_elisa_data.csv')

library(ggplot2)
library(multcomp)
library(agricolae)

elis_data$Peptide <- factor(elis_data$Peptide, levels = c("Pa22", "Pa20", "1186", '1410', "5001", 
                                                          "5002", '5003', "5004", "5005"))

conf_aov <- aov(A60_fold ~ Peptide, elis_data)
summary(conf_aov)
t <- glht(conf_aov, linfct=mcp(Peptide="Dunnett"))
summary(t)

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

tukey_test <- HSD.test(conf_aov, trt="Peptide")

ggplot(elis_data, aes(Peptide, A60_fold))+
  geom_boxplot()+
  geom_point()+
  labs(y="A650 Fold Induction")+
  theme_bw(base_size = 16)
