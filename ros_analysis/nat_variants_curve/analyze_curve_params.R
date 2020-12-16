#this script will be used to perform a pca plot. I will be performing a pca to identify the important measurements in the curves.

library(FactoMineR)
curve_params <- read.csv("/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/nat_variants_curve/combined.curve.params.csv", header = T)

res <- MFA(curve_params, group = c(1,1,1,1,1,1,1,1,1,1), type = c("n", "n", "s", "s", "s", "s", "s", "s", "s", "s"), ncp = 3,
           name.group = c("Sample", "Plate", "P1", "P2", "P3", "P4", "P5", "P6", "FullAuc", "ScaledAuc"))
curve_params$FullAuc <- log10(curve_params$FullAuc)
res.pca <- PCA(curve_params[,c(2,5:(ncol(curve_params)-2), ncol(curve_params))], quali.sup = c(1), graph = T, ncp =4, scale.unit = T)
plot.PCA(res.pca, choix = c("ind"), label = "quali", invisible = "ind")
dimdesc(res.pca)
plotellipses(res.pca)
res.pca$var
library(corrplot)
cor.mat <- cor(curve_params[,c(5:ncol(curve_params))])
corrplot(cor.mat)

library(ggbiplot)
ggbiplot(res.pca)
library(reshape2)
params_long <- melt(curve_params[,c(2,5:(ncol(curve_params)))], id.vars = c("Name"), variable.name = "Param", value.name = "Param_Val")

median_param_vals <- with(params_long, tapply(Param_Val,list(Name,Param), mean, na.rm=T))
res.pca <- prcomp(median_param_vals[,c(1,2,3,4,5,7)], center = T, scale. = T)
summary(res.pca)
ggbiplot(res.pca, labels = row.names(median_param_vals), ellipse = T, choices = c(1,3))
#only get the peptides that produce a signal above PtaDA
median_param_vals <- median_param_vals[order(median_param_vals[,7], decreasing = T),]
active_param_vals <- median_param_vals[1:(which(row.names(median_param_vals) == "PtaDA")-1),]
top.pca <- prcomp(active_param_vals[,c(1,2,3,4,5,7)], center = T, scale. = T)
ggbiplot(top.pca, choices = c(1,2))+
  geom_text(aes(label=row.names(active_param_vals), color=active_param_vals[,7]), hjust = 0, vjust=0)+
  scale_colour_gradient(low = "Gold", high = "firebrick")

top.pca$x
