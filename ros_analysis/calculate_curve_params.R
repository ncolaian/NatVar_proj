#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) # What you need to pass in arguements

#this code will be used to find the parameters of the curves in the data.
library(getopt)
library(BayesianTools)
library(tidyr)
library(Bolstad2) #used to find the area under the curve
library(coda)

params = matrix(c(
  "input_file", "i", 1, "character",
  "out_dir", "o", 1, "character",
  "plate_name", "n", 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(params)

### TESTING ###
#opt$input_file <- "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/ros_8_27-30/8_27_18/141414.txt"
#opt$out_dir <- "/Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/"
#opt$plate_name <- "141414"
#################

format_data <- function(df_path) {
  df <- read.csv(df_path, stringsAsFactors = F, header = T)
  df <- gather(df, Sample, Photons, A1:A12 )
  return(df)
}

full_likelihood <- function(params) {
  theta1 <- params[1]
  theta2 <- params[2]
  theta3 <- params[3]
  theta4 <- params[4]
  theta5 <- params[5]
  
  #find change point
  pred <- c()
  for (i in x) {
    if ( i < theta4) {
      pred <- append(pred, theta1*exp(-exp(theta3*(i+theta2))))
    }
    else{
      #constant to link the two curves (theta1*exp(-1*exp(theta3*(theta4+theta2)))*exp(theta5*(theta4)))
      pred <- append(pred, exp(-i*theta5)*(theta1*exp(-1*exp(theta3*(theta4+theta2)))*exp(theta5*(theta4))))
    }
  }
  singlelikelihoods = dnorm(y, mean = pred, sd = 1/params[6]^2, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

get_parameter_values_2 <- function(iter = 500000) {
  #get rid of time 0
  x <- x[2:length(x)]
  y <- y[2:length(y)]
  
  #nrchains can be changed for testing purposes
  #mcmc settings
  settings <- list(iterations = iter,nrChains=1, adapt = F, DRlevels = 1, gibbsProbabilities = NULL, temperingFunction = NULL, optimize = F, message = FALSE)
  #mins and max are taken from paper
  sprior <- createUniformPrior(lower = c(0,-10, -10,2,0,0), upper = c(3.2, 0, 0, 30, 1,50))
  setup <- createBayesianSetup(full_likelihood, prior =sprior)
  out <- runMCMC(bayesianSetup = setup, sampler = "Metropolis", settings = settings)
  return(out)
}

### MAIN ###
data <- format_data(opt$input_file)

curve_params <- c()
for ( i in unique(data$Name) ) {
  for ( j in unique(data$Sample) ) {
    x <- data$Time[data$Name == i & data$Sample == j]
    y <- data$Photons[data$Name == i & data$Sample == j]
    max_photon <- max(y)
    auc_1 <- sintegral(x, y)$int
    #rescale the data
    x <- x/60 #time into minutes
    y <- y/max_photon #rescale to max ie 0-1
    out <- get_parameter_values_2()
    info <- MAP(out)
    auc_2 <- sintegral(x, y)$int
    curve_params <- append(curve_params, c(i,j,opt$plate_name, info$parametersMAP,auc_1, auc_2))
  }
}

curve_mat <- matrix(curve_params, ncol = 11, byrow = T)
colnames(curve_mat) <- c("Name", "Sample", "Plate", "P1", "P2", "P3", "P4", "P5", "P6", "FullAuc", "ScaledAuc" )

write.csv(curve_mat, paste(opt$out_dir, paste(opt$plate_name, ".curve_params.csv", sep = "" ), sep = "/"), row.names = F, quote = F)

