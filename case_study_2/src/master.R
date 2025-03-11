#.rs.restartR()
library(Rcpp)
library(tools)
#compileAttributes("/Users/evsi8432/Documents/Research/momentuHMM")
#package_native_routine_registration_skeleton("/Users/evsi8432/Documents/Research/momentuHMM")
#install.packages("/Users/evsi8432/Documents/Research/momentuHMM",
#                 repos=NULL,
#                 type="source")
#library(momentuHMM)
#detach("package:momentuHMM", unload=TRUE)
library(momentuHMM)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(GGally)
library(signal)
library(oce)
library(pROC)
library(latex2exp)
library(randomForest)
library(caret)
library(e1071)
library(logistf)

directory <- "/Users/evsi8432/Documents/Research/PHMM/case_study_2"
setwd(directory)

# set options

args = commandArgs(trailingOnly=TRUE)
 
K <- 4 #as.numeric(args[1]) # number of cross-validations (one means just do all the data)
lambda <- 0.01 #as.numeric(args[2]) # lambda for paper
num_seeds <- 10 #as.numeric(args[3]) # number of random seeds

# create directories
dir.create(directory, showWarnings = FALSE)
dir.create(paste0(directory,"/params"), showWarnings = FALSE)
dir.create(paste0(directory,"/plt"), showWarnings = FALSE)

# set seed
set.seed(1)

# load in data
df <- data.frame(fread("dat/case_study_2_data.csv"))
df_echo <- data.frame(fread("../../dat/Final_Data_fine.csv"))
df_echo <- df_echo[df_echo$divenum %in% unique(df$divenum),]
df <- merge(df,df_echo[,c("stime","divenum","echo.steady","echo.rapid","crunch")], 
            all.x=T,all.y=F)

# create cross-validation groups
print("seperating train and test set...")
source("src/make_test_train.R")

# plot the data
print("plotting data...")
source("src/plot_data.R")

# initialize model lists
models_base <- list()
models_rf <- list()
models_svm <- list()
models_lr <- list()
models_PHMM <- list()

probs_PHMM <- list()
probs_rf <- list()
probs_svm <- list()
probs_lr <- list()
probs_base <- list()

AUCs_base <- rep(0,K)
AUCs_rf <- rep(0,K)
AUCs_svm <- rep(0,K)
AUCs_lr <- rep(0,K)
AUCs_PHMM <- rep(0,K)

conf_matrices_base <- list()
conf_matrices_rf <- list()
conf_matrices_svm <- list()
conf_matrices_lr <- list()
conf_matrices_PHMM <- list()

for(k in 1:K){
  
  train_dives <- train_sets[[k]]
  test_dives <- test_sets[[k]]
  
  # fit baseline
  print("fitting baseline...")
  source("src/fit_base.R") 
  models_base[[k]] <- base_model
  models_rf[[k]] <- rf_model
  models_svm[[k]] <- svm_model
  models_lr[[k]] <- lr_model
  
  # evaluate baseline
  print("evaluating baseline...")
  source("src/eval_base.R")  
  
  # fit the PHMM
  print("fitting PHMM...")
  best_hmm <- NULL
  max_ll <- -Inf
  for(rand_seed in 1:num_seeds){
    print(paste("working on fold",k,"of",K,"and seed",rand_seed,"of",num_seeds))
    source("src/fit_PHMM.R")
    if(-hmm$mod$minimum > max_ll){
      best_hmm <- hmm
      max_ll <- -hmm$mod$minimum
      print("new best hmm")
    }
  }
  hmm <- best_hmm
  models_PHMM[[k]] <- hmm

  # evaluate PHMM
  print("evaluating PHMM...")
  source("src/eval_PHMM.R")
  
  # plot hmm results
  print("plotting PHMM...")
  source("src/plot_PHMM.R")
}

# summarize cross-validation results
print("summarizing results...")
source("src/summarize_results.R")

# plot AUCs
source("src/plot_AUCs.R")
