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
library(mclust)
library(data.table)
library(RcppHungarian)
library(mvtnorm)
library(randomForest)
library(caret)
library(e1071)
library(nnet)

setwd("/Users/evsi8432/Documents/Research/PHMM/case_study_1/src/HMM")
set.seed(1)

# get options
opt_file <- "logMDDD_1-1-1_dd-30_2023-10-23.R"
source(paste0('../opt/',opt_file))

# define whales
whales <- c("none","A100a","A100b","A113a","A113b",
            "D21a","D21b","D26a","D26b",
            "I107a","I107b","I129","I145a","I145b",
            "L87","L88","R48a","R48b","R58a","R58b")

make_title <- function(start,end){
  title <- paste0(start,statesPerBehaviour[1])
  for(nstates in statesPerBehaviour[2:3]){
    title <- paste0(title,"-",nstates)
  }
  for(feature in names(dist)){
    title <- paste0(title,"-",feature)
  }
  if(length(sex) > 1){
    title <- paste0(title,"_all")
  } else {
    title <- paste0(title,"_",sex)
  }
  title <- paste0(title,"_",end)
  return(title)
}

for(holdout_whale in whales){
  
  # get data and remove held-out whale
  Data <- data.frame(fread("../../dat/case_study_1_data.csv"))
  if(holdout_whale != "none"){
    train_data <- Data[!(Data$ID %in% holdout_whale) & Data$knownState != 4,]
    test_data <- Data[Data$ID %in% holdout_whale & Data$knownState != 4,]
  }
  
  # get weights
  class_counts <- table(train_data$knownState)
  class_weights <- 1 / class_counts
  weights <- as.numeric(class_weights[as.character(train_data$knownState)])  # Assign weights based on class
  weights <- weights / sum(weights) * length(weights)  # Normalize weights
  
  # Train Random Forest model
  rf_model <- randomForest(
    factor(knownState) ~ logMDDD.x + logMDDD.y, 
    data = train_data,
    ntree=500, 
    importance = TRUE,
    weights = weights
  )
  conf_matrix <- matrix(0, nrow = 3, ncol = 3)
  rownames(conf_matrix) <- c("True Resting", "True Travelling", "True Foraging")
  colnames(conf_matrix) <- c("Predicted Resting", "Predicted Travelling", "Predicted Foraging")
  
  prob_labs <- matrix(0,ncol = 4)
  colnames(prob_labs) <- c("prob_resting","prob_travelling","prob_foraging","knownState")
  
  predictions = matrix(predict(rf_model, test_data, type = "prob"),ncol=3)
  if(nrow(test_data) > 0){
    for (i in 1:nrow(test_data)) {
      conf_matrix[test_data$knownState[i], ] <- conf_matrix[test_data$knownState[i], ] + predictions[i, ]
      prob_labs <- rbind(prob_labs,c(predictions[i, ],test_data$knownState[i]))
    }
  }
  write.csv(conf_matrix,
            make_title(paste0(directory,"/params/"),
                       paste0("rf-rf-",holdout_whale,"-",
                              "confusion_matrix.csv")))
  write.csv(prob_labs,
            make_title(paste0(directory,"/params/"),
                       paste0("rf-rf-",holdout_whale,"-",
                              "probs_labs.csv")))
  
  # Train SVM model
  svm_model <- svm(
    factor(knownState) ~ logMDDD.x + logMDDD.y, 
    data = train_data,
    kernel = "radial", 
    weights = weights,
    probability = TRUE  # Enable probability predictions
  )
  conf_matrix <- matrix(0, nrow = 3, ncol = 3)
  rownames(conf_matrix) <- c("True Resting", "True Travelling", "True Foraging")
  colnames(conf_matrix) <- c("Predicted Resting", "Predicted Travelling", "Predicted Foraging")

  prob_labs <- matrix(0,ncol = 4)
  colnames(prob_labs) <- c("prob_resting","prob_travelling","prob_foraging","knownState")
  
  if(nrow(test_data) > 0){
    predictions <- predict(svm_model, test_data, decision.values = TRUE, probability = TRUE)
    predictions = matrix(attr(predictions, "probabilities"),ncol=3)
    for (i in 1:nrow(test_data)) {
      conf_matrix[test_data$knownState[i], ] <- conf_matrix[test_data$knownState[i], ] + predictions[i, ]
      prob_labs <- rbind(prob_labs,c(predictions[i, ],test_data$knownState[i]))
    }
  }
  write.csv(conf_matrix,
            make_title(paste0(directory,"/params/"),
                       paste0("svm-svm-",holdout_whale,"-",
                              "confusion_matrix.csv")))
  write.csv(prob_labs,
            make_title(paste0(directory,"/params/"),
                       paste0("svm-svm-",holdout_whale,"-",
                              "probs_labs.csv")))

  # Train logistic regression
  lr_model <- multinom(
    factor(knownState) ~ logMDDD.x + logMDDD.y, 
    data = train_data,
    weights = weights
  )
  conf_matrix <- matrix(0, nrow = 3, ncol = 3)
  rownames(conf_matrix) <- c("True Resting", "True Travelling", "True Foraging")
  colnames(conf_matrix) <- c("Predicted Resting", "Predicted Travelling", "Predicted Foraging")
  
  prob_labs <- matrix(0,ncol = 4)
  colnames(prob_labs) <- c("prob_resting","prob_travelling","prob_foraging","knownState")
  
  if(nrow(test_data) > 0){
    predictions <- matrix(predict(lr_model, test_data, type = "probs"),ncol=3)
    for (i in 1:nrow(test_data)) {
      conf_matrix[test_data$knownState[i], ] <- conf_matrix[test_data$knownState[i], ] + predictions[i, ]
      prob_labs <- rbind(prob_labs,c(predictions[i, ],test_data$knownState[i]))
    }
  }
  write.csv(conf_matrix,
            make_title(paste0(directory,"/params/"),
                       paste0("lr-lr-",holdout_whale,"-",
                              "confusion_matrix.csv")))
  write.csv(prob_labs,
            make_title(paste0(directory,"/params/"),
                       paste0("lr-lr-",holdout_whale,"-",
                              "probs_labs.csv")))
}
