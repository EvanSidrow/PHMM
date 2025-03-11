library(momentuHMM)
library(ggplot2)
library(GGally)
library(latex2exp)
library(pROC)
library(tidyr)
library(purrr)
library(markovchain)

setwd("/Users/evsi8432/Documents/Research/PHMM/sim_study")

feat_to_change = commandArgs(trailingOnly=TRUE)

dof <- 4 # 2.5, 4, 100
sep <- 1 # 0.5, 1, 2
g <- 0.95 # 0.9, 0.95, 0.99
l <- 0.01 # 0.001, 0.01, 0.1 prop of labels
N <- 2000 # 200, 2000, 20000
n_exp <- 100

if(feat_to_change == "dof"){
  feats = c(2.5, 4, 100)
} else if(feat_to_change == "sep"){
  feats = c(0.5, 1, 2)
} else if(feat_to_change == "g"){
  feats = c(0.75,0.95,0.99)
} else if(feat_to_change == "l"){
  feats = c(0.001, 0.01, 0.1)
} else if(feat_to_change == "N"){
  feats = c(200, 2000, 20000)
} else{
  print("unknown feature to change")
  stop()
}

# sample data
sample_data <- function(N,l,sep,dof,g){
  
  # define delta and Gamma
  Gamma <- matrix(c(0.99, 0.01,
                    1-g, g), 
                  byrow=T,nrow = 2, ncol = 2)
  delta <- steadyStates(new("markovchain", transitionMatrix = Gamma))
  
  # set mean and sd
  sd <- sqrt(dof/(dof-2))
  mu <- c(-sep*sd/2,sep*sd/2)
  
  X <- rep(1,N)
  Y1 <- rt(n=N,df=dof)
  Y2 <- rt(n=N,df=dof)
  
  X[1] <- 2#sample(c(1,2), size=1, prob = delta)
  Y1[1] <- Y1[1] + mu[X[1]]
  Y2[1] <- Y2[1] + mu[X[1]]
    
  for (i in 2:N){
    X[i] <- sample(c(1,2), size=1, prob = Gamma[X[i-1],])
    Y1[i] <- Y1[i] + mu[X[i]]
    Y2[i] <- Y2[i] + mu[X[i]]
  }
  
  Z <- rep(NA,N)
  inds <- c(sample(which(X == 1), size = min(sum(X == 1),
                                             max(2,round(delta[1]*l*N)))),
            sample(which(X == 2), size = min(sum(X == 2),
                                             max(2,round(delta[2]*l*N)))))
  Z[inds] <- X[inds]
  return(data.frame("X"=X,"Y1"=Y1,"Y2"=Y2,"Z"=Z))
}

hmms <- list()
conf_matrices <- list()
AUCs <- list()
sens <- list()
spec <- list()

# initialize results
for(feat in feats){
  feat_str <- as.character(feat)
  hmms[[feat_str]] <- list()
  conf_matrices[[feat_str]] <- list()
  for(lamb in c(0.0,0.001,0.01,0.1,0.5,1.0)){
    lamb_str <- as.character(lamb)
    hmms[[feat_str]][[lamb_str]] <- list()
    conf_matrices[[feat_str]][[lamb_str]] <- list()
  }
}

total_exps <- length(feats)*5*length(n_exp)

AUCs <- rep(NA,total_exps)
sens <- rep(NA,total_exps)
spec <- rep(NA,total_exps)
feat_vec <- rep(NA,total_exps)
lamb_vec <- rep(NA,total_exps)
n <- 1
prog_exp <- 0

for(feat in feats){
  print(feat)
  assign(feat_to_change, feat)
  feat_str <- as.character(feat)
  for(exp in 1:n_exp){
    
    print(prog_exp/(3*n_exp))
    prog_exp <- prog_exp + 1
    
    exp_str <- as.character(exp)
    
    # sample data
    train_data <- sample_data(N,l,sep,dof,g)
    test_data <- sample_data(N,l,sep,dof,g)
    
    # fit PHMM with multiple values of alpha
    sd <- sqrt(dof/(dof-2))
    mu <- c(-sep*sd/2,sep*sd/2)
    
    Par0 <- list()
    Par0[["Y1"]] <- matrix(rep(NA,4), nrow = 2)
    Par0[["Y1"]][,1] <- mu
    Par0[["Y1"]][,2] <- sd
    
    for(lamb in c(0.0,0.001,0.01,0.1,0.5,1.0)){
      lamb_str <- as.character(lamb)
      
      hmms[[feat_str]][[lamb_str]][[exp_str]] <- tryCatch(
        {
          fitHMM(data=prepData(train_data,coordNames=NULL),
                                                nbStates=2,
                                                Par0 = Par0,
                                                dist=list("Y1"="norm"),
                                                knownStates=train_data$Z,
                                                lambda=lamb,
                                                retryFits = 0)
        },
        error = function(e){ 
          return(NA)
        }
      )
      
      if(typeof(hmms[[feat_str]][[lamb_str]][[exp_str]]) == "logical"){
        feat_vec[n] = feat_str
        lamb_vec[n] = lamb_str
        n <- n+1
        next
      }
      
      # run test data without labels
      hmm_pars <- getPar0(hmms[[feat_str]][[lamb_str]][[exp_str]])
      eps <- 1e-4
      
      hmm0 <- tryCatch(
        {
          fitHMM(data=prepData(test_data,coordNames=NULL),
                     nbStates=2,
                     Par0=hmm_pars$Par,
                     beta0=hmm_pars$beta,
                     delta0=(1-eps)*(hmm_pars$delta)+eps*rep(1/2,2),
                     dist=list("Y1"="norm"),#,"Y2"="norm"),
                     nlmPar = list('stepmax'=1e-100,
                                   'iterlim'=1))
        },
        error = function(e){ 
          return(NA)
        }
      )
      
      # skip if the test set is too large
      if(typeof(hmm0) == "logical"){
        conf_matrices[[feat_str]][[lamb_str]][[exp_str]] <- matrix(c(NA,NA,NA,NA),nrow=2,ncol=2)
        feat_vec[n] = feat_str
        lamb_vec[n] = lamb_str
        n <- n+1
        next
      }
      if(hmm0$mod$minimum > 1e+308){
        conf_matrices[[feat_str]][[lamb_str]][[exp_str]] <- matrix(c(NA,NA,NA,NA),nrow=2,ncol=2)
        feat_vec[n] = feat_str
        lamb_vec[n] = lamb_str
        n <- n+1
        next
      }
      
      # evaluate model
      probs0 <- stateProbs(hmm0)
      
      conf_matrices[[feat_str]][[lamb_str]][[exp_str]] <- matrix(c(0,0,0,0),nrow=2,ncol=2)
      for(i in 1:2){
        for(j in 1:2){
          conf_matrices[[feat_str]][[lamb_str]][[exp_str]][i,j] = sum(probs0[test_data$X == i,j])
        }
      }
      
      rownames(conf_matrices[[feat_str]][[lamb_str]][[exp_str]]) <- c("True 1", "True 2")
      colnames(conf_matrices[[feat_str]][[lamb_str]][[exp_str]]) <- c("Predicted 1", "Predicted 2")
      
      # record metrics
      AUCs[n] = roc(response  = test_data$X == 1, 
                    predictor = probs0[,1] ,
                    direction = "<")$auc
      spec[n] = conf_matrices[[feat_str]][[lamb_str]][[exp_str]][1,1] / sum(conf_matrices[[feat_str]][[lamb_str]][[exp_str]][,1])
      sens[n] = conf_matrices[[feat_str]][[lamb_str]][[exp_str]][2,2] / sum(conf_matrices[[feat_str]][[lamb_str]][[exp_str]][,2])
      feat_vec[n] = feat_str
      lamb_vec[n] = lamb_str
      n <- n+1
    }
  }
}

save(hmms,conf_matrices,AUCs,spec,sens,feat_vec,lamb_vec,
     file = paste0("results/",feat_to_change,".Rdata"))