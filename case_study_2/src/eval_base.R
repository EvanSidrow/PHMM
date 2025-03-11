Data_for_model <- df[df$ID %in% test_dives & 
                     df$ID %in% c(pos_dives,neg_dives),]

# get max depths
dive_maxs <- Data_for_model %>% 
  group_by(ID) %>%
  summarize(max_depth = max(ad),
            dive_duration = max(stime)-min(stime))

Data_for_model <- left_join(Data_for_model,
                            dive_maxs[c("ID","max_depth","dive_duration")])

# get bottom phase
Data_for_model$bottom <- Data_for_model$ad > 0.7*Data_for_model$max_depth
Data_for_model <- Data_for_model[Data_for_model$bottom,]

# get jerk peak and heading total variation
dive_sums <- Data_for_model %>% 
  group_by(ID) %>%
  summarize(avg_htv_bot = mean(htv),
            jp_normed = max(jp_normed),
            max_depth = max(ad),
            dive_duration = mean(dive_duration))

# get roll at jerk peak
dive_sums <- left_join(dive_sums,
                       Data_for_model[c("ID","jp_normed","rajp")])

dive_sums$label <- dive_sums$ID %in% pos_dives

probs <- c()
labs <- c()

conf_matrices_base[[k]] = matrix(0, nrow = 2, ncol = 2)
rownames(conf_matrices_base[[k]]) <- c("True Positive", "True Negative")
colnames(conf_matrices_base[[k]]) <- c("Pred Positive", "Pred Negative")

for(ID in test_dives){
  
  if(ID %in% pos_dives){
    rownum <- 1
    labs <- c(labs,TRUE)
  } else if (ID %in% neg_dives){
    rownum <- 2
    labs <- c(labs,FALSE)
  }
  
  htv_check <- dive_sums$avg_htv_bot[dive_sums$ID == ID] >= base_model$htv
  jp_check <- dive_sums$jp_normed[dive_sums$ID == ID] >= base_model$jp_normed
  rajp_check <- abs(dive_sums$rajp[dive_sums$ID == ID]) >= base_model$rajp
  
if(htv_check & jp_check & rajp_check){
    colnum <- 1
    probs <- c(probs,1)
  } else {
    colnum <- 2
    probs <- c(probs,0)
  }
  
  conf_matrices_base[[k]][rownum,colnum] = conf_matrices_base[[k]][rownum,colnum] + 1
}

probs_base[[k]] <- probs
AUCs_base[k] <- roc(response=labs, predictor=probs, direction = "<")$auc
plot(roc(response = labs, predictor=probs))

# do prediction for random forest
print("rf")
probs = predict(rf_model, dive_sums, type = "prob")[,2]
probs_rf[[k]] <- probs
AUCs_rf[k] <- roc(response=dive_sums$label, predictor=probs, direction = "<")$auc
plot(roc(response = dive_sums$label, predictor=probs))

# do prediction for svm
print("svm")
probs <- predict(svm_model, dive_sums, decision.values = TRUE, probability = TRUE)
probs = attr(probs, "probabilities")[,2]
probs_svm[[k]] <- probs
AUCs_svm[k] <- roc(response = dive_sums$label, predictor=probs, direction = "<")$auc
plot(roc(response = dive_sums$label, predictor=probs))

# do prediction for logistic regression
print("lr")
probs = predict(lr_model, dive_sums, type = "response")
probs_lr[[k]] <- probs
AUCs_lr[k] <- roc(response = dive_sums$label, predictor=probs, direction = "<")$auc
plot(roc(response = dive_sums$label, predictor=probs))

print(dive_sums)
