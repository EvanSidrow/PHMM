Data_for_model <- df[df$ID %in% train_dives & 
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

# summarize results in list (this defines the model)
base_model <- list("htv" = min(dive_sums$avg_htv_bot[dive_sums$label]),
                   "jp_normed" = min(dive_sums$jp_normed[dive_sums$label]),
                   "rajp" = min(abs(dive_sums$rajp[dive_sums$label])))

# Train Random Forest model
rf_model <- randomForest(
  factor(label) ~ avg_htv_bot + jp_normed + rajp, 
  data = dive_sums,
  ntree=500, 
  importance = TRUE,
)

# Train SVM model
svm_model <- svm(
  factor(label) ~ avg_htv_bot + jp_normed + rajp, 
  data = dive_sums,
  kernel = "radial",  
  probability = TRUE  # Enable probability predictions
)

# Train logistic regression
lr_model <- logistf(factor(label) ~ avg_htv_bot + jp_normed + rajp, 
                    data = dive_sums)
                    