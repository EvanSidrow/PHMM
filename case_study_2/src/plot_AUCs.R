# make plots of AUCs
means <- c()
sds <- c()
AUCs <- c()

param_directory <- "params"#"params_unshared_label_1"
plt_directory <- "plt"#"plt_unshared_label_1"

# load in base AUCS
filename <- paste0(param_directory,"/AUC_base_",K,".csv")
AUCs <- c(AUCs,data.frame(fread(filename))[,2])
means <- c(means,mean(data.frame(fread(filename))[,2]))
sds <- c(sds,sd(data.frame(fread(filename))[,2]))

# load in lr
filename <- paste0(param_directory,"/AUC_lr_",K,".csv")
AUCs <- c(AUCs,data.frame(fread(filename))[,2])
means <- c(means,mean(data.frame(fread(filename))[,2]))
sds <- c(sds,sd(data.frame(fread(filename))[,2]))

# load in svm
filename <- paste0(param_directory,"/AUC_svm_",K,".csv")
AUCs <- c(AUCs,data.frame(fread(filename))[,2])
means <- c(means,mean(data.frame(fread(filename))[,2]))
sds <- c(sds,sd(data.frame(fread(filename))[,2]))

# load in rf
filename <- paste0(param_directory,"/AUC_rf_",K,".csv")
AUCs <- c(AUCs,data.frame(fread(filename))[,2])
means <- c(means,mean(data.frame(fread(filename))[,2]))
sds <- c(sds,sd(data.frame(fread(filename))[,2]))

# load in other AUCs
for(log_alpha in c("-4","-3","-2","-1","0")){
  filename <- paste0(param_directory,"/AUC_delt_d_htv_jp_normed_",log_alpha,"_",K,".csv")
  AUCs <- c(AUCs,data.frame(fread(filename))[,2])
  means <- c(means,mean(data.frame(fread(filename))[,2]))
  sds <- c(sds,sd(data.frame(fread(filename))[,2]))
}

# make dataframe
df <- data.frame(alpha = rep(c("Baseline","LR","SVM","RF","0.0001","0.001","0.01","0.1","1"),each = K),
                 baseline = rep(c(1,2,3,4,F,F,F,F,F),each = K),
                 fold = rep(1:K,9),
                 AUC = AUCs)

df_means <- data.frame(alpha = c("Baseline","LR","SVM","RF","0.0001","0.001","0.01","0.1","1"),
                       baseline = c(1,2,3,4,F,F,F,F,F),
                       mean = means,
                       sd = sds)
# ggplot
plot0 <- ggplot(df, aes(x = alpha, y = AUC)) + 
  geom_jitter(size = 4, alpha = 0.5,
             width = 0.05, height = 0.00) +
  geom_point(data = df_means, 
             aes(x = alpha, y = mean),
             shape = "*", size = 20) + 
  geom_line(data = df_means, 
             aes(x = alpha, y = mean, group = baseline)) + 
  geom_vline(aes(xintercept = 5.5)) +
  labs(x=TeX("\\alpha"),y="AUC") +
  theme_classic() +
  theme(text = element_text(size=16))

ggsave(paste0(directory,"/plt/AUCs_",
              paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
              round(log10(lambda),3),"_",
              K,".png"),
       plot0,
       width = 6, height = 4)
