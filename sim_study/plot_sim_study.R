library(momentuHMM)
library(ggplot2)
library(GGally)
library(latex2exp)
library(pROC)
library(tidyr)
library(purrr)
library(markovchain)

setwd("/Users/evsi8432/Documents/Research/PHMM/sim_study")

feat_to_change = "g"
feats = c(0.75, 0.95, 0.99)

dof <- 4 # 2.5, 4, 100
sep <- 1 # 0.5, 1, 2
g <- 0.95 # 0.75, 0.95, 0.99
l <- 0.01 # 0.001, 0.01, 0.1 prop of labels
N <- 2000 # 200, 2000, 20000

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

# load old data
load(file = paste0("results/",feat_to_change,".Rdata"))

# plot histograms
hist_df <- data.frame(feat = c(),
                      alpha = c(), 
                      X = c(), 
                      Y = c())

density_df <- data.frame(feat = c(),
                         alpha = c(),
                         X = c(),
                         Y = c(), 
                         density = c())

nbins = 100

for(feat in feats){
  feat_str <- as.character(feat)
  assign(feat_to_change, feat)
  hist_data <- sample_data(N,l,sep,dof,g)
  hist_df_feat <- data.frame(feat = rep(feat,3*N),
                             alpha = rep(c(0.0,l,1.0),each=N), 
                             X = factor(rep(hist_data$X,3)), 
                             Y = rep(hist_data$Y1,3))
  hist_df <- rbind(hist_df,hist_df_feat)
  
  Y_vals <- seq(min(hist_df_feat$Y), max(hist_df_feat$Y), length.out = 250)
  
  for(lamb in c(0.0,l,1.0)){
    lamb_str <- as.character(lamb)
    hmm <- hmms[[feat_str]][[lamb_str]][[1]]
    
    w <- steadyStates(new("markovchain", transitionMatrix = hmm$mle$gamma))
    mean1 = hmm$mle$Y1[1,1]
    sd1 = hmm$mle$Y1[2,1]
    mean2 = hmm$mle$Y1[1,2]
    sd2 = hmm$mle$Y1[2,2]
    
    density_df <- rbind(density_df, 
                        data.frame(
                          feat = feat_str,
                          alpha = lamb_str,
                          X = factor(1),
                          Y = Y_vals,
                          density = dnorm(Y_vals, mean = mean1, sd = sd1) * ((max(Y_vals)-min(Y_vals))/nbins)*w[1]*N
                        ),
                        data.frame(
                          feat = feat_str,
                          alpha = lamb_str,
                          X = factor(2),
                          Y = Y_vals,
                          density = dnorm(Y_vals, mean = mean2, sd = sd2) * ((max(Y_vals)-min(Y_vals))/nbins)*w[2]*N
                        ))
  }
}

hist_df$feat = factor(hist_df$feat, levels=feats)
hist_df$alpha = factor(hist_df$alpha, levels=c(0,l,1))
density_df$feat = factor(density_df$feat, levels=feats)
density_df$alpha = factor(density_df$alpha, levels=c(0,l,1))

alpha_labels = setNames(paste("α =",c(0,l,1)), c(0,l,1))
feat_labels = setNames(paste(feat_to_change,"=",feats), feats)

# plot the data itself
plot_data <- sample_data(2000,0.01,1,4,0.95)
plot_data$t <- 1:N
plot_data$X <- factor(plot_data$X)
p <- ggplot(plot_data, aes(x = t, y = Y1)) +
  geom_line(color="black",alpha=0.1) + 
  geom_point(aes(color = X),size = 0.5) +
  labs(title = "", x = TeX("$t$"), y = TeX("$Y_t$"), 
       color = TeX("$X_t$")) + 
  theme_classic() + 
  geom_vline(data = plot_data[!is.na(plot_data$Z),],
             aes(xintercept = t, color = X), 
             linetype = "dashed", size = 0.25, alpha = 0.75,
             show.legend = FALSE)
print(p)

ggsave(paste0("plt/sim_data.tiff"),
       width = 6,
       height = 3,
       units = "in",
       dpi=600)

ggsave(paste0("plt/sim_data.png"),
       width = 6,
       height = 3,
       units = "in",
       dpi=600)

# Create histograms of the state-dependent distributions data
p <- ggplot(hist_df, aes(x = Y, fill = X)) +
  geom_histogram(bins = nbins, color = "black", alpha = 0.4, size = 0.25, position = "identity") +
  geom_line(data = density_df, aes(x = Y, y = density, color = X)) +  # Theoretical normal densities
  labs(title = "", x = TeX("$Y_t$"), y = "Count", 
       fill = TeX("True $X_t$"), color = TeX("Estimated $X_t$")) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=16),
        plot.title = element_text(hjust = 0.5)) + 
  ggh4x::facet_grid2(feat ~ alpha, scales = "free", independent = "x",
                     labeller = labeller(
                       alpha = as_labeller(alpha_labels),
                       feat = as_labeller(feat_labels)
                     )) 

print(p)

ggsave(paste0("plt/hists_",feat_to_change,".tiff"),
       width = 8,
       height = 6,
       units = "in",
       dpi=300)

plot_data <- data.frame("AUC" = AUCs,
                        "Sensitivity" = sens,
                        "Specificity" = spec,
                        feat_to_change = feat_vec, 
                        "alpha" = lamb_vec)

plot_data <- plot_data %>%
  pivot_longer(cols = c("AUC","Sensitivity","Specificity"), 
               names_to = "Variable", 
               values_to = "Value")

plot_data$feat_to_change = factor(plot_data$feat_to_change, levels=feats)

# Create the box-and-whiskers plots using facet_grid to structure by two variables
plot0 <- ggplot(plot_data, aes(x = alpha, y = Value, fill = alpha)) +
  geom_boxplot(alpha = 0.6) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank()) +
  scale_fill_grey(start = 0.8, end = 0.2) +
  labs(title = "", x = "", y = "") +
  ggh4x::facet_grid2(Variable ~  feat_to_change, scales = "free_y", independent = "y",
                     labeller = labeller(
                       feat_to_change = as_labeller(feat_labels)
                     ))

print(plot0)
ggsave(paste0("plt/boxplots_",feat_to_change,".tiff"),
       width = 8,
       height = 6,
       units = "in",
       dpi=300)