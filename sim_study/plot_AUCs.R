library(momentuHMM)
library(ggplot2)
library(GGally)
library(latex2exp)
library(pROC)
library(tidyr)
library(purrr)
library(markovchain)

setwd("/Users/evsi8432/Documents/Research/PHMM/sim_study")

# Define parameter ranges
variables <- list(
  "dof" = c(2.5, 4, 100),
  "sep" = c(0.5, 1, 2),
  "g" = c(0.75, 0.95, 0.99),
  "l" = c(0.001, 0.01, 0.1),
  "N" = c(200, 2000, 20000)
)

# Data storage
plot_data <- data.frame()
value_levels <- c()

for (var_name in names(variables)) {
  
  print(var_name)
  load(file = paste0("results/",var_name,".Rdata"))
  plot_data <- plot_data <- rbind(plot_data,
                                  data.frame("AUC" = AUCs,
                                             "Sensitivity" = sens,
                                             "Specificity" = spec,
                                             "FeatureChanged" = var_name,
                                             "FeatureValue" = feat_vec,
                                             "alpha" = lamb_vec))
  value_levels <- c(value_levels,unique(feat_vec))
}

plot_data$FeatureChanged <- factor(plot_data$FeatureChanged, levels = names(variables))
levels(plot_data$FeatureChanged) = c("dof","sep","γ","ℓ","T")

plot_data$FeatureValue <- factor(plot_data$FeatureValue, levels = value_levels)
plot_data$logit_AUC <- -log10(1.0 - plot_data$AUC)

#plot_data <- plot_data %>%
#  pivot_longer(cols = c("logit_AUC","AUC","Sensitivity","Specificity"), 
#               names_to = "Variable", 
#               values_to = "Value")

plot_data$experiment <- interaction(plot_data$FeatureChanged, 
                                    plot_data$FeatureValue, sep = " = ")

inv_logit <- function(x) {
  return(1 - 10^(-x))
}

# Create the box-and-whiskers plots using facet_grid to structure by two variables
plot0 <- ggplot(plot_data, 
                aes(x = alpha, y = logit_AUC, fill = alpha)) +
  geom_boxplot(alpha = 0.6) +
  #theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank()) +
  labs(title = "", x = "", y = "", fill="α") +
  scale_y_continuous(
    breaks = c(0,1,2,3,4),  # Set breaks for logit scale
    labels = function(x) scales::label_number(accuracy = 0.0001)(inv_logit(x))  # Apply inverse logit to labels
  ) +
  facet_wrap(~ experiment, 
             scales = "free_x",
             nrow = 5, ncol = 3)

print(plot0)
ggsave(paste0("plt/sim_boxplots_log.tiff"),
       width = 8,
       height = 10,
       units = "in",
       dpi=600)

ggsave(paste0("plt/sim_boxplots_log.png"),
       width = 8,
       height = 10,
       units = "in",
       dpi=600)

plot0 <- ggplot(plot_data, 
                aes(x = alpha, y = AUC, fill = alpha)) +
  geom_boxplot(alpha = 0.6) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=16),
        plot.title = element_text(hjust = 0.5),
        axis.text.x=element_blank()) +
  labs(title = "", x = "", y = "", fill="α") +
  #scale_fill_grey(start = 1.0, end = 0.0) +
  facet_wrap(~ interaction(FeatureChanged, FeatureValue, sep = " = "), 
             scales = "free_x",
             nrow = 5, ncol = 3) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5, alpha = 0.5) + 
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.5) + 
  geom_hline(yintercept = 01.0, linetype = "dashed", color = "black", size = 0.5, alpha = 0.5)

print(plot0)
ggsave(paste0("plt/sim_boxplots.tiff"),
       width = 8,
       height = 10,
       units = "in",
       dpi=600)

ggsave(paste0("plt/sim_boxplots.png"),
       width = 8,
       height = 10,
       units = "in",
       dpi=600)
