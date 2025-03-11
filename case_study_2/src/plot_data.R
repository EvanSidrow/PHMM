colors <- hcl(h = seq(15, 375, length = 6 + 1),
              l = 65, c = 100)[1:6]

colors <- c("1" = colors[1],
            "2" = colors[2],
            "3" = colors[3],
            "4" = colors[4],
            "5" = colors[5],
            "6" = colors[6])

labs <- c(Elevation = "Depth (m)",
          maxDepth = "Maximum Depth (m)",
          diveDuration = "Dive Duration (s)",
          w_low = "Wiggliness (Low Frequency)",
          w_high = "Wiggliness (High Frequency)",
          logWLow = "Wiggliness (Low Frequency) (log10)",
          logWHigh = "Wiggliness (High Frequency) (log10)",
          logWTotal = "Wiggliness (Total) (log10)",
          postDiveInt = "Post Dive Interval (s)",
          htv = "Heading Total Variation (rad / s)",
          logHtv = "Heading Total Variation (rad) (log10)",
          logJpNorm = "Jerk Peak (normalized, log10)",
          jp_normed = "Jerk Peak (normalized)",
          delt_d = "Change in Depth (m)",
          roll = "Roll (rad)",
          rajp = "Roll @ Peak Jerk (rad)",
          p_catch = "Probability of Catch",
          VeDBA = "VeDBA")

df$knownState <- as.factor(df$knownState)


p1 <- ggplot() +
  geom_density_2d(data = df, 
                  aes(x = htv, y = jp_normed), 
                  color = "black", linewidth = 0.25, alpha = 0.5,
                  bins = 20) +
  geom_point(data = subset(df, knownState != 7), 
             aes(x = htv, y = jp_normed, color = knownState)) +
  theme_classic() + 
  scale_color_manual(limits = c("1","4","5","6"),
                     labels = c("1" = "Descent",
                                "4" = "Capture",
                                "5" = "Ascent w/o fish",
                                "6" = "Ascent w/ fish"),
                     values = colors) + 
  labs(x = "Heading Total Variation (rad/s)",
       y = "Jerk Peak (normalized)",
       color = "") +
  scale_x_log10() +
  scale_y_log10() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=12))

p2 <- ggplot() +
  geom_density_2d(data = df, 
                  aes(x = htv, y = delt_d), 
                  color = "black", linewidth = 0.25, alpha = 0.5,
                  bins = 20) +
  geom_point(data = subset(df, knownState != 7), 
             aes(x = htv, y = delt_d, color = knownState)) +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_color_manual(limits = c("1","4","5","6"),
                     labels = c("1" = "Descent",
                                "4" = "Capture",
                                "5" = "Ascent w/o fish",
                                "6" = "Ascent w/ fish"),
                     values = colors) + 
  labs(x = "Heading Total Variation (rad/s)",
       y = "Change in Depth (m)",
       color = "") +
  scale_x_log10() +
  scale_y_reverse() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=12))

p3 <- ggplot() +
  geom_density_2d(data = df, 
                  aes(x = jp_normed, y = delt_d), 
                  color = "black", linewidth = 0.25, alpha = 0.5,
                  bins = 20) +
  geom_point(data = subset(df, knownState != 7), 
             aes(x = jp_normed, y = delt_d, color = knownState)) +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_color_manual(limits = c("1","4","5","6"),
                     labels = c("1" = "Descent",
                                "4" = "Capture",
                                "5" = "Ascent w/o fish",
                                "6" = "Ascent w/ fish"),
                     values = colors) + 
  labs(x = "Jerk Peak (normalized)",
       y = "Change in Depth (m)",
       color = "") +
  scale_x_log10() +
  scale_y_reverse() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=12))

ggsave("plt/data_htv_jp.png", 
       p1, 
       width = 5.5, height = 4)
ggsave("plt/data_delt_d_htv.png", 
       p2, 
       width = 4, height = 4)
ggsave("plt/data_delt_d_jp.png", 
       p3, 
       width = 4, height = 4)

### plot the dive-level statistics
dive_df <- df

# get max depths

dive_maxs <- dive_df %>% 
  group_by(ID) %>%
  summarize(max_depth = max(ad))
dive_df <- left_join(dive_df,
                            dive_maxs[c("ID","max_depth")])

# get bottom phase
dive_df$bottom <- dive_df$ad > 0.7*dive_df$max_depth
dive_df <- dive_df[dive_df$bottom,]

# get jerk peak and heading total variation
dive_sums <- dive_df %>% 
  group_by(ID) %>%
  summarize(avg_htv_bot = mean(htv),
            jp_normed = max(jp_normed),
            max_depth = max(ad))

# get roll at jerk peak
dive_sums <- left_join(dive_sums,
                       dive_df[c("ID","jp_normed","rajp")])

dive_sums$label <- NA
dive_sums$label[dive_sums$ID %in% pos_dives] <- T
dive_sums$label[dive_sums$ID %in% neg_dives] <- F

# plot data
p1 <- ggplot(dive_sums,aes(x = jp_normed, y = avg_htv_bot, color = label)) +
  geom_point() +
  theme_classic() + 
  labs(x = "Jerk Peak (normalized)",
       y = "Average Bottom Heading total Variation (rad/s)",
       color = "") +
  scale_color_manual(values = c("red","blue","grey"),
                     labels = c("No capture","Capture","Unknown")) +
  scale_x_log10() +
  scale_y_log10() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=12))
ggsave("plt/dive_data_jp_htv.png", width = 5.5, height = 4)

p2 <- ggplot(dive_sums,aes(x = rajp, y = avg_htv_bot, color = label)) +
  geom_point() +
  theme_classic() + 
  labs(x = "Roll at Jerk Peak",
       y = "Average Bottom Heading total Variation (rad/s)",
       color = "Prey Capture") +
  scale_color_manual(values = c("red","blue","grey"),
                     labels = c("No capture","Capture","Unknown")) +
  scale_y_log10() +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=12))
ggsave("plt/dive_data_rajp_htv.png", width = 4, height = 4)

p3 <- ggplot(dive_sums,aes(x = rajp, y = jp_normed, color = label)) +
  geom_point() +
  theme_classic() + 
  labs(x = "Roll at Jerk Peak",
       y = "Jerk Peak (normalized)",
       color = "") + 
  scale_color_manual(values = c("red","blue","grey"),
                     labels = c("No capture","Capture","Unknown")) +
  scale_y_log10() +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        text = element_text(size=12))
ggsave("plt/dive_data_rajp_jp.png", width = 4, height = 4)
