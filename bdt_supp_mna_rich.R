## 
# supplementary figure
##

library(tidyverse)

# read data
mammal_community_df <- read_rds("processed_data/mammal_community_df_2025-10-16.rds") 

## filter and select
mammal_mna_rich <- mammal_community_df %>%
  mutate(region = case_when(domain_id %in% c("D01", "D02") ~ "Northeast",
                            domain_id == "D05" ~ "Upper Midwest")) %>%
  select(plot_session, region, richness, total_mna)

## plot
ggplot(mammal_mna_rich, aes(x = factor(richness), y = total_mna)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.2, color = "lightgray") +
  geom_boxplot(fill = NA, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_grid(region ~ .) +
  labs(x = "Species Richness", y = "Minimum Number Alive (All Species)") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16))
