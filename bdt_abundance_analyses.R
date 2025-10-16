###
# compare species abundance/proportion to richness and overall abundance
###

## loading packages
library(tidyverse)


## read mammal community data
mammal_community_df <- read_rds("processed_data/mammal_community_df_2025-10-16.rds")

## mammals of interest
mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## for each species, what proportion of the community does it make up when it is present?
## binomial regression on richness and total_mna
species_prop_obs <- mammal_community_df %>%
  filter(taxon_id %in% mamms, mna > 0) %>%
  group_by(taxon_id) %>%
  nest() %>%
  mutate(
    rich_mod = map(data, ~ glm(prop_mna ~ richness, data = ., weights = total_mna, family = "binomial")),
    abund_mod = map(data, ~ glm(prop_mna ~ total_mna, data = ., weights = total_mna, family = "binomial")),
    rich_coef = map_dbl(rich_mod, ~ coef(.)[2]),
    abund_coef = map_dbl(abund_mod, ~ coef(.)[2]),
    rich_p = map_dbl(rich_mod, ~ summary(.)$coefficients[2,4]),
    abund_p = map_dbl(abund_mod, ~ summary(.)$coefficients[2,4]),
    rich_preds = map(rich_mod, ~ {
      richness_seq <- seq(0, max(mammal_community_df$richness), by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    }),
    abund_preds = map(abund_mod, ~ {
      abund_seq <- seq(1, max(mammal_community_df$total_mna), by = 1)
      pred_df <- tibble(total_mna = abund_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    })
  ) %>%
  select(-data, -rich_mod, -abund_mod)

rich_preds_obs <- species_prop_obs %>%
  select(taxon_id, rich_preds) %>%
  unnest(cols = c(rich_preds))

abund_preds_obs <- species_prop_obs %>%
  select(taxon_id, abund_preds) %>%
  unnest(cols = c(abund_preds))

## plot richness predictions by taxon
ggplot(rich_preds_obs, aes(x = richness, y = fit, color = taxon_id)) +
  geom_line(size = 1) +
  geom_line(aes(x = richness, y = 1/richness), linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = seq(1, max(rich_preds_obs$richness), by = 1),
                     limits = c(1, max(rich_preds_obs$richness))) +
  labs(x = "Species Richness", y = "Predicted Proportion of Community") +
  theme_bw()

## plot abundance predictions by taxon
ggplot(abund_preds_obs, aes(x = total_mna, y = fit, color = taxon_id)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(0, max(abund_preds_obs$total_mna), by = 5)) +
  labs(x = "Total MNA", y = "Predicted Proportion of Community") +
  theme_bw()


## now we need to bootstrap random communities to compare the observed model results to...
