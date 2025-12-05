###
# compare species presence to richness
###

## loading packages
library(tidyverse)


## read mammal community data
mammal_community_df <- read_rds("processed_data/mammal_community_df_2025-10-16.rds")

## mammals of interest
mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## nlcd to exclude
developed <- c("cultivatedCrops", "pastureHay")

## filter mammal community data for nlcds to exclude
mamms_presence_df <- mammal_community_df %>%
  filter(!nlcd_class %in% developed) %>%
  ungroup() %>%
  mutate(site_id = str_sub(plot_id, 1, 4),
         plot_id = factor(plot_id))

## generate species pools by nlcd and site
site_nlcd_species_pools <- mamms_presence_df %>%
  filter(presence == 1) %>%
  distinct(site_id, nlcd_class, taxon_id) %>%
  group_by(site_id, nlcd_class) %>%
  summarize(species_pool = list(taxon_id),
            n_species = n(),
            .groups = "drop")

## create a function to randomly sample n species from the pool of species
## observed at a site and land cover class combination
sample_species <- function(plot, reps = 1000){
  ## grab site id and nlcd class for the plot
  site <- str_sub(string = plot, start = 1, end = 4) ## first four characters of plot_id
  nlcd <- mamms_presence_df %>%
    filter(plot_id == plot) %>%
    select(nlcd_class) %>%
    distinct() %>%
    pull(nlcd_class) 
  
  ## get the species pool for the site and nlcd class
  species_pool <- site_nlcd_species_pools %>%
    filter(site_id == site, nlcd_class == nlcd) %>%
    pull(species_pool) %>%
    unlist()
  
  ## get the observed richnesses for the plot's sessions
  session_richness <- mamms_presence_df %>%
    filter(plot_id == plot, richness > 0) %>%
    select(domain_id, plot_session, richness, mean_yday, mean_month) %>%
    distinct() %>%
    arrange(plot_session, richness)
  
  ## create a data frame of the species samples
  sample_comms <- expand_grid(session_richness, rep = 1:reps) %>%
    mutate(species_sample = map(richness, ~sample(species_pool, size = .x, replace = FALSE))) 
  
  ## create a data frame of the species samples with presence information
  out <- expand_grid(session_richness, taxon_id = species_pool) %>%
    uncount(weights = reps) %>%
    mutate(rep = rep(1:reps, times = nrow(session_richness) * length(species_pool))) %>%
    arrange(rep, taxon_id) %>%
    left_join(sample_comms, by = c("domain_id", "rep", "plot_session", "richness", "mean_yday", "mean_month")) %>%
    group_by(rep, plot_session) %>%
    mutate(presence = as.numeric(taxon_id %in% unlist(species_sample)),
           plot_id = plot,
           nlcd = nlcd)
  
  return(out)
}

test <- sample_species(plot = "BART_001", rep = 2)

## simulate communities for all plots
simulated_comms <- map(levels(mamms_presence_df$plot_id), sample_species, reps = 1000, .progress = TRUE)

## and bind the list of data frames into one data frame
simulated_comms_df <- bind_rows(simulated_comms)
rm(simulated_comms)

## filter for mammals of interest
simulated_mamms_df <- simulated_comms_df %>%
  filter(taxon_id %in% mamms)

## fit logistic regression models to simulated data
logistic_sim_results <- simulated_mamms_df %>%
  filter(domain_id %in% c("D05")) %>%
  group_by(taxon_id, rep) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(presence ~ richness, data = ., family = binomial)),
    intercept = map_dbl(model, ~ coef(.x)[["(Intercept)"]]),
    coef  = map_dbl(model, ~ coef(.x)[["richness"]]),
    type = "sim",
    preds = map(model, ~ {
      richness_seq <- seq(0, max(mamms_presence_filtered_df$richness), by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    })
  ) %>%
  select(taxon_id, rep, coef, intercept, type, preds)

## fit logistic regression models to observed data, only for rows where taxon is in species pool
mamms_presence_filtered_df <- mamms_presence_df %>%
  filter(taxon_id %in% mamms, richness > 0) %>%
  left_join(site_nlcd_species_pools, by = c("site_id", "nlcd_class")) %>%
  rowwise() %>%
  filter(taxon_id %in% unlist(species_pool)) %>%
  ungroup() %>%
  select(-species_pool, -n_species)
logistic_obs_results <- mamms_presence_filtered_df %>%
  filter(domain_id %in% c("D05")) %>%
  group_by(taxon_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(presence ~ richness, data = ., family = binomial)),
    intercept = map_dbl(model, ~ coef(.x)[["(Intercept)"]]),
    coef  = map_dbl(model, ~ coef(.x)[["richness"]]),
    type = "obs",
    preds = map(model, ~ {
      richness_seq <- seq(0, max(mamms_presence_filtered_df$richness), by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    })
  ) %>%
  select(taxon_id, coef, intercept, type, preds)

## plot observed preds
logistic_obs_preds <- logistic_obs_results %>%
  unnest(preds)

ggplot(logistic_obs_preds, aes(x = richness, y = fit, color = taxon_id)) +
  geom_line() +
  labs(x = "Species Richness", y = "Probability of Presence") +
  scale_x_continuous(breaks = seq(0, max(mamms_presence_filtered_df$richness), by = 1)) +
  theme_bw()

# plot sim preds
logistic_sim_preds <- logistic_sim_results %>%
  unnest(preds)

ggplot() +
  geom_line(data = logistic_sim_preds, aes(x = richness, y = fit, color = taxon_id, group = interaction(taxon_id, rep)), alpha = 0.01) +
  geom_line(data = logistic_obs_preds, aes(x = richness, y = fit, color = taxon_id)) +
  facet_wrap(~ taxon_id) +
  labs(x = "Species Richness", y = "Probability of Presence") +
  scale_x_continuous(breaks = seq(1, max(mamms_presence_filtered_df$richness), by = 1),
                     limits = c(1, max(mamms_presence_filtered_df$richness))) +
  theme_bw()

mamms_pres_compare <- mamms_presence_filtered_df %>%
  filter(taxon_id %in% "PELE") %>%
  distinct(plot_session, taxon_id)

mna_by_sess_compare <- mna_by_session_corrected %>%
  filter(taxon_id %in% "PELE") %>%
  distinct(plot_session, taxon_id)

plot_sessions_missing <- anti_join(mamms_pres_compare, mna_by_sess_compare, by = c("plot_session", "taxon_id"))  
