###
# compare species presence to richness
###

## loading packages
library(tidyverse)


## read mammal community data
mammal_community_df <- read_rds("processed_data/mammal_community_df_2025-10-16.rds") %>%
  mutate(region = if_else(domain_id == "D05", "Upper Midwest", "Northeast")) %>%
  group_by(region) %>%
  mutate(max_richness = max(richness),
         max_mna = max(total_mna)) %>%
  ungroup()

## mammals of interest
mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## nlcd to exclude
developed <- c("cultivatedCrops", "pastureHay")

plots <- mammal_community_df %>%
  distinct(plot_id, nlcd_class) %>%
  arrange(nlcd_class, plot_id)

years <- mammal_community_df %>%
  distinct(year) %>%
  arrange(year)

months <- mammal_community_df %>%
  distinct(mean_month) %>%
  arrange(mean_month)

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
    select(region, domain_id, plot_session, richness, max_richness, mean_yday, mean_month) %>%
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
    left_join(sample_comms, by = c("region", "domain_id", "rep", "plot_session", "richness", "max_richness", "mean_yday", "mean_month")) %>%
    group_by(rep, plot_session) %>%
    mutate(presence = as.numeric(taxon_id %in% unlist(species_sample)),
           plot_id = plot,
           nlcd = nlcd) %>%
    ungroup()
  
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
rm(simulated_comms_df)

mamms_presence_filtered_df <- mamms_presence_df %>%
  filter(taxon_id %in% mamms, richness > 0) %>%
  left_join(site_nlcd_species_pools, by = c("site_id", "nlcd_class")) %>%
  rowwise() %>%
  filter(taxon_id %in% unlist(species_pool)) %>%
  ungroup() %>%
  select(-species_pool, -n_species)

## fit logistic regression models to simulated data
logistic_sim_results <- simulated_mamms_df %>%
  group_by(region, taxon_id, rep, max_richness) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(presence ~ I(1/richness), data = ., family = binomial)),
    # intercept = map_dbl(model, ~ coef(.x)[["(Intercept)"]]),
    # coef  = map_dbl(model, ~ coef(.x)[["richness"]]),
    type = "sim",
    preds = map(model, ~ {
      richness_seq <- seq(0, max_richness, by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    })
  )

## fit logistic regression models to observed data, only for rows where taxon is in species pool

logistic_obs_results <- mamms_presence_filtered_df %>%
  group_by(region, taxon_id, max_richness) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(presence ~ I(1/richness), data = ., family = binomial)),
    # intercept = map_dbl(model, ~ coef(.x)[["(Intercept)"]]),
    # coef  = map_dbl(model, ~ coef(.x)[["richness"]]),
    type = "obs",
    preds = map(model, ~ {
      richness_seq <- seq(0, max_richness, by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    })
  )

## unnest observed predictions
logistic_obs_preds <- logistic_obs_results %>%
  unnest(preds) %>%
  select(region, taxon_id, richness, max_richness, fit, type) %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

ggplot() +
  geom_line(data = logistic_obs_preds, linewidth = 1.5,
            aes(x = richness, y = fit, color = taxon_id, linetype = "Observed")) +
  facet_grid(rows = vars(region), cols = vars(taxon_id),
             labeller = labeller(.cols = mamm_labels)) +
  labs(x = "Species Richness", y = "Predicted Probability of Presence") +
  scale_colour_viridis_d(guide = "none") +
  scale_fill_viridis_d(guide = "none") +
  scale_linetype_manual(name = "Data Type", values = c("solid", "dashed"), breaks = c("Observed", "Simulated")) +
  scale_x_continuous(breaks = seq(1, max(logistic_obs_preds$richness), by = 1),
                     limits = c(1, max(logistic_obs_preds$richness))) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.position = "bottom")

## create ribbon of max and min fits data for simulated predictions
logistic_sim_pred_ribbon <- logistic_sim_results %>%
  unnest(preds) %>%
  group_by(region, taxon_id, richness, max_richness) %>%
  summarize(max = max(fit),
            min = min(fit),
            .groups = "drop") %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

mamm_labels <- c(
  PELE = "White-footed\nMouse",
  PEMA = "Deer Mouse",
  BLBR = "Short-tailed\nShrew",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)
  
## combine obs and sim
#logistic_all_preds <- bind_rows(logistic_obs_preds, logistic_sim_pred_ribbon)

ggplot() +
  geom_line(data = logistic_obs_preds, linewidth = 1.5,
            aes(x = richness, y = fit, color = taxon_id, linetype = "Observed")) +
  geom_ribbon(data = logistic_sim_pred_ribbon,
              aes(ymin = min, ymax = max, x = richness, fill = taxon_id, 
                  color = taxon_id, linetype = "Simulated"), 
              alpha = 0.25) +
  geom_line(data = logistic_sim_pred_ribbon,
            aes(x = richness, y = max, color = taxon_id, linetype = "Simulated"),
            alpha = 0.5) +
  facet_grid(rows = vars(region), cols = vars(taxon_id),
             labeller = labeller(.cols = mamm_labels)) +
  labs(x = "Species Richness", y = "Predicted Probability of Presence") +
  scale_colour_viridis_d(guide = "none") +
  scale_fill_viridis_d(guide = "none") +
  scale_linetype_manual(name = "Data Type", values = c("solid", "dashed"), breaks = c("Observed", "Simulated")) +
  scale_x_continuous(breaks = seq(1, max(logistic_obs_preds$richness), by = 1),
                     limits = c(1, max(logistic_obs_preds$richness))) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.position = "bottom")


## direct region comparison for pema and pele

# compare preds
logistic_pred_compare <- logistic_sim_results %>%
  unnest(preds) %>%
  rename(sim_fit = fit) %>%
  select(region, taxon_id, rep, richness, sim_fit) %>%
  group_by(region, taxon_id, richness) %>%
  mutate(mean_fit = mean(sim_fit),
         sd_fit = sd(sim_fit),
         scaled_fit = (sim_fit - mean_fit)/sd_fit) %>%
  left_join(logistic_obs_preds %>%
              rename(obs_fit = fit),
            by = c("region", "taxon_id", "richness")) %>%
  filter(taxon_id %in% c("PELE", "PEMA")) %>%
  mutate(fit_diff = obs_fit - sim_fit,
         scaled_fit_diff = (obs_fit - mean_fit)/sd_fit)

pred_compare_dists <- logistic_pred_compare %>%
  group_by(region, taxon_id, richness) %>%
  summarize(mean_scaled_diff = mean(scaled_fit_diff),
            upper_scaled_diff = quantile(scaled_fit_diff, probs = 0.95),
            lower_scaled_diff = quantile(scaled_fit_diff, probs = 0.05),
            mean_diff = mean(fit_diff),
            upper_diff = quantile(fit_diff, probs = 0.95),
            lower_diff = quantile(fit_diff, probs = 0.05),
            .groups = "drop")

ggplot(data = pred_compare_dists %>% filter(richness != 0),
       aes(x = richness, color = region)) +
  geom_point(aes(y = mean_scaled_diff)) +
  geom_errorbar(aes(ymin = lower_scaled_diff, ymax = upper_scaled_diff), width = 0.2) +
  facet_wrap(~ taxon_id, labeller = labeller(.cols = mamm_labels)) +
  labs(x = "Species Richness", y = "Mean Scaled Difference (Obs - Sim)") +
  theme_bw()

ggplot(data = pred_compare_dists %>% filter(richness != 0),
       aes(x = richness, color = region)) +
  geom_point(aes(y = mean_diff, shape = region), size = 3) +
  geom_line(aes(y = mean_diff, group = region), linetype = "dashed") +
  geom_errorbar(aes(ymin = lower_diff, ymax = upper_diff), width = 0.3, linewidth = 0.5) +
  facet_wrap(~ taxon_id, labeller = labeller(.cols = mamm_labels)) +
  scale_color_viridis_d(option = "turbo", name = "Region") +
  scale_shape_discrete(name = "Region") +
  labs(x = "Species Richness", y = "Difference in Predicted Presence \nProbability (Observed - Simulated)") +
  scale_x_continuous(breaks = seq(1, max(mamms_presence_filtered_df$richness), by = 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18))
  






ggplot(logistic_obs_preds, aes(x = richness, y = fit, color = taxon_id)) +
  geom_line() +
  facet_wrap(~region) +
  labs(x = "Species Richness", y = "Probability of Presence") +
  scale_x_continuous(breaks = seq(0, max(mamms_presence_filtered_df$richness), by = 1)) +
  theme_bw()


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
