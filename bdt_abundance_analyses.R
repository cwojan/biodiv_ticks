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
    mna_mod = map(data, ~ MASS::glm.nb(mna ~ richness, data = .)),
    rich_intercept = map_dbl(rich_mod, ~ coef(.)[1]),
    abund_intercept = map_dbl(abund_mod, ~ coef(.)[1]),
    mna_intercept = map_dbl(mna_mod, ~ coef(.)[1]),
    rich_coef = map_dbl(rich_mod, ~ coef(.)[2]),
    abund_coef = map_dbl(abund_mod, ~ coef(.)[2]),
    mna_coef = map_dbl(mna_mod, ~ coef(.)[2]),
    rich_intercept_upper = map_dbl(rich_mod, ~ confint(.x, "(Intercept)", level = 0.95)[2]),
    rich_intercept_lower = map_dbl(rich_mod, ~ confint(.x, "(Intercept)", level = 0.95)[1]),
    rich_coef_upper = map_dbl(rich_mod, ~ confint(.x, "richness", level = 0.95)[2]),
    rich_coef_lower = map_dbl(rich_mod, ~ confint(.x, "richness", level = 0.95)[1]),
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
    }),
    mna_preds = map(mna_mod, ~ {
      richness_seq <- seq(0, max(mammal_community_df$richness), by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    })
  ) %>%
  select(-data, -rich_mod, -abund_mod, -mna_mod)

rich_preds_obs <- species_prop_obs %>%
  select(taxon_id, rich_preds) %>%
  unnest(cols = c(rich_preds))

abund_preds_obs <- species_prop_obs %>%
  select(taxon_id, abund_preds) %>%
  unnest(cols = c(abund_preds))

mna_preds_obs <- species_prop_obs %>%
  select(taxon_id, mna_preds) %>%
  unnest(cols = c(mna_preds))

mamm_labels <- c(
  PELE = "White-footed\nMouse",
  PEMA = "Deer Mouse",
  BLBR = "Short-tailed\nShrew",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)

## quick look at total mna versus richness
mammal_community_summary <- mammal_community_df %>%
  group_by(plot_session) %>%
  summarize(
    richness = unique(richness),
    total_mna = unique(total_mna),
    .groups = "drop"
  )

ggplot(mammal_community_summary %>% filter(richness > 0), aes(x = richness, y = total_mna)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "red") +
  scale_x_continuous(breaks = seq(1, max(mammal_community_summary$richness), by = 1)) +
  labs(x = "Species Richness", y = "Total # of Mammals") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16))

## plot richness predictions by taxon
ggplot(rich_preds_obs, aes(x = richness, y = fit, color = taxon_id)) +
  geom_line(linewidth = 1) +
  geom_line(aes(x = richness, y = 1/richness), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#ccebc5", "BLBR" = "#a8ddb5",
                                "MYGA" = "#7bccc4", "TAST" = "#43a2ca", "NAIN" = "#0868ac"),
                     labels = mamm_labels,
                     name = "Species") +
  scale_x_continuous(breaks = seq(1, max(rich_preds_obs$richness), by = 1),
                     limits = c(1, max(rich_preds_obs$richness))) +
  labs(x = "Species Richness", y = "Predicted Proportion of Community") +
  theme_bw() + 
  theme(legend.position = c(0.9,0.7),
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

## plot abundance predictions by taxon
ggplot(abund_preds_obs, aes(x = total_mna, y = fit, color = taxon_id)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0, max(abund_preds_obs$total_mna), by = 5)) +
  labs(x = "Total MNA", y = "Predicted Proportion of Community") +
  theme_bw()

## plot mna predictions by taxon
ggplot(mna_preds_obs, aes(x = richness, y = fit, color = taxon_id)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(1, max(mna_preds_obs$richness), by = 1),
                     limits = c(1, max(mna_preds_obs$richness))) +
  labs(x = "Species Richness", y = "Predicted MNA per Species") +
  theme_bw()


## now we need to bootstrap random communities to compare the observed model results to...


## test in pieces
session_id <- "TREE_028_10"
species_observed <- mammal_community_df %>%
  filter(plot_session == session_id,
         presence == 1) %>%
  pull(taxon_id)

total_mna <- mammal_community_df %>%
  filter(plot_session == session_id) %>%
  pull(total_mna) %>%
  unique()

reps <- 5

sim_comms <- tibble(rep = rep(1:reps, each = total_mna),
                    species = as.vector(replicate(reps,
                                                  sample(species_observed, size = total_mna, replace = TRUE)))) %>%
  group_by(rep, species) %>%
  summarize(sim_mna = n(), .groups = "drop")

sample_community <- function(session_id, reps = 1000){
  ## get the processed species observed in the session
  species_observed <- mammal_community_df %>%
    filter(plot_session == session_id,
           presence == 1) %>%
    pull(taxon_id)
  
  ## get total number of processed individuals observed in the session
  total_mna <- mammal_community_df %>%
    filter(plot_session == session_id) %>%
    pull(total_mna) %>%
    unique()
  
  ## simulate community compositions for given session
  sim_comms <- tibble(rep = rep(1:reps, each = total_mna),
                      taxon_id = replicate(reps,
                                           sample(species_observed, 
                                                  size = total_mna, replace = TRUE)) %>%
                        as.vector()) %>%
    group_by(rep, taxon_id) %>%
    summarize(sim_mna = n(), .groups = "drop") %>%
    mutate(session_id = session_id, total_mna = total_mna)
  
  
  return(sim_comms)
}

sim_comms_list <- map(mammal_community_df %>% 
                        filter(total_mna > 0) %>%
                        pull(plot_session) %>%
                        unique(), 
                      sample_community, reps = 1000, .progress = TRUE)

sim_comms_df <- bind_rows(sim_comms_list)
rm(sim_comms_list)

sim_comms_join <- sim_comms_df %>%
  left_join(mammal_community_df %>%
              select(plot_session, richness) %>%
              distinct(),
            by = c("session_id" = "plot_session")) %>%
  mutate(prop_mna = sim_mna / total_mna)

summary(sim_comms_join %>% filter(taxon_id %in% mamms, rep == 1) %>% pull(taxon_id) %>% as.factor())

species_prop_sim <- sim_comms_join %>%
  filter(taxon_id %in% mamms) %>%
  group_by(taxon_id, rep) %>%
  nest() %>%
  arrange(rep, taxon_id) %>%
  mutate(
    rich_mod = map(data, ~ glm(prop_mna ~ richness, data = ., weights = total_mna, family = "binomial")),
    abund_mod = map(data, ~ glm(prop_mna ~ total_mna, data = ., weights = total_mna, family = "binomial")),
    mna_mod = map(data, ~ MASS::glm.nb(sim_mna ~ richness, data = .)),
    rich_intercept = map_dbl(rich_mod, ~ coef(.)[1]),
    abund_intercept = map_dbl(abund_mod, ~ coef(.)[1]),
    mna_intercept = map_dbl(mna_mod, ~ coef(.)[1]),
    rich_coef = map_dbl(rich_mod, ~ coef(.)[2]),
    abund_coef = map_dbl(abund_mod, ~ coef(.)[2]),
    mna_coef = map_dbl(mna_mod, ~ coef(.)[2])
  ) %>%
  select(-data, -rich_mod, -abund_mod, -mna_mod)

## summarize simulation results
species_prop_sim_summary <- species_prop_sim %>%
  group_by(taxon_id) %>%
  summarize(
    rich_intercept_mean = mean(rich_intercept),
    rich_intercept_upper = quantile(rich_intercept, 0.975),
    rich_intercept_lower = quantile(rich_intercept, 0.025),
    rich_coef_mean = mean(rich_coef),
    rich_coef_upper = quantile(rich_coef, 0.975),
    rich_coef_lower = quantile(rich_coef, 0.025),
    abund_intercept_mean = mean(abund_intercept),
    abund_intercept_upper = quantile(abund_intercept, 0.975),
    abund_intercept_lower = quantile(abund_intercept, 0.025),
    abund_coef_mean = mean(abund_coef),
    abund_coef_upper = quantile(abund_coef, 0.975),
    abund_coef_lower = quantile(abund_coef, 0.025),
    mna_intercept_mean = mean(mna_intercept),
    mna_intercept_upper = quantile(mna_intercept, 0.975),
    mna_intercept_lower = quantile(mna_intercept, 0.025),
    mna_coef_mean = mean(mna_coef),
    mna_coef_upper = quantile(mna_coef, 0.975),
    mna_coef_lower = quantile(mna_coef, 0.025),
    .groups = "drop"
  ) 

## reorder taxon_id factor levels
species_prop_obs <- species_prop_obs %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

species_prop_sim <- species_prop_sim %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

species_prop_sim_summary <- species_prop_sim_summary %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

## compare observed to simulated
ggplot(species_prop_sim, aes(x = taxon_id, color = taxon_id, fill = taxon_id)) +
  geom_jitter(alpha = 0.05,
             aes(y = rich_intercept, shape = "Simulated"),
             width = 0.1) +
  geom_errorbar(data = species_prop_sim_summary, 
                aes(ymin = rich_intercept_lower, 
                    ymax = rich_intercept_upper), 
                width = 0.2, linewidth = 1) +
  geom_errorbar(data = species_prop_obs, 
                aes(ymin = rich_intercept_lower, 
                    ymax = rich_intercept_upper), 
                width = 0.2, linewidth = 1.5,
                color = "black") +
  geom_errorbar(data = species_prop_obs, 
                aes(ymin = rich_intercept_lower, 
                    ymax = rich_intercept_upper), 
                width = 0.2, linewidth = 1) +
  geom_point(data = species_prop_obs,
             aes(y = rich_intercept, shape = "Observed"), color = "black", size = 4) +
  labs(x = "Taxon ID", y = "Species Proportion ~ Richness \nBinomial Intercept") +
  scale_shape_manual(name = "Type", values = c("Observed" = 21, "Simulated" = 16),
                     guide = guide_legend(override.aes = list(alpha = 1, fill = "grey"))) +
  scale_alpha_manual(values = c(1, 0.1)) +
  scale_size_manual(values = c(2, 0.5), guide = "none") +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#008080", "BLBR" = "#008080",
                                "MYGA" = "#008080", "TAST" = "#008080", "NAIN" = "#008080"),
                     guide = "none") +
  scale_fill_manual(values = c("PELE" = "#E66101", "PEMA" = "#008080", "BLBR" = "#008080",
                               "MYGA" = "#008080", "TAST" = "#008080", "NAIN" = "#008080"),
                    guide = "none") +
  scale_x_discrete(labels = mamm_labels) +
  theme_bw() +
  theme(legend.position = c(0.9,0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))


ggplot(species_prop_sim, aes(x = taxon_id, color = taxon_id, fill = taxon_id)) +
  geom_jitter(alpha = 0.05,
              aes(y = rich_coef, shape = "Simulated"),
              width = 0.1) +
  geom_errorbar(data = species_prop_sim_summary, 
                aes(ymin = rich_coef_lower, 
                    ymax = rich_coef_upper), 
                width = 0.2, linewidth = 1) +
  geom_errorbar(data = species_prop_obs, 
                aes(ymin = rich_coef_lower, 
                    ymax = rich_coef_upper), 
                width = 0.2, linewidth = 1.5,
                color = "black") +
  geom_errorbar(data = species_prop_obs, 
                aes(ymin = rich_coef_lower, 
                    ymax = rich_coef_upper), 
                width = 0.2, linewidth = 1) +
  geom_point(data = species_prop_obs,
             aes(y = rich_coef, shape = "Observed"), color = "black", size = 4) +
  labs(x = "Taxon ID", y = "Species Proportion ~ Richness \nBinomial Coefficient") +
  scale_shape_manual(name = "Type", values = c("Observed" = 21, "Simulated" = 16),
                     guide = guide_legend(override.aes = list(alpha = 1, fill = "grey"))) +
  scale_alpha_manual(values = c(1, 0.1)) +
  scale_size_manual(values = c(2, 0.5), guide = "none") +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#008080", "BLBR" = "#008080",
                                "MYGA" = "#008080", "TAST" = "#008080", "NAIN" = "#008080"),
                     guide = "none") +
  scale_fill_manual(values = c("PELE" = "#E66101", "PEMA" = "#008080", "BLBR" = "#008080",
                               "MYGA" = "#008080", "TAST" = "#008080", "NAIN" = "#008080"),
                    guide = "none") +
  scale_x_discrete(labels = mamm_labels) +
  theme_bw() +
  theme(legend.position = c(0.9,0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))


ggplot() +
  geom_point(data = species_prop_obs, aes(x = taxon_id, y = abund_intercept), color = "red", size = 3) +
  geom_errorbar(data = species_prop_sim_summary, 
                aes(x = taxon_id, 
                    ymin = abund_intercept_lower, 
                    ymax = abund_intercept_upper), 
                width = 0.2) +
  geom_point(data = species_prop_sim_summary, 
             aes(x = taxon_id, y = abund_intercept_mean), 
             color = "blue", size = 2) +
  labs(x = "Taxon ID", y = "Abundance Intercept") +
  theme_bw()

ggplot() +
  geom_point(data = species_prop_obs, aes(x = taxon_id, y = abund_coef), color = "red", size = 3) +
  geom_errorbar(data = species_prop_sim_summary, 
                aes(x = taxon_id, 
                    ymin = abund_coef_lower, 
                    ymax = abund_coef_upper), 
                width = 0.2) +
  geom_point(data = species_prop_sim_summary, 
             aes(x = taxon_id, y = abund_coef_mean), 
             color = "blue", size = 2) +
  labs(x = "Taxon ID", y = "Abundance Coefficient") +
  theme_bw()

ggplot() +
  geom_point(data = species_prop_obs, aes(x = taxon_id, y = mna_intercept), color = "red", size = 3) +
  geom_errorbar(data = species_prop_sim_summary, 
                aes(x = taxon_id, 
                    ymin = mna_intercept_lower, 
                    ymax = mna_intercept_upper), 
                width = 0.2) +
  geom_point(data = species_prop_sim_summary, 
             aes(x = taxon_id, y = mna_intercept_mean), 
             color = "blue", size = 2) +
  labs(x = "Taxon ID", y = "MNA Intercept") +
  theme_bw()

ggplot() +
  geom_point(data = species_prop_obs, aes(x = taxon_id, y = mna_coef), color = "red", size = 3) +
  geom_errorbar(data = species_prop_sim_summary, 
                aes(x = taxon_id, 
                    ymin = mna_coef_lower, 
                    ymax = mna_coef_upper), 
                width = 0.2) +
  geom_point(data = species_prop_sim_summary, 
             aes(x = taxon_id, y = mna_coef_mean), 
             color = "blue", size = 2) +
  labs(x = "Taxon ID", y = "MNA Coefficient") +
  theme_bw()


###

ggplot(data = sim_comms_join %>% filter(taxon_id %in% mamms, rep == 3),
       aes(x = richness, y = sim_mna)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~taxon_id) +
  scale_x_continuous(breaks = seq(1, max(sim_comms_join$richness), by = 1)) +
  labs(x = "Species Richness", y = "Simulated MNA") +
  theme_bw()

ggplot(data = sim_comms_join %>% filter(taxon_id %in% mamms, rep == 3)) +
  geom_histogram(aes(x = sim_mna), bins = 30) +
  facet_wrap(~taxon_id) +
  theme_bw()

ggplot(data = mammal_community_df %>% filter(taxon_id %in% mamms, mna > 0),
       aes(x = richness, y = mna)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~taxon_id) +
  scale_x_continuous(breaks = seq(1, max(sim_comms_join$richness), by = 1)) +
  labs(x = "Species Richness", y = "Simulated MNA") +
  theme_bw()

sim_mna_summary <- sim_comms_join %>%
  filter(taxon_id %in% mamms) %>%
  group_by(taxon_id, rep, richness) %>%
  summarize(
    mean_mna = mean(sim_mna),
    sd_mna = sd(sim_mna),
    lower_mna = quantile(sim_mna, 0.025),
    upper_mna = quantile(sim_mna, 0.975),
    .groups = "drop"
  )

ggplot(data = sim_mna_summary,
       aes(x = richness, y = mean_mna, color = taxon_id)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  facet_wrap(~taxon_id) +
  labs(x = "Species Richness", y = "MNA", color = "Species") +
  theme_bw()

obs_mna_summary <- mammal_community_df %>%
  filter(taxon_id %in% mamms, mna > 0) %>%
  group_by(taxon_id, richness) %>%
  summarize(
    mean_mna = mean(mna),
    upper_mna = quantile(mna, 0.975),
    lower_mna = quantile(mna, 0.025),
    .groups = "drop"
  )

## compare mna summaries
ggplot() +
  geom_jitter(data = sim_mna_summary, width = 0.1,
              aes(x = richness, y = mean_mna), alpha = 0.5) +
  geom_errorbar(data = obs_mna_summary, 
                aes(x = richness, 
                    ymin = lower_mna, 
                    ymax = upper_mna), 
                width = 0.2, alpha = 0.5) +
  geom_point(data = obs_mna_summary, 
             aes(x = richness, y = mean_mna), 
             color = "red", size = 2) +
  facet_wrap(~taxon_id) +
  scale_x_continuous(breaks = seq(1, max(sim_mna_summary$richness), by = 1)) +
  labs(x = "Species Richness", y = "Mean MNA") +
  theme_bw()






  
  

  
  
    