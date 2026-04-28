###
# file: bdt_abundance_analyses.R
# author: chris wojan
# description: analyze how the proportion of the community made up by each species changes with richness and compare to null model expectations
###

## loading packages
library(tidyverse)

## read mammal community data, make region variable, find max richness by region
mammal_community_df <- read_rds("processed_data/mammal_community_df_2026-04-24.rds") %>%
  mutate(region = if_else(domain_id == "D05", "Upper Midwest", "Northeast")) %>%
  group_by(region) %>%
  mutate(max_richness = max(richness),
         max_mna = max(total_mna)) %>%
  ungroup()

## mammals of interest, top 6 most common species  
mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## pretty labels for the mammals
mamm_labels <- c(
  PELE = "White-footed\nMouse",
  PEMA = "Deer Mouse",
  BLBR = "Short-tailed\nShrew",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)


## for each species, what proportion of the community does it make up when it is present?
## binomial regression of proportinal MNA on richness, with total_mna as weights
## using 1/richness as predictor to capture expected relationship of decreasing proportion with increasing richness
species_prop_obs <- mammal_community_df %>%
  filter(taxon_id %in% mamms, mna > 0) %>%
  group_by(taxon_id, region, max_richness, max_mna) %>%
  nest() %>%
  mutate(
    rich_mod = map(data, ~ glm(prop_mna ~ I(1/richness), data = ., weights = total_mna, family = "binomial")),
    rich_preds = map(rich_mod, ~ {
      richness_seq <- seq(0, max_richness, by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    })
  )

## unpack the predictions of the binomial models for plotting
rich_preds_obs <- species_prop_obs %>%
  select(taxon_id, rich_preds) %>%
  unnest(cols = c(rich_preds)) %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))


## now we need to bootstrap random communities to compare the observed model results to...

## create a function to randomly sample individuals from the species observed until the community size is reached
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

## create bootstrapped communities, transform to data frame
sim_comms_list <- map(mammal_community_df %>% 
                        filter(total_mna > 0) %>%
                        pull(plot_session) %>%
                        unique(), 
                      sample_community, reps = 1000, .progress = TRUE)
sim_comms_df <- bind_rows(sim_comms_list)
rm(sim_comms_list)

## join in richness and region info to simulated communities
sim_comms_join <- sim_comms_df %>%
  left_join(mammal_community_df %>%
              select(plot_session, richness, region, max_richness, max_mna) %>%
              distinct(),
            by = c("session_id" = "plot_session")) %>%
  mutate(prop_mna = sim_mna / total_mna)


## fit models to simulated data
## in three pieces (the last two take a while)
nested_sim_comms <- sim_comms_join %>%
  filter(taxon_id %in% mamms) %>%
  group_by(taxon_id, rep, region, max_richness) %>%
  nest() %>%
  arrange(rep, taxon_id)


species_prop_sim <- nested_sim_comms %>%
  mutate(
    rich_mod = map(data, ~ glm(prop_mna ~ I(1/richness), data = ., weights = total_mna, family = "binomial"),
                          .progress = TRUE)
  )

species_prop_sim <- species_prop_sim %>%
  mutate(
    rich_preds = map(rich_mod, ~ {
      richness_seq <- seq(0, max_richness, by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    }, .progress = TRUE)
  )



## unpack the predictions from the models for each simulated community, and create ranges for each richness value
rich_preds_sim <- species_prop_sim %>%
  select(-data) %>%
  unnest(cols = c(rich_preds)) %>%
  group_by(taxon_id, region, richness) %>%
  summarize(
    mean_fit = mean(fit),
    max = max(fit),
    min = min(fit),
    .groups = "drop"
  ) %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

## plot observed and simulated richness preds
abundance_plot <- ggplot() +
  geom_line(data = rich_preds_obs, linewidth = 1.5,
            aes(x = richness, y = fit, color = taxon_id, linetype = "Observed")) +
  geom_ribbon(data = rich_preds_sim,
              aes(ymin = min, ymax = max, x = richness, fill = taxon_id, 
                  color = taxon_id, linetype = "Simulated"), 
              alpha = 0.25) +
  geom_line(data = rich_preds_sim,
            aes(x = richness, y = max, color = taxon_id, linetype = "Simulated"),
            alpha = 0.5) +
  facet_grid(rows = vars(region), cols = vars(taxon_id),
             labeller = labeller(.cols = mamm_labels)) +
  labs(x = "Species Richness", y = "Predicted Proportion \nof Community") +
  scale_colour_viridis_d(guide = "none") +
  scale_fill_viridis_d(guide = "none") +
  scale_linetype_manual(name = "Data Type", values = c("solid", "dashed"), breaks = c("Observed", "Simulated")) +
  scale_x_continuous(breaks = seq(1, max(rich_preds_obs$richness), by = 2),
                     limits = c(1, max(rich_preds_obs$richness))) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.position = "bottom")

ggsave(abundance_plot, filename = "figures/bdt_fig3_abundance.pdf", width = 12, height = 6, units = "in", dpi = 300)


