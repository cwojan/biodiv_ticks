###
# streamlined presence by richness analysis by bootstrapping
##

## loading packages
library(tidyverse)
library(neonUtilities)
library(janitor)
library(lubridate)

## read reference data tables
## load data product ids
data_products <- read_csv("neon_data_ids.csv")
## reformat data product names
data_products$data_name <- make_clean_names(data_products$data_category)
## load sites of interest
sites <- read_csv("neon_sites.csv") %>%
  clean_names()

## load the data for each site with the map functional,
## with the filename format: "raw_data/SITE/SITE_small_mammal_box_trapping_2012-01-01_2024-12-31/mam_pertrapnight.csv"
## and the variable file format: "raw_data/SITE/SITE_small_mammal_box_trapping_2012-01-01_2024-12-31/variables_10072.csv"
mammal_trap_data <- map(.x = sites$site_code,
                        .f = function(site){
                          readTableNEON(dataFile = paste0("raw_data/", site, "/", site, "_small_mammal_box_trapping_2012-01-01_2024-12-31/mam_pertrapnight.csv"),
                                        varFile = paste0("raw_data/", site, "/", site, "_small_mammal_box_trapping_2012-01-01_2024-12-31/variables_10072.csv"))
                        })

## first changing the list we loaded into a data frame, and adding some helpful columns
mammal_trap_df <- bind_rows(mammal_trap_data) %>%
  clean_names() %>%
  mutate(trap_status_code = str_sub(string = trap_status, start = 1, end = 1),
         trapping_date = as.Date(collect_date),
         year = year(trapping_date),
         month = month(trapping_date),
         yday = yday(trapping_date))

## then we will identify trapping sessions for each plot using the difference of days among trapping nights
trapping_sessions <- select(mammal_trap_df, plot_id, trapping_date, yday, month) %>%
  distinct() %>%
  arrange(plot_id, trapping_date) %>%
  group_by(plot_id) %>%
  mutate(plot_date = str_c(plot_id, "_", trapping_date),
         date_diff = c(0, diff(trapping_date)),
         session = cumsum(date_diff > 10),
         plot_session = str_c(plot_id, "_", session)) %>%
  group_by(plot_session) %>%
  mutate(mean_yday = mean(yday),
         mean_month = round(mean(month))) %>%
  select(-month, - yday)

## and then the session data can be joined into the main data frame
mammal_session_df <- left_join(mammal_trap_df, trapping_sessions, 
                               by = c("plot_id", "trapping_date")) %>%
  mutate(taxon_id = as.factor(taxon_id),
         plot_session = as.factor(plot_session)) %>%
  filter(!str_detect(scientific_name, "sp."), taxon_id != "PELEPEMA")

## now we can create a df with columns for session and taxon, with every possible combo
taxon_by_session <- mammal_session_df %>%
  distinct(year, mean_month, mean_yday, plot_session, nlcd_class) %>%
  # repeat each row by the number of unique taxa
  slice(rep(1:n(), each = length(levels(mammal_session_df$taxon_id)))) %>%
  # add in taxon ids
  mutate(taxon_id = rep(levels(mammal_session_df$taxon_id), times = n_distinct(mammal_session_df$plot_session)))

## now summarize the minimum number alive (number of unique captures) for each taxon in each session
mna_summary <-  mammal_session_df %>%
  filter(!is.na(tag_id)) %>% # tagged captures
  distinct(tag_id, plot_session, taxon_id) %>% # exclude recaptures
  arrange(plot_session, taxon_id) %>%
  group_by(plot_session, taxon_id) %>%
  summarize(mna = n(), .groups = "drop")
  

## this can be joined into the taxon by session df
mna_by_session <- left_join(taxon_by_session, mna_summary, 
                            by = c("plot_session", "taxon_id")) %>%
  mutate(mna = replace_na(mna, 0)) %>%
  group_by(plot_session) %>%
  mutate(presence = as.numeric(mna > 0), 
         total_mna = sum(mna),
         richness = sum(mna > 0),
         prop_mna = if_else(total_mna > 0, mna / total_mna, 0)) %>%
  separate_wider_delim(cols = plot_session, names = c("site_id", "plot_num", "session"), 
                       delim = "_", cols_remove = FALSE) %>%
  mutate(session = as.numeric(session),
         plot_id = as.factor(str_c(site_id, "_", plot_num))) %>%
  ungroup()


## which species in which sites/nlcd_class
site_nlcd_summary <- mna_by_session %>%
  ungroup() %>%
  filter(mna > 0) %>%
  select(site_id, nlcd_class, taxon_id, mna) %>%
  group_by(site_id, nlcd_class, taxon_id) %>%
  summarize(total_mna = sum(mna)) %>%
  ungroup()

## species pools for each site and nlcd
site_nlcd_species_pools <- site_nlcd_summary %>%
  group_by(site_id, nlcd_class) %>%
  summarize(species_pool = list(taxon_id),
            n_species = n(),
            .groups = "drop")

## average pool size
ave_pool <- mean(site_nlcd_species_pool_sizes$n_species)

## create a function to randomly sample n species from the pool of species
## observed at a site and land cover class combination
sample_species <- function(plot, reps = 1000){
  ## grab site id and nlcd class for the plot
  site <- str_sub(string = plot, start = 1, end = 4) ## first four characters of plot_id
  nlcd <- mna_by_session %>%
    filter(plot_id == plot) %>%
    select(nlcd_class) %>%
    distinct() %>%
    pull(nlcd_class) 
  
  ## get the species pool for the site and nlcd class
  species_pool <- site_nlcd_summary %>%
    filter(site_id == site, nlcd_class == nlcd) %>%
    pull(taxon_id)
  
  ## get the observed richnesses for the plot's sessions
  session_richness <- mna_by_session %>%
    filter(plot_id == plot, richness > 0) %>%
    select(session, richness) %>%
    distinct() %>%
    arrange(session, richness)
  
  ## create a data frame of the species samples
  sample_comms <- expand_grid(session_richness, rep = 1:reps) %>%
    mutate(species_sample = map(richness, ~sample(species_pool, size = .x, replace = FALSE))) 
  
  ## create a data frame of the species samples with presence information
  out <- expand_grid(session_richness, species = species_pool) %>%
    uncount(weights = reps) %>%
    mutate(rep = rep(1:reps, times = nrow(session_richness) * length(species_pool))) %>%
    arrange(rep, session, species) %>%
    left_join(sample_comms, by = c("rep", "session", "richness")) %>%
    group_by(rep, session) %>%
    mutate(presence = as.numeric(species %in% unlist(species_sample)),
           plot_id = plot)
  
  return(out)
}

test <- sample_species(plot = "BART_001", rep = 1)

distinct(mna_by_session, plot_id)

## simulate communities for all plots
simulated_comms <- map(levels(mna_by_session$plot_id), sample_species, reps = 1000, .progress = TRUE)

## and bind the list of data frames into one data frame
simulated_comms_df <- bind_rows(simulated_comms)

rm(simulated_comms)

## next step run logistic regression for each species, each rep, of presence ~ richness
## generate distribution of effect sizes

mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## nest and mutate method suggested by chatgpt

logistic_sim_results <- simulated_comms_df %>%
  filter(species %in% mamms) %>%
  group_by(species, rep) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(presence ~ richness, data = ., family = binomial)),
    intercept = map_dbl(model, ~ coef(.x)[["(Intercept)"]]),
    coef  = map_dbl(model, ~ coef(.x)[["richness"]]),
    type = "sim"
  ) %>%
  select(species, rep, coef, intercept, type)

## now run the same logistic regression on the observed data
# first correct mna_by_session to only include rows where the species is in the relevant species pool
mna_by_session_corrected <- mna_by_session %>%
  left_join(site_nlcd_species_pools, by = c("site_id", "nlcd_class")) %>%
  rowwise() %>%
  filter(taxon_id %in% species_pool) %>%
  ungroup() %>%
  select(-species_pool, -n_species)


logistic_obs_results <- mna_by_session_corrected %>%
  filter(taxon_id %in% mamms) %>%
  group_by(taxon_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(presence ~ richness, data = ., family = binomial)),
    intercept = map_dbl(model, ~ coef(.x)[["(Intercept)"]]),
    coef  = map_dbl(model, ~ coef(.x)[["richness"]]),
    type = "obs",
    rep = "NA"
  ) %>%
  select(species = taxon_id, rep, coef, intercept, type)

## bind the simulated and observed results
logistic_sim_results <- logistic_sim_results %>%
  group_by(species) %>%
  mutate(rep = as.character(rep),
         prob_at_1 = plogis(intercept + coef * 1),
         species = as.factor(species),
         species = fct_relevel(species, mamms))
logistic_obs_results <- logistic_obs_results %>%
  group_by(species) %>%
  mutate(prob_at_1 = plogis(intercept + coef * 1),
         species = as.factor(species),
         species = fct_relevel(species, mamms))
logistic_results <- bind_rows(logistic_sim_results, logistic_obs_results) %>%
  mutate(species = as.factor(species),
         species = fct_relevel(species, mamms))

logistic_intervals <- logistic_sim_results %>%
  group_by(species) %>%
  summarize(lower_coef = quantile(coef, probs = 0.025),
            upper_coef = quantile(coef, probs = 0.975),
            lower_prob = quantile(prob_at_1, probs = 0.025),
            upper_prob = quantile(prob_at_1, probs = 0.975),
            .groups = "drop") %>%
  mutate(species = as.factor(species),
         species = fct_relevel(species, mamms))

## visualize the results

mamm_labels <- c(
  PELE = "White-footed\nMouse",
  PEMA = "Deer Mouse",
  BLBR = "Short-tailed\nShrew",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)

ggplot(logistic_sim_results, aes(x = species, color = species, fill = species)) +
  geom_jitter(aes(y = coef, shape = "Simulated"), alpha = 0.05,
              width = 0.1) +
  geom_errorbar(data = logistic_intervals, 
                aes(ymin = lower_coef, ymax = upper_coef),
                width = 0.2,
                linewidth = 1) +
  geom_point(data = logistic_obs_results,
             aes(y = coef, shape = "Observed"), color = "black", size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
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
  scale_y_continuous(limits = c(-0.1, max(logistic_results$coef))) +
  labs(y = "Logistic Regression Coefficient (Log Odds Ratio)", x = "Species") +
  theme_bw() +
  theme(legend.position = c(0.9,0.85),
        legend.background = element_rect(fill = "white", color = "black"))

ggplot(logistic_sim_results, aes(x = species, color = species, fill = species)) +
  geom_jitter(aes(y = prob_at_1 * 100, shape = "Simulated"), alpha = 0.05,
              width = 0.1) +
  geom_errorbar(data = logistic_intervals, 
                aes(ymin = lower_prob * 100, ymax = upper_prob * 100),
                width = 0.2,
                linewidth = 1) +
  geom_point(data = logistic_obs_results,
             aes(y = prob_at_1 * 100, shape = "Observed"), color = "black", size = 3) +
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
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     limits = c(0, 100)) +
  labs(y = "Presence Probability at Richness = 1", x = "Species") +
  theme_bw() +
  theme(legend.position = c(0.9,0.85),
        legend.background = element_rect(fill = "white", color = "black"))


logistic_obs_preds <- mna_by_session_corrected %>%
  filter(taxon_id %in% mamms) %>%
  group_by(taxon_id) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(presence ~ richness, data = ., family = binomial)),
    ## save fitted values with corresponding richness values
    preds = map(model, ~ {
      richness_seq <- seq(0, max(mna_by_session$richness), by = 1)
      pred_df <- tibble(richness = richness_seq)
      pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
      return(pred_df)
    })
  ) %>%
  select(species = taxon_id, preds) %>%
  unnest(cols = c(preds)) %>%
  mutate(species = as.factor(species),
         species = fct_relevel(species, mamms))



reference_preds <- tibble(richness = seq(0, max(mna_by_session$richness), by = 1)) %>%
  mutate(fit = richness / ave_pool)

ggplot() +
  geom_line(data = reference_preds,
            aes(x = richness, y = fit * 100), color = "black", linetype = "dashed", size = 1) +
  geom_line(data = logistic_obs_preds,
            aes(x = richness, y = fit * 100, color = species), size = 1) +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#ccebc5", "BLBR" = "#a8ddb5",
                                "MYGA" = "#7bccc4", "TAST" = "#43a2ca", "NAIN" = "#0868ac"),
                     labels = mamm_labels,
                     name = "Species") +
  scale_x_continuous(breaks = seq(0, max(mna_by_session$richness), by = 1)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     limits = c(0, 100)) +
  labs(y = "Predicted Presence Probability", x = "Species Richness") +
  theme_bw() +
  theme(legend.position = c(0.9,0.25),
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.background = element_rect(fill = "white", color = "black"))

## check observed presence of PELE across richness, visually
pele_obs <- mna_by_session_corrected %>%
  filter(taxon_id == "PELE") %>%
  group_by(richness) %>%
  summarize(presence_rate = mean(presence),
            n_sessions = n(),
            n_presences = sum(presence),
            .groups = "drop")






