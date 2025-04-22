## explore infected tick prevalence for teaching Biodiv and Health / Research Computing Exhibition

# load libraries
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

# load my processed data
mammal_comms <- read_csv("processed_data/mammal_community_effects.csv")


## load the tick pathogen data for each site with the map functional,
## with the filename format: "raw_data/SITE/SITE_tick_pathogen_status_2012-01-01_2024-12-31/tck_pathogen.csv"
## and the variable file format: "raw_data/SITE/SITE_small_mammal_box_trapping_2012-01-01_2024-12-31/variables_10092.csv"
tick_path_data <- map(.x = sites$site_code,
                      .f = function(site){
                        readTableNEON(dataFile = paste0("raw_data/", site, "/", site, "_tick_pathogen_status_2012-01-01_2024-12-31/tck_pathogen.csv"),
                                      varFile = paste0("raw_data/", site, "/", site, "_tick_pathogen_status_2012-01-01_2024-12-31/variables_10092.csv"))
                      })

## rewrite the above code with error catching if the directory doesn't exist
tick_path_data <- map(.x = sites$site_code,
                      .f = safely(function(site){
                        readTableNEON(dataFile = paste0("raw_data/", site, "/", site, "_tick_pathogen_status_2012-01-01_2024-12-31/tck_pathogen.csv"),
                                      varFile = paste0("raw_data/", site, "/", site, "_tick_pathogen_status_2012-01-01_2024-12-31/variables_10092.csv"))
                      },
                      otherwise = NULL,
                      quiet = TRUE
                      ))

## clean tick data
tick_path_data_clean <- map_df(tick_path_data, function(x){
  if(!is.null(x$result)){
    x$result %>%
      clean_names()
  } else {
    NULL
  }
})

## filter out only test_pathogen_name values that start with "Borrelia"
borrelia_tests <- tick_path_data_clean %>%
  filter(str_starts(test_pathogen_name, "Borrelia"))

## check individual counts
summary(borrelia_tests$individual_count)

## summarize proportion of positive tests by plot_id and test_pathogen_name
borrelia_tests_summary <- borrelia_tests %>%
  group_by(site_id, nlcd_class, test_pathogen_name) %>%
  filter(!is.na(test_result)) %>%
  summarize(n_positive = sum(test_result == "Positive", na.rm = TRUE),
            n_samples = n(),
            proportion_positive = n_positive / n_samples) %>%
  ungroup() %>%
  select(site_id, nlcd_class, test_pathogen_name, n_positive, n_samples, proportion_positive)

borrelia_sp_summary <- borrelia_tests_summary %>%
  filter(test_pathogen_name == "Borrelia sp.")

## filter PELE observations in mammal_comms
pele_comms <- mammal_comms %>%
  filter(taxon_id == "PELE",
         !is.na(obs_prop_eff_comp),
         !is.na(obs_pres_eff_comp))

pele_presence <- pele_comms %>%
  select(site_id, nlcd_class, obs_pres_effect, obs_pres_eff_comp) %>%
  unique() %>%
  rename(pele_pres = obs_pres_eff_comp)

## join the two dataframes
pele_borrelia <- pele_presence %>%
  left_join(borrelia_sp_summary, by = c("site_id", "nlcd_class"))

## summarize by nlcd_class
pele_borrelia_summary <- pele_borrelia %>%
  group_by(site_id, nlcd_class, n_samples, proportion_positive) %>%
  summarize(mean_pele_deviation = mean(obs_pres_effect, na.rm = TRUE))

ggplot(pele_borrelia_summary, aes(x = mean_pele_deviation, y = proportion_positive)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Average Deviation from Expected Proportion of Trapping Session PELE Present",
       y = "Proportion of Borrelia sp. positive tests",
       title = "Borrelia sp. positive tests vs PELE presence") +
  coord_cartesian(xlim = c(0, 0.7), ylim = c(0, 1)) +
  theme_bw()

quick_lm <- glm(proportion_positive ~ mean_pele_deviation,
                data = pele_borrelia_summary,
                family = binomial,
                weights = n_samples)
summary(quick_lm)

