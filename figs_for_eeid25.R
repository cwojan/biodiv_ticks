## create figures for EEID 2025 poster


###
# Setup
###

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


###
# ticks
###

## read tick pathogen data with error catching if the directory doesn't exist
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

###
# mammals
###

# load my processed data
mammal_comms <- read_csv("processed_data/mammal_community_effects.csv") %>%
  mutate(region = case_when(site_id == "BART" ~ "New Hampshire",
                            site_id == "HARV" ~ "Massachusetts",
                            site_id == "STEI" ~ "Northwest WI",
                            site_id == "TREE" ~ "Northeast WI",
                            site_id == "UNDE" ~ "Upper MI"))

## check most commonly present taxon_ids
mammal_comms %>%
  filter(!is.na(taxon_id)) %>%
  group_by(taxon_id) %>%
  summarise(sum = sum(presence)) %>%
  arrange(desc(sum)) %>%
  slice_head(n = 10)

## create pele only data
pele_pres <- mammal_comms %>%
  select(region, site_id, plot_id, nlcd_class, taxon_id, obs_pres_effect, sim_upper_pres, sim_lower_pres, obs_pres_eff_comp) %>%
  unique() %>%
  filter(taxon_id == "PELE") %>%
  arrange(obs_pres_effect) %>%
  group_by(site_id) %>%
  mutate(pres_rank = row_number()) %>%
  ungroup()

###
# mammal plot
###

## plot pele presence
ggplot(data = pele_pres, aes(x = pres_rank, color = obs_pres_eff_comp)) +
  geom_point(aes(y = obs_pres_effect)) +
  geom_errorbar(aes(ymin = sim_lower_pres, ymax = sim_upper_pres),
                alpha = 0.5, width = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(vars(region)) +
  scale_x_continuous(name = "Plot") +
  labs(y = "Deviation from Expected Proportion of \nTrapping Sessions Species is Present",
       color = "Observed Presence Compared \nto Bootstrapped Communities:") +
  scale_color_manual(values = c("springgreen4", "black"), labels = c("Higher", "As Expected")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8,0.2))

## create tiff
pele_pres_plot <- ggplot(data = pele_pres, aes(x = pres_rank)) +
  geom_errorbar(aes(ymin = sim_lower_pres, ymax = sim_upper_pres),
                alpha = 0.5, width = 0.5, linewidth = 1.5) +
  geom_point(aes(y = obs_pres_effect, color = obs_pres_eff_comp), size = 8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  facet_wrap(vars(region)) +
  scale_x_continuous(name = "Trapping Grid") +
  labs(y = expression(atop("Deviation from Expected Proportion of", "Trapping Sessions " * italic("P. leucopus") * " is Present")),
       color = "Observed Presence \nCompared to \nSimulated  Communities:\n") +
  scale_color_manual(values = c("#b73779", "black"), labels = c("Higher", "As Expected")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 32),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 24),
        legend.key.height = unit(0.5, "in"),
        legend.position = "inside",
        legend.position.inside = c(0.84,0.28))

timestamp <- format(Sys.time(), format = "%Y%m%d_%H%M%S")
ggsave(str_c("figures/eeid25_pele_pres", timestamp, ".tiff"), plot = pele_pres_plot,
       units = "in", width = 16, height = 14)

###
# tick plots
###

## summarize mean pele presence, richness by site and nlcd class
site_nlcd_summary <- mammal_comms %>%
  filter(taxon_id == "PELE") %>%
  group_by(site_id, region, nlcd_class) %>%
  summarize(mean_pele_pres = mean(obs_pres_effect, na.rm = TRUE),
            mean_richness = mean(precise_richness, na.rm = TRUE)) %>%
  ungroup()

## join with borrelia data
borrelia_mamm_data <- borrelia_sp_summary %>%
  left_join(site_nlcd_summary, by = c("site_id", "nlcd_class")) %>%
  mutate(habitat = case_when(nlcd_class == "deciduousForest" ~ "Deciduous Forest",
                             nlcd_class == "evergreenForest" ~ "Evergreen Forest",
                             nlcd_class == "mixedForest" ~ "Mixed Forest",
                             nlcd_class == "woodyWetlands" ~ "Woody Wetlands"))


## plot borrelia prevalence by mean pele presence
ggplot(data = borrelia_mamm_data, aes(x = mean_pele_pres, y = proportion_positive)) +
  geom_point(aes(color = habitat, size = n_samples)) +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.5)) +
  labs(x = "Mean Proportion of Trapping Sessions P. leucopus is Present",
       y = "Proportion of Borrelia sp. Positive Tests") +
  theme_bw()

## plot borrelia prevalence by mean richness
ggplot(data = borrelia_mamm_data, aes(x = mean_richness, y = proportion_positive)) +
  geom_point(aes(color = habitat, size = n_samples)) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "Mean Species Richness",
       y = "Proportion of Borrelia sp. Positive Tests") +
  theme_bw()

## combo plot
ggplot(data = borrelia_mamm_data, aes(x = mean_richness, y = mean_pele_pres)) +
  geom_point(aes(fill = proportion_positive, size = n_samples), shape = 21) +
  geom_label(aes(label = str_c(site_id, nlcd_class, sep = "-")), size = 5, nudge_x = 0.05, nudge_y = 0.05) +
  labs(y = expression(atop("Mean Proportion of Trapping Sessions", italic("P. leucopus") * " is Present")),
       x = "Mean Species Richness") +
  scale_fill_viridis_c(option = "magma", name = expression(atop("Proportion of " * italic("Borrelia sp. "), "Positive Tests"))) +
  scale_size_continuous(name = "Number of Ticks Tested", range = c(2,7), breaks = c(1, 100, 500), 
                        guide = guide_legend(nrow = 2)) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.spacing = unit(1, "in"),
        legend.key.width = unit(0.5, "in"))

## create tiff
borr_plot <- ggplot(data = borrelia_mamm_data, aes(x = mean_richness, y = mean_pele_pres)) +
  geom_point(aes(fill = proportion_positive, size = n_samples), shape = 21) +
  labs(y = expression(atop("Mean Deviation from Expected Proportion of", "Trapping Sessions " * italic("P. leucopus") * " is Present")),
       x = "Mean Species Richness") +
  scale_fill_viridis_c(option = "magma", name = expression(atop("Proportion of " * italic("Borrelia") * " sp.", "Positive Tests"))) +
  scale_size_continuous(name = "Number of Ticks Tested", range = c(5,25),
                        breaks = c(1, 100, 500)) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.spacing = unit(1, "in"),
        legend.key.width = unit(1, "in"),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 32),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 24),
        plot.margin = unit(c(1, 1, 1.5, 1.5), "cm"))

ggsave(str_c("figures/eeid25_borr_mamm", timestamp, ".tiff"), plot = borr_plot,
       units = "in", width = 16, height = 14)
