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

# load my processed data
mammal_comms <- read_csv("processed_data/mammal_community_effects.csv") %>%
  mutate(region = case_when(site_id == "BART" ~ "New Hampshire",
                            site_id == "HARV" ~ "Massachussetts",
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
  geom_point(aes(y = obs_pres_effect, color = obs_pres_eff_comp), size = 6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  facet_wrap(vars(region)) +
  scale_x_continuous(name = "Trapping Grid") +
  labs(y = "Deviation from Expected Proportion of \nTrapping Sessions Species is Present",
       color = "Observed Presence \nCompared to \nBootstrapped Communities:\n") +
  scale_color_manual(values = c("springgreen4", "black"), labels = c("Higher", "As Expected")) +
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

