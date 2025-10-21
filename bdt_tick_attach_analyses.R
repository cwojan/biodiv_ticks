###
#  tick attachment by richness and mna
###

## loading packages
library(tidyverse)
library(glmmTMB)

## read mammal community data
mammal_community_df <- read_rds("processed_data/mammal_community_df_2025-10-16.rds") %>%
  ungroup()

## read tick attachment data
tick_attachment_df <- read_rds("processed_data/mammal_tick_captures_df_2025-10-21.rds")

## mammals of interest
mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## summarize tick data into proportion of species with ticks
tick_attach_summary <- tick_attachment_df %>%
  group_by(plot_session, taxon_id) %>%
  summarize(
    num_w_ticks = sum(ticks == 1),
    num_total = n(),
    prop_w_ticks = num_w_ticks / num_total,
    .groups = "drop"
  ) %>%
  group_by(plot_session) %>%
  mutate(
    all_sp_total = sum(num_total),
    all_sp_w_ticks = sum(num_w_ticks),
  ) %>%
  ungroup()

## join in community variables 
tick_attach_join <- tick_attach_summary %>%
  left_join(mammal_community_df %>% select(plot_session, year, mean_month,
                                           nlcd_class, taxon_id, mna, total_mna,
                                           richness, prop_mna), 
            by = c("plot_session", "taxon_id")) %>%
  mutate(site_id = str_sub(plot_session, 1, 4))

check <- tick_attach_join %>%
  filter(taxon_id == "PELE")

## visualize tick attachment by richness and mna
ggplot(tick_attach_join) +
  geom_jitter(width = 0.2, height = 0, 
              aes(x = richness, y = num_total, color = num_w_ticks)) +
  facet_wrap(~taxon_id) +
  scale_color_viridis_c() +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw()

ggplot(tick_attach_join %>% filter(taxon_id == "PELE")) +
  geom_jitter(width = 0.2, height = 0, 
              aes(x = richness, y = num_total, color = prop_w_ticks)) +
  facet_wrap(~site_id) +
  scale_color_viridis_c() +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw()

## test model
pele_mod <- glmmTMB(prop_w_ticks ~ richness + num_total + (1|nlcd_class),
                     family = "binomial", weights = num_total,
                     data = tick_attach_join %>%
                       filter(taxon_id == "PELE"))
summary(pele_mod)

pele_mod <- glm(prop_w_ticks ~ richness + num_total,
                family = "binomial", weights = num_total,
                data = tick_attach_join %>%
                  filter(taxon_id == "PELE"))
summary(pele_mod)

## generate prediction heat map
newdata <- expand.grid(
  richness = seq(0, 10, by = 1),
  num_total = seq(0, 50, by = 1),
  nlcd_class = NA
)
newdata$predicted_prop <- predict(pele_mod, newdata = newdata, type = "response")
newdata$predicted_num <- newdata$predicted_prop * newdata$num_total
ggplot(newdata, aes(x = richness, y = num_total, fill = predicted_num)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Predicted\nProportion\nwith Ticks") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  labs(x = "Mammal Species Richness", y = "Number of Individuals Sampled") +
  theme_bw()

## bootstrap ticks randomly assembling on species by plot_session







  
