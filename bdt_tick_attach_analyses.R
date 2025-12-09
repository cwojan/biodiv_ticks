###
#  tick attachment by richness and mna
###

## loading packages
library(tidyverse)
library(glmmTMB)

## read mammal community data
mammal_community_df <- read_rds("processed_data/mammal_community_df_2025-10-16.rds") %>%
  mutate(region = if_else(domain_id == "D05", "Upper Midwest", "Northeast")) %>%
  group_by(region) %>%
  mutate(max_richness = max(richness),
         max_mna = max(total_mna)) %>%
  ungroup()

## read tick attachment data
tick_attachment_df <- read_rds("processed_data/mammal_tick_captures_df_2025-10-21.rds") %>%
  mutate(region = if_else(site_id %in% c("TREE", "STEI", "UNDE"), "Upper Midwest", "Northeast"))

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
  left_join(mammal_community_df %>% select(region, plot_session, year, mean_month,
                                           nlcd_class, taxon_id, mna, total_mna,
                                           richness, prop_mna), 
            by = c("plot_session", "taxon_id")) %>%
  mutate(site_id = str_sub(plot_session, 1, 4))

check <- tick_attach_join %>%
  filter(taxon_id == "PELE")

## visualize tick attachment by richness and mna
ggplot(tick_attach_join %>% filter(taxon_id == "PELE")) +
  geom_jitter(width = 0, height = 0, 
              aes(x = total_mna, y = num_w_ticks,)) +
  geom_smooth(method = "lm",
              aes(x = total_mna, y = num_w_ticks), color = "black") +
  facet_wrap(~region) +
  scale_color_viridis_c() +
  #scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw()

ggplot(tick_attach_join %>% filter(taxon_id == "PELE")) +
  geom_jitter(width = 0.2, height = 0, 
              aes(x = richness, y = mna, color = prop_w_ticks)) +
  scale_color_viridis_c() +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw()

## test model
pele_mod <- glmmTMB(prop_w_ticks ~ richness + (1|nlcd_class) + (1|year),
                     family = "binomial", weights = num_total,
                     data = tick_attach_join %>%
                       filter(taxon_id == "PELE"))
summary(pele_mod)

pele_mod <- glm(prop_w_ticks ~ total_mna,
                family = "binomial", weights = num_total,
                data = tick_attach_join %>%
                  filter(taxon_id == "PELE"))
summary(pele_mod)

## generate prediction heat map
newdata <- expand.grid(
  richness = seq(0, 10, by = 1),
  num_total = seq(0, 50, by = 1),
  nlcd_class = NA,
  site_id = NA,
  year = NA
)
newdata$predicted_prop <- predict(pele_mod, newdata = newdata, type = "response")
newdata$predicted_num <- newdata$predicted_prop * newdata$num_total
ggplot(newdata, aes(x = richness, y = num_total, fill = predicted_prop)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Predicted\nProportion\nwith Ticks") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  labs(x = "Mammal Species Richness", y = "Number of Individuals Sampled") +
  theme_bw()

## bootstrap ticks randomly assembling on species by plot_session

## shuffle ticks column in tick attachment df 
tick_attach_sims <- replicate(n = 1000, simplify = FALSE,
                              expr = tick_attachment_df %>%
                                group_by(plot_session) %>%
                                mutate(ticks_shuffled = sample(ticks)) %>%
                                ungroup() %>%
                                group_by(region, plot_session, taxon_id) %>%
                                summarize(
                                  num_w_ticks = sum(ticks_shuffled == 1),
                                  num_total = n(),
                                  prop_w_ticks = num_w_ticks / num_total,
                                  .groups = "drop"
                                ) %>%
                                group_by(plot_session) %>%
                                mutate(
                                  all_sp_total = sum(num_total),
                                  all_sp_w_ticks = sum(num_w_ticks),
                                ) %>%
                                ungroup()) %>%
  bind_rows(.id = "sim_id")

## calculate max and mins for each plot session and taxon
tick_attach_sim_dists <- tick_attach_sims %>%
  group_by(region, plot_session, taxon_id) %>%
  summarize(
    upper_prop = quantile(prop_w_ticks, probs = 0.975),
    lower_prop = quantile(prop_w_ticks, probs = 0.025),
    upper_num = quantile(num_w_ticks, probs = 0.975),
    lower_num = quantile(num_w_ticks, probs = 0.025),
    min_prop = min(prop_w_ticks),
    max_prop = max(prop_w_ticks),
    min_num = min(num_w_ticks),
    max_num = max(num_w_ticks),
    .groups = "drop"
  )

## join community variables and observed values to sim dists
tick_attach_sim_join <- tick_attach_sim_dists %>%
  left_join(mammal_community_df %>% select(plot_session, year, mean_month,
                                           nlcd_class, taxon_id, mna, total_mna,
                                           richness, prop_mna), 
            by = c("plot_session", "taxon_id")) %>%
  left_join(tick_attach_summary %>%
              select(plot_session, taxon_id, prop_w_ticks, num_w_ticks, num_total),
            by = c("plot_session", "taxon_id")) %>%
  mutate(site_id = str_sub(plot_session, 1, 4),
         prop_v_sim = case_when(prop_w_ticks > upper_prop ~ "higher",
                                prop_w_ticks < lower_prop ~ "lower",
                                TRUE ~ "within"),
         num_v_sim = case_when(num_w_ticks > upper_num ~ "higher",
                               num_w_ticks < lower_num ~ "lower",
                               TRUE ~ "within")
         )

## alternate percentile approach
tick_attach_sim_perc <- tick_attach_sims %>%
  group_by(plot_session, taxon_id) %>%
  nest() %>%
  left_join(mammal_community_df %>% select(region, plot_session, year, mean_month,
                                           nlcd_class, taxon_id, mna, total_mna,
                                           richness, prop_mna), 
            by = c("plot_session", "taxon_id")) %>%
  left_join(tick_attach_summary %>%
              select(plot_session, taxon_id, prop_w_ticks, num_w_ticks, num_total),
            by = c("plot_session", "taxon_id")) %>%
  mutate(
    prop_mean = map_dbl(data, ~ mean(.x$prop_w_ticks)),
    num_mean = map_dbl(data, ~ mean(.x$num_w_ticks)),
    prop_sd = map_dbl(data, ~ sd(.x$prop_w_ticks)),
    num_sd = map_dbl(data, ~ sd(.x$num_w_ticks)),
    prop_dist = (prop_w_ticks - prop_mean) / prop_sd,
    num_dist = (num_w_ticks - num_mean) / num_sd,
    prop_percentile = map_dbl(data, ~ ecdf(.x$prop_w_ticks)(prop_w_ticks)),
    num_percentile = map_dbl(data, ~ ecdf(.x$num_w_ticks)(num_w_ticks))
    ) %>%
  select(-data)

check <- tick_attach_sim_perc %>%
  filter(prop_percentile == 1)

## visualize percentiles for all taxons across richness

mamm_labels <- c(
  PELE = "White-footed\nMouse",
  PEMA = "Deer Mouse",
  BLBR = "Short-tailed\nShrew",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)

tick_attach_sim_perc <- tick_attach_sim_perc %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

ggplot(tick_attach_sim_perc %>% filter(taxon_id %in% mamms, prop_sd > 0),
       aes(x = taxon_id, y = prop_dist, color = taxon_id, fill = taxon_id)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(rows = vars(region)) +
  scale_color_viridis_d(guide = "none") +
  scale_fill_viridis_d(guide = "none") +
  scale_x_discrete(labels = mamm_labels) +
  labs(x = "Mammal Species", y = "Standardized Difference in Proportion \nw/ Ticks to Simulated Mean") +
  theme_bw() +
  theme(legend.position = c(0.9,0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

## visualize results for all taxons across richness

tick_attach_sim_join <- tick_attach_sim_join %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

ggplot() +
  geom_jitter(data = tick_attach_sim_join %>% filter(taxon_id %in% mamms, prop_v_sim == "within"),
              width = 0.2, height = 0, 
              aes(x = richness, y = prop_w_ticks, 
                  color = "Within", alpha = "Within")) +
  geom_jitter(data = tick_attach_sim_join %>% filter(taxon_id %in% mamms, prop_v_sim == "higher"),
              width = 0.2, height = 0, size = 2,
              aes(x = richness, y = prop_w_ticks, 
                  color = "Higher", alpha = "Higher")) +
  geom_jitter(data = tick_attach_sim_join %>% filter(taxon_id %in% mamms, prop_v_sim == "lower"),
              width = 0.2, height = 0, size = 2,
              aes(x = richness, y = prop_w_ticks, 
                  color = "Lower", alpha = "Lower")) +
  scale_color_manual(name = "Observed vs \n95% Simulated CI",
                     values = c("Within" = "gray",
                                "Higher" = "red",
                                "Lower" = "blue")) +
  scale_alpha_manual(name = "Observed vs \n95% Simulated CI",
                     values = c("Within" = 0.5,
                                "Higher" = 1,
                                "Lower" = 1)) +
  labs(x = "Mammal Species Richness", y = "Proportion of Species Pop. w/ Ticks") +
  facet_wrap(~taxon_id, labeller = as_labeller(mamm_labels)) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw() +
  theme(legend.position = c(0.85,0.25),
        legend.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))


ggplot() +
  geom_jitter(data = tick_attach_sim_join %>% filter(taxon_id %in% mamms, prop_v_sim == "within"),
              width = 0.2, height = 0, 
              aes(x = richness, y = num_w_ticks, 
                  color = "Within", alpha = "Within")) +
  geom_jitter(data = tick_attach_sim_join %>% filter(taxon_id %in% mamms, prop_v_sim == "higher"),
              width = 0.2, height = 0, size = 2,
              aes(x = richness, y = num_w_ticks, 
                  color = "Higher", alpha = "Higher")) +
  geom_jitter(data = tick_attach_sim_join %>% filter(taxon_id %in% mamms, prop_v_sim == "lower"),
              width = 0.2, height = 0, size = 2,
              aes(x = richness, y = num_w_ticks, 
                  color = "Lower", alpha = "Lower")) +
  scale_color_manual(name = "Observed vs \n95% Simulated CI",
                     values = c("Within" = "gray",
                                "Higher" = "red",
                                "Lower" = "blue")) +
  scale_alpha_manual(name = "Observed vs \n95% Simulated CI",
                     values = c("Within" = 0.5,
                                "Higher" = 1,
                                "Lower" = 1)) +
  labs(x = "Mammal Species Richness", y = "Number of Mammals w/ Ticks") +
  facet_grid(cols = vars(taxon_id), rows = vars(region), labeller = labeller(.cols = mamm_labels)) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))







## look at raw sim data
tick_attach_sim_join_raw <- tick_attach_sims %>%
  left_join(mammal_community_df %>% select(plot_session, year, mean_month,
                                           nlcd_class, taxon_id, mna, total_mna,
                                           richness, prop_mna), 
            by = c("plot_session", "taxon_id")) %>%
  mutate(site_id = str_sub(plot_session, 1, 4))

ggplot(tick_attach_sim_join_raw %>% filter(taxon_id == "PELE")) +
  geom_boxplot(aes(x = richness, y = num_w_ticks, group = richness)) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw()

## compare to observed
ggplot(tick_attach_join %>% filter(taxon_id == "PELE")) +
  geom_boxplot(aes(x = richness, y = num_w_ticks, group = richness)) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw()

  
