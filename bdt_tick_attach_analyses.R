###
#  tick attachment by richness and mna
###

## loading packages
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(gt)
library(officer)

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

## mammals of interest (top 6 most common species across all sites)
mamms <- c("PELE", "PEMA", "MYGA", "TAST", "NAIN")

## summarize tick data into proportion of each species with ticks by plot and session
## then total number of mammals with tick data and number of those with ticks
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

## join in community variables like richness and minimum number alive
tick_attach_join <- tick_attach_summary %>%
  left_join(mammal_community_df %>% select(region, plot_session, year, mean_month,
                                           nlcd_class, taxon_id, mna, total_mna,
                                           richness, prop_mna), 
            by = c("plot_session", "taxon_id")) %>%
  mutate(site_id = str_sub(plot_session, 1, 4))


## now model the number of animals with ticks for each taxon and region
## using negative binomial GLMMs with richness as a fixed effect, NLCD as fixed due to limited classes, and random intercepts for year and mean_month
tick_attach_num_mods <- tick_attach_join %>%
  filter(taxon_id %in% mamms) %>%
  group_by(taxon_id, region) %>%
  nest() %>%
  mutate(
    num_mod = map(data, ~ glmmTMB(num_w_ticks ~ richness + nlcd_class + (1|site_id) + (1|year) + (1|mean_month),
                                  family = nbinom2, data = .x)),
    intercept = map_dbl(num_mod, ~ fixef(.x)$cond[1]),
    slope = map_dbl(num_mod, ~ fixef(.x)$cond[2]),
    p_value = map_dbl(num_mod, ~ summary(.x)$coefficients$cond[2,4]),
  ) %>%
  select(-data, -num_mod)



## modify model output table for plotting
tick_mod_table <- tick_attach_num_mods %>%
  ungroup() %>%
  # add nice names for plotting
  mutate(taxon_id = factor(taxon_id, levels = mamms),
         taxon_id = fct_recode(taxon_id, "White-footed Mouse" = "PELE",
                               "Deer Mouse" = "PEMA",
                               "Red-backed Vole" = "MYGA",
                               "Eastern Chipmunk" = "TAST",
                               "W. Jumping Mouse" = "NAIN"),
         signif = case_when(p_value < 0.001 ~ "***",
                            p_value < 0.01 ~ "**",
                            p_value < 0.05 ~ "*",
                            p_value < 0.1 ~ ".",
                            p_value > 0.05 ~ ""),
         p_value = round(p_value, 3),
         p_value = if_else(p_value < 0.001, "<0.001", as.character(p_value)),
         p_value = str_c(p_value, " ", signif))
  

## create table with model output for each species, stacked by region
mod_tbl <- tick_mod_table %>%
  arrange(region, taxon_id) %>%
  select(region, taxon_id, intercept, slope, p_value) %>%
  gt(groupname_col = "region") %>%
  fmt_number(columns = c(intercept, slope), decimals = 3) %>%
  cols_align(columns = c(intercept, slope, p_value), align = "center") %>%
  cols_label(
    taxon_id = md("**Species**"),
    intercept = md("**Intercept**"),
    slope = md("**Richness Slope**"),
    p_value = md("**Richness P-value**")
  ) %>%
  as_word()

doc <- read_docx() %>%
  body_add_xml(mod_tbl)

print(doc, target = "figures/bdt_table1_tick_mods.docx")

## bootstrap ticks randomly assembling on species by plot_session

## shuffle ticks column in tick attachment df 1000 times,
## essentially "re-attaching" ticks to random mammals
## takes a little while
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

## calculate max and min animals with ticks for each plot session and taxon based on shuffling
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
## for bootstrap comparisons (observed values to 95% CI of simulated values)
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

## summarize number and proportion of observed values that are higher, lower, or within 
## the 95% CI of simulated values for each taxon and region
tick_attach_sim_sum <- tick_attach_sim_join %>%
  filter(taxon_id %in% mamms) %>%
  group_by(taxon_id, region) %>%
  summarize(
    num_higher = sum(num_v_sim == "higher"),
    num_lower = sum(num_v_sim == "lower"),
    num_within = sum(num_v_sim == "within"),
    prop_higher = sum(prop_v_sim == "higher") / n(),
    prop_lower = sum(prop_v_sim == "lower") / n(),
    prop_within = sum(prop_v_sim == "within") / n(),
    .groups = "drop"
  )

## re order the taxon_id factor for plotting
tick_attach_sim_join <- tick_attach_sim_join %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

## create nice names for plotting
mamm_labels <- c(
  PELE = "White-footed\nMouse",
  PEMA = "Deer Mouse",
  BLBR = "Short-tailed\nShrew",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)

## visualize observed number of mammals with ticks by richness for each taxon and region, 
## colored by whether observed values are higher, lower, or within the 95% CI of simulated values
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


## join community variables and observed values, but to a nested data frame this time
## so we can calculate simulated means, std dev, and 
## standardized differences between observed data and simulated means for each plot session and taxon
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


## re order the taxon_id factor for plotting
tick_attach_sim_perc <- tick_attach_sim_perc %>%
  mutate(taxon_id = factor(taxon_id, levels = mamms))

## visualize standardized differences between observed and simulated proportions of mammals with ticks for each taxon and region
tick_plot <- ggplot(tick_attach_sim_perc %>% filter(taxon_id %in% mamms, prop_sd > 0),
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
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

ggsave(tick_plot, filename = "figures/bdt_fig4_ticks.pdf", width = 10, height = 8, units = "in", dpi = 300)





  
