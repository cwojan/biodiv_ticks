##
# analyze tick attachment by species and richness
##

## loading packages
library(tidyverse)
library(neonUtilities)
library(janitor)
library(lubridate)

## read mammal data
mammal_session_df <- readRDS("processed_data/mammal_session_df.rds")

mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## create any ticks attached column, remove any capture with ambiguous tick results
mammal_session_df <- mammal_session_df %>%
  filter(taxon_id %in% mamms) %>%
  mutate(ticks = if_any(ends_with("ticks_attached"), ~ . == "Y"),
         no_ticks = if_all(ends_with("ticks_attached"), ~ . == "N"),
         unk_ticks = if_any(ends_with("ticks_attached"), ~ . == "U")) %>%
  filter(!is.na(ticks), !is.na(no_ticks), !is.na(unk_ticks), unk_ticks == FALSE) %>%
  mutate(ticks = as.numeric(ticks))

## for recaptures, take the max ticks attached value
mammal_ticks <- mammal_session_df %>%
  select(site_id, plot_id, nlcd_class, tag_id, taxon_id, year, mean_month, mean_yday, plot_session, ticks) %>%
  group_by(plot_session, tag_id) %>%
  summarize(ticks = max(ticks),
            site_id = first(site_id),
            plot_id = first(plot_id),
            nlcd_class = first(nlcd_class),
            taxon_id = first(taxon_id),
            year = first(year),
            mean_month = first(mean_month),
            mean_yday = first(mean_yday),
            .groups = "drop")

tick_results <- mammal_ticks %>%
  group_by(taxon_id) %>%
  nest() %>%
  mutate(model = map(data, ~ glm(ticks ~ 1, data = ., family = "binomial")),
         predicted_value = map_dbl(model, ~ plogis(coef(.)[1])),
         upper_ci = map_dbl(model, ~ plogis(confint(.)[2])),
         lower_ci = map_dbl(model, ~ plogis(confint(.)[1]))) %>%
  select(taxon_id, predicted_value, upper_ci, lower_ci) %>%
  mutate(taxon_id = as.factor(taxon_id),
         taxon_id = fct_relevel(taxon_id, mamms))

tick_results_filt <- mammal_ticks %>%
  filter(!nlcd_class %in% c("pastureHay", "cultivatedCrops"), mean_yday >= 125) %>%
  group_by(taxon_id) %>%
  nest() %>%
  mutate(model = map(data, ~ glm(ticks ~ 1, data = ., family = "binomial")),
         predicted_value = map_dbl(model, ~ plogis(coef(.)[1])),
         upper_ci = map_dbl(model, ~ plogis(confint(.)[2])),
         lower_ci = map_dbl(model, ~ plogis(confint(.)[1]))) %>%
  select(taxon_id, predicted_value, upper_ci, lower_ci) %>%
  mutate(taxon_id = as.factor(taxon_id),
         taxon_id = fct_relevel(taxon_id, mamms))

## visualize the results

mamm_labels <- c(
  PELE = "White-footed\nMouse",
  PEMA = "Deer Mouse",
  BLBR = "Short-tailed\nShrew",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)

ggplot(tick_results, aes(x = taxon_id, y = predicted_value, color = taxon_id)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#008080", "BLBR" = "#008080",
                                "MYGA" = "#008080", "TAST" = "#008080", "NAIN" = "#008080"),
                     guide = "none") +
  scale_x_discrete(labels = mamm_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Species", y = "Proportion of Individuals\nwith Ticks Attached") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

ggplot(tick_results_filt, aes(x = taxon_id, y = predicted_value, color = taxon_id)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#008080", "BLBR" = "#008080",
                                "MYGA" = "#008080", "TAST" = "#008080", "NAIN" = "#008080"),
                     guide = "none") +
  scale_x_discrete(labels = mamm_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Species", y = "Proportion of Individuals\nwith Ticks Attached") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

blbr_ticks <- mammal_ticks %>%
  filter(taxon_id == "BLBR")


## tick attachment by richness

## read mna session data
mna_session_df <- readRDS("processed_data/mna_by_session_corrected.rds") %>%
  filter(taxon_id %in% mamms) %>%
  select(plot_session, taxon_id, mna, total_mna, prop_mna, richness)

mamm_tick_comm <- left_join(mammal_ticks, mna_session_df, by = c("plot_session", "taxon_id"))

tick_by_rich_preds <- mamm_tick_comm %>%
  group_by(taxon_id) %>%
  nest() %>%
  mutate(model = map(data, ~ glm(ticks ~ richness, data = ., family = "binomial")),
         preds = map(model, ~ {
           richness_seq <- seq(0, max(mamm_tick_comm$richness), by = 1)
           pred_df <- tibble(richness = richness_seq)
           pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
           pred_df$coef <- summary(.x)$coefficients["richness", "Estimate"]
           pred_df$pval <- summary(.x)$coefficients["richness", "Pr(>|z|)"]
           return(pred_df)
         })) %>%
  select(taxon_id, preds) %>%
  unnest(cols = c(preds)) %>%
  mutate(taxon_id = as.factor(taxon_id),
         taxon_id = fct_relevel(taxon_id, mamms),
         significance = if_else(pval < 0.05, "p < 0.05", "p >= 0.05"))

tick_by_rich_preds_filt <- mamm_tick_comm %>%
  filter(!nlcd_class %in% c("pastureHay", "cultivatedCrops"), mean_yday >= 125) %>%
  group_by(taxon_id) %>%
  nest() %>%
  mutate(model = map(data, ~ glm(ticks ~ richness, data = ., family = "binomial")),
         preds = map(model, ~ {
           richness_seq <- seq(0, max(mamm_tick_comm$richness), by = 1)
           pred_df <- tibble(richness = richness_seq)
           pred_df$fit <- predict(.x, newdata = pred_df, type = "response")
           pred_df$coef <- summary(.x)$coefficients["richness", "Estimate"]
           pred_df$pval <- summary(.x)$coefficients["richness", "Pr(>|z|)"]
           return(pred_df)
         })) %>%
  select(taxon_id, preds) %>%
  unnest(cols = c(preds)) %>%
  mutate(taxon_id = as.factor(taxon_id),
         taxon_id = fct_relevel(taxon_id, mamms),
         significance = if_else(pval < 0.05, "p < 0.05", "p >= 0.05"))

## visualize tick attachment by richness

ggplot(tick_by_rich_preds, aes(x = richness, y = fit, color = taxon_id)) +
  geom_line(aes(linetype = significance), linewidth = 1) +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#ccebc5", "BLBR" = "#a8ddb5",
                                "MYGA" = "#7bccc4", "TAST" = "#43a2ca", "NAIN" = "#0868ac"),
                     labels = mamm_labels,
                     name = "Species") +
  scale_linetype_manual(values = c("p < 0.05" = "solid", "p >= 0.05" = "dashed"),
                        name = "Significance") +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(1, max(mamm_tick_comm$richness))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Small Mammal Species Richness", y = "Predicted Probability of\nTick Attachment") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.background = element_rect(fill = "white", color = "black"))

ggplot(tick_by_rich_preds_filt, aes(x = richness, y = fit, color = taxon_id)) +
  geom_line(aes(linetype = significance), linewidth = 1) +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#ccebc5", "BLBR" = "#a8ddb5",
                                "MYGA" = "#7bccc4", "TAST" = "#43a2ca", "NAIN" = "#0868ac"),
                     labels = mamm_labels,
                     name = "Species") +
  scale_linetype_manual(values = c("p < 0.05" = "solid", "p >= 0.05" = "dashed"),
                        name = "Significance") +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(1, max(mamm_tick_comm$richness))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Small Mammal Species Richness", y = "Predicted Probability of\nTick Attachment") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.background = element_rect(fill = "white", color = "black"))

## tick "preference" analysis

# first let's check if the total mna column is accurate for mammals with tick info
mna_check <- mamm_tick_comm %>%
  group_by(plot_session) %>%
  summarize(total_mna_check = n(),
            total_mna_reported = first(total_mna),
            .groups = "drop") %>%
  filter(total_mna_check != total_mna_reported)

## many mismatches, so make a new total mna column based on mammal ticks data
mamm_tick_pref <- mamm_tick_comm %>%
  group_by(plot_session) %>%
  mutate(total_tick_mna = n()) %>%
  group_by(plot_session, taxon_id) %>%
  mutate(mna_tick = n(),
         prop_mna_tick = mna_tick / total_tick_mna,
         n_w_ticks = sum(ticks == 1),
         total_w_ticks = sum(ticks),
         prop_w_ticks = n_w_ticks / total_w_ticks,
         tick_pref = prop_w_ticks / prop_mna_tick) %>%
  filter(total_w_ticks > 0) %>%
  ungroup()

tick_pref_summary <- mamm_tick_pref %>%
  group_by(taxon_id) %>%
  summarize(mean_tick_pref = mean(tick_pref, na.rm = TRUE),
            upper = quantile(tick_pref, 0.975, na.rm = TRUE),
            lower = quantile(tick_pref, 0.025, na.rm = TRUE),
            n = n(),
            .groups = "drop") %>%
  mutate(taxon_id = as.factor(taxon_id),
         taxon_id = fct_relevel(taxon_id, mamms))

ggplot(tick_pref_summary, aes(x = taxon_id, y = mean_tick_pref, color = taxon_id)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#008080", "BLBR" = "#008080",
                                "MYGA" = "#008080", "TAST" = "#008080", "NAIN" = "#008080"),
                     guide = "none") +
  scale_x_discrete(labels = mamm_labels) +
  scale_y_continuous(limits = c(0, max(tick_pref_summary$upper)), breaks = seq(0, max(tick_pref_summary$upper), by = 1)) +
  labs(x = "Species", y = "Mean Tick Preference Index") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16))

tick_pref_sum_by_rich <- mamm_tick_pref %>%
  group_by(taxon_id, richness) %>%
  summarize(mean_tick_pref = mean(tick_pref, na.rm = TRUE),
            upper = quantile(tick_pref, 0.975, na.rm = TRUE),
            lower = quantile(tick_pref, 0.025, na.rm = TRUE),
            n = n(),
            .groups = "drop") %>%
  mutate(taxon_id = as.factor(taxon_id),
         taxon_id = fct_relevel(taxon_id, mamms))

ggplot(tick_pref_sum_by_rich, aes(x = richness, y = mean_tick_pref, color = taxon_id)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = taxon_id), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("PELE" = "#E66101", "PEMA" = "#ccebc5", "BLBR" = "#a8ddb5",
                                "MYGA" = "#7bccc4", "TAST" = "#43a2ca", "NAIN" = "#0868ac"),
                     labels = mamm_labels,
                     name = "Species") +
  scale_fill_manual(values = c("PELE" = "#E66101", "PEMA" = "#ccebc5", "BLBR" = "#a8ddb5",
                               "MYGA" = "#7bccc4", "TAST" = "#43a2ca", "NAIN" = "#0868ac"),
                    guide = "none") +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(1, max(mamm_tick_comm$richness))) +
  scale_y_continuous(limits = c(0, max(tick_pref_sum_by_rich$upper)), breaks = seq(0, max(tick_pref_sum_by_rich$upper), by = 1)) +
  labs(x = "Small Mammal Species Richness", y = "Mean Tick Preference Index") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.position = "bottom",
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.background = element_rect(fill = "white", color = "black"))

## ticks on pele analysis (bad method)

ticks_on_pele <- mamm_tick_comm %>%
  group_by(plot_session) %>%
  mutate(species_present = list(taxon_id)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(pele_present = ("PELE" %in% species_present)) %>%
  filter(pele_present == TRUE) %>%
  group_by(plot_session) %>%
  summarize(n_mamms = n(),
            n_mamms_w_ticks = sum(ticks),
            n_pele_w_ticks = sum(taxon_id == "PELE" & ticks == 1),
            prop_pele_w_ticks = n_pele_w_ticks / n_mamms_w_ticks,
            richness = first(richness),
            .groups = "drop")

check <- ticks_on_pele %>%
  filter(n_mamms_w_ticks > 0, richness == 1)

pele_richness_model <- glm(prop_pele_w_ticks ~ richness, data = ticks_on_pele %>% filter(n_mamms_w_ticks > 0), 
                           family = "binomial", weights = n_mamms_w_ticks)
summary(pele_richness_model)
str(pele_richness_model)
pele_richness_preds <- tibble(
  richness = seq(0, max(ticks_on_pele$richness), by = 1),
  fit = predict(pele_richness_model, newdata = tibble(richness = seq(0, max(ticks_on_pele$richness), by = 1)), type = "response"),
) 

pele_richness_preds <- predict(pele_richness_model, newdata = tibble(richness = seq(0, max(ticks_on_pele$richness), by = 1)), type = "response", se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(richness = seq(0, max(ticks_on_pele$richness), by = 1),
         lower_ci = fit - (1.96 * se.fit),
         upper_ci = fit + (1.96 * se.fit))
pele_richness_obs_points <- ticks_on_pele %>%
  filter(n_mamms_w_ticks > 0) %>%
  select(richness, prop_pele_w_ticks)
summary(pele_richness_obs_points)

reference_pele_preds <- tibble(
  richness = seq(0, max(ticks_on_pele$richness), by = 1),
  pred = 1 / richness
 )


ggplot(pele_richness_preds, aes(x = richness, y = fit)) +
  geom_point(data = pele_richness_obs_points, aes(x = richness, y = prop_pele_w_ticks), color = "#E66101", size = 1, alpha = 0.2) +
  geom_line(data = reference_pele_preds, aes(x = richness, y = pred), color = "black", linetype = "dashed") +
  geom_line(color = "#E66101", size = 1) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "#E66101") +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(1, max(ticks_on_pele$richness))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Small Mammal Species Richness", y = "Predicted Percent of Mammals with\nTicks that are White-footed Mice") +
  theme_bw() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 16))

###

ticks_on_pele_filt <- mamm_tick_comm %>%
  filter(!nlcd_class %in% c("pastureHay", "cultivatedCrops"), mean_yday >= 125) %>%
  group_by(plot_session) %>%
  mutate(species_present = list(taxon_id)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(pele_present = ("PELE" %in% species_present)) %>%
  filter(pele_present == TRUE) %>%
  group_by(plot_session) %>%
  summarize(n_mamms = n(),
            n_mamms_w_ticks = sum(ticks),
            n_pele_w_ticks = sum(taxon_id == "PELE" & ticks == 1),
            prop_pele_w_ticks = n_pele_w_ticks / n_mamms_w_ticks,
            richness = first(richness),
            .groups = "drop")

check <- ticks_on_pele_filt %>%
  filter(n_mamms_w_ticks > 0, richness == 1)

pele_richness_model_filt <- glm(prop_pele_w_ticks ~ richness, data = ticks_on_pele_filt %>% filter(n_mamms_w_ticks > 0), 
                           family = "binomial", weights = n_mamms_w_ticks)


pele_richness_preds_filt <- predict(pele_richness_model_filt, newdata = tibble(richness = seq(0, max(ticks_on_pele_filt$richness), by = 1)), 
                               type = "response", se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(richness = seq(0, max(ticks_on_pele_filt$richness), by = 1),
         lower_ci = fit - (1.96 * se.fit),
         upper_ci = fit + (1.96 * se.fit))
pele_richness_obs_points_filt <- ticks_on_pele_filt %>%
  filter(n_mamms_w_ticks > 0) %>%
  select(richness, prop_pele_w_ticks)
summary(pele_richness_obs_points_filt)


ggplot(pele_richness_preds_filt, aes(x = richness, y = fit)) +
  geom_point(data = pele_richness_obs_points_filt, aes(x = richness, y = prop_pele_w_ticks), color = "#E66101", size = 1, alpha = 0.2) +
  geom_line(color = "#E66101", size = 1) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "#E66101") +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(1, max(ticks_on_pele$richness))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(x = "Small Mammal Species Richness", y = "Predicted Percent of Mammals with\nTicks that are White-footed Mice") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

