##
# explore mna relationships
##

## loading packages
library(tidyverse)
library(neonUtilities)
library(janitor)

## read mammal data
mammal_session_df <- readRDS("processed_data/mammal_session_df.rds")
mna_by_session <- readRDS("processed_data/mna_by_session_corrected.rds")
mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## first, check richness and total mna relationships
mna_richness_df <- mna_by_session %>%
  select(-c(taxon_id, mna, presence, prop_mna)) %>%
  distinct()

mammal_session_df %>%
  filter(plot_session == "BART_001_0", taxon_id == "BLBR")

## viz total mna by rich
ggplot(mna_richness_df, aes(x = richness, y = total_mna, color = nlcd_class)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = mean(mna_richness_df$richness), linetype = "dashed", color = "black") +
  geom_hline(yintercept = mean(mna_richness_df$total_mna), linetype = "dashed", color = "black") +
  labs(x = "Species Richness", y = "Total MNA") +
  scale_x_continuous(breaks = seq(0, max(mna_by_session$richness), by = 1)) +
  theme_bw()

mna_rich_lm <- lm(total_mna ~ richness, data = mna_richness_df)
## add residuals to df
mna_richness_df <- mna_richness_df %>%
  mutate(residual = resid(mna_rich_lm))

## now by species
mna_rich_sp <- mna_by_session %>%
  filter(taxon_id %in% mamms) %>%
  left_join(mna_richness_df %>% select(plot_session, residual), by = "plot_session")

## viz sp mna by rich
ggplot(mna_rich_sp, aes(x = richness, y = mna, color = taxon_id)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Species Richness", y = "MNA per Species", color = "Species") +
  facet_wrap(~taxon_id) +
  scale_x_continuous(breaks = seq(0, max(mna_by_session$richness), by = 1)) +
  theme_bw()

## viz sp mna by total mna
ggplot(mna_rich_sp, aes(x = total_mna, y = mna, color = taxon_id)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Total MNA", y = "MNA per Species", color = "Species") +
  facet_wrap(~taxon_id) +
  theme_bw()

## viz sp prop by rich
ggplot(mna_rich_sp, aes(x = richness, y = prop_mna, color = taxon_id)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  labs(x = "Species Richness", y = "Prop", color = "Species") +
  facet_wrap(~taxon_id) +
  scale_x_continuous(breaks = seq(0, max(mna_by_session$richness), by = 1)) +
  theme_bw()

## viz sp prop by total mna
ggplot(mna_rich_sp, aes(x = total_mna, y = prop_mna, color = taxon_id)) +
  geom_point(alpha = 0.5) +
  labs(x = "Total MNA", y = "Prop", color = "Species") +
  facet_wrap(~taxon_id) +
  theme_bw()

## experimental...
ggplot(mna_rich_sp, aes(x = richness, y = log(prop_mna/total_mna), color = taxon_id)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  labs(x = "Species Richness", y = "MNA per Species", color = "Species") +
  facet_wrap(~taxon_id) +
  scale_x_continuous(breaks = seq(0, max(mna_by_session$richness), by = 1)) +
  theme_bw()

ggplot(mna_rich_sp, aes(x = residual, y = prop_mna, color = taxon_id)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  labs(x = "Deviation from MNA:Rich", color = "Species") +
  facet_wrap(~taxon_id) +
  theme_bw()
