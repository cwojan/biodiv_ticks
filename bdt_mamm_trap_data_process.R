###
# mammal trap data processing script
###

## loading packages
library(tidyverse)
library(neonUtilities)
library(janitor)

## read reference data tables
## load data product ids
data_products <- read_csv("logistics_data/neon_data_ids.csv")
## reformat data product names
data_products$data_name <- make_clean_names(data_products$data_category)
## load sites of interest
sites <- read_csv("logistics_data/neon_sites.csv") %>%
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
  mutate(
    trap_status_code = str_sub(string = trap_status, start = 1, end = 1),
    trap_num = str_split_fixed(string = trap_coordinate, pattern = "[:alpha:]", n = 2)[,2] %>%
      str_pad(width = 2, side = "left", pad = "0"),
    trap_coordinate = str_c(str_sub(string = trap_coordinate, start = 1, end = 1),
                            trap_num),
    trapping_date = as.Date(collect_date),
    year = year(trapping_date),
    month = month(trapping_date),
    yday = yday(trapping_date),
    tag_suffix = str_split_fixed(string = tag_id, pattern = "\\.", n = 4)[,4],
    # ear tag or pit tag = TRUE
    true_tag = case_when(str_starts(string = tag_suffix, pattern = "O") == TRUE ~ FALSE,
                         is.na(tag_id) ~ NA,
                         .default = TRUE),
    # identified to species = TRUE
    true_id = case_when(taxon_id == "PELEPEMA" ~ FALSE,
                        str_ends(string = scientific_name, pattern = "sp.") ~ FALSE,
                        is.na(taxon_id) ~ NA,
                        .default = TRUE)
  )

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
  arrange(domain_id, plot_session, trapping_date, trap_coordinate)

## grab current date
current_date <- Sys.Date()

## save rds of full session data
write_rds(mammal_session_df, file = str_c("processed_data/mammal_session_df_", current_date, ".rds", sep = ""))

## now set up for richness and mna calculations

## save only actual species taxon_ids
taxon_ids <- mammal_session_df %>%
  filter(true_id == TRUE) %>%
  distinct(taxon_id) %>%
  pull()

## create a df with columns for session and taxon, with every possible combo
taxon_by_session <- mammal_session_df %>%
  distinct(domain_id, plot_id, year, mean_month, mean_yday, plot_session, nlcd_class) %>%
  # repeat each row by the number of unique taxa
  slice(rep(1:n(), each = length(taxon_ids))) %>%
  # add in taxon ids
  mutate(taxon_id = rep(taxon_ids, times = n_distinct(mammal_session_df$plot_session)))

## calculate mna for O tags (same individual can have two tags in a session)
mna_o_tags <- mammal_session_df %>%
  filter(true_id == TRUE, true_tag == FALSE, !is.na(taxon_id)) %>%
  group_by(plot_session, plot_date, taxon_id) %>%
  summarise(mna_o = n_distinct(tag_id),
            .groups = "drop") %>%
  group_by(plot_session, taxon_id) %>%
  summarise(one_night_max = max(mna_o),
            .groups = "drop")

## join those maxes in
session_mna_o <- left_join(taxon_by_session, mna_o_tags,
                               by = c("plot_session", "taxon_id"))

## now calculate mna for true tag ids (one tag per individual)
mna_true_tags <-  mammal_session_df %>%
  filter(!is.na(tag_id), true_id == TRUE, true_tag == TRUE) %>% # true tagged captures
  distinct(tag_id, plot_session, taxon_id) %>% # exclude recaptures
  arrange(plot_session, taxon_id) %>%
  group_by(plot_session, taxon_id) %>%
  summarize(mna_over_session = n(), .groups = "drop")

## join those mnas in
session_mna_all <- left_join(session_mna_o, mna_true_tags,
                                 by = c("plot_session", "taxon_id")) %>%
  mutate(mna_over_session = replace_na(mna_over_session, 0),
         one_night_max = replace_na(one_night_max, 0),
         mna = pmax(mna_over_session, one_night_max, na.rm = TRUE))

## now create df with richness and total mna per session
session_community_df <- session_mna_all %>%
  group_by(plot_session) %>%
  mutate(presence = as.numeric(mna > 0), 
         total_mna = sum(mna),
         richness = sum(mna > 0),
         prop_mna =  mna / total_mna)

## save that community data
write_rds(session_community_df, file = str_c("processed_data/mammal_community_df_", current_date, ".rds", sep = ""))


## now filter down a ticks data frame

## only for species of interest
mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")

## only true tick data (no unknowns), and summarize to "max" tick presence on any one night in a session
tick_captures_df <- mammal_session_df %>%
  mutate(ticks = if_any(ends_with("ticks_attached"), ~ . == "Y"),
         no_ticks = if_all(ends_with("ticks_attached"), ~ . == "N"),
         unk_ticks = if_any(ends_with("ticks_attached"), ~ . == "U")) %>%
  filter(!is.na(ticks), !is.na(no_ticks), !is.na(unk_ticks), unk_ticks == FALSE,
         true_id == TRUE, true_tag == TRUE) %>%
  mutate(ticks = as.numeric(ticks)) %>%
  select(site_id, plot_id, nlcd_class, tag_id, taxon_id, year, mean_month, mean_yday, plot_session, ticks) %>%
  group_by(plot_session, tag_id) %>%
  summarize(site_id = first(site_id),
            plot_id = first(plot_id),
            nlcd_class = first(nlcd_class),
            year = first(year),
            mean_month = first(mean_month),
            mean_yday = first(mean_yday),
            taxon_id = first(taxon_id),
            ticks = max(ticks),
            .groups = "drop")

## save that ticks data
write_rds(tick_captures_df, file = str_c("processed_data/mammal_tick_captures_df_", current_date, ".rds", sep = ""))
