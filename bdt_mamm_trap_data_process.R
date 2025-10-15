###
# mammal trap data processing script
###

## loading packages
library(tidyverse)
library(neonUtilities)
library(janitor)

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
saveRDS(mammal_session_df, file = paste0("processed_data/mammal_session_df_", current_date, ".rds"))


## now filter down to only captures for another df



