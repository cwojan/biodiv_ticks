##
#' file: mammal_trap_path_data_explore.R
#' description: explore the neon data on small mammal trapping and pathogen status
#' particularly how mammal diversity metrics relate to tick presence and pathogen status
##

## load packages
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

## try loading data from one site
mammal_trap_data_tree <- readTableNEON(dataFile = "raw_data/TREE/TREE_small_mammal_box_trapping_2012-01-01_2024-12-31/mam_pertrapnight.csv",
                                       varFile = "raw_data/TREE/TREE_small_mammal_box_trapping_2012-01-01_2024-12-31/variables_10072.csv",)

## that works, so we can try to load the data for each site with the map functional,
## with the filename format: "raw_data/SITE/SITE_small_mammal_box_trapping_2012-01-01_2024-12-31/mam_pertrapnight.csv"
## and the variable file format: "raw_data/SITE/SITE_small_mammal_box_trapping_2012-01-01_2024-12-31/variables_10072.csv"
mammal_trap_data <- map(.x = sites$site_code,
                        .f = function(site){
                          readTableNEON(dataFile = paste0("raw_data/", site, "/", site, "_small_mammal_box_trapping_2012-01-01_2024-12-31/mam_pertrapnight.csv"),
                                        varFile = paste0("raw_data/", site, "/", site, "_small_mammal_box_trapping_2012-01-01_2024-12-31/variables_10072.csv"))
                        })
## change the list of data into one data frame, clean the names
mammal_trap_df <- bind_rows(mammal_trap_data) %>%
  clean_names() %>%
  mutate(trap_status_code = str_sub(string = trap_status, start = 1, end = 1) #first character of trap_status
        )
## filter out only traps with captures
mammal_captures <- filter(mammal_trap_df, trap_status_code %in% c("4","5"))

