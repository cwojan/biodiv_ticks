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
  mutate(trap_status_code = str_sub(string = trap_status, start = 1, end = 1))

## filter out only traps with captures
mammal_captures <- filter(mammal_trap_df, trap_status_code %in% c("4","5"))

## pivot tick presence data longer so that there is a tick_presence column and a tick_life_stage column,
## based one the columns "larval_ticks_attached", "nymphal_ticks_attached", "adult_ticks_attached",
## but cut out the "ticks_attached" part of the column names
mammal_captures_long <- mammal_captures %>%
  pivot_longer(cols = c(larval_ticks_attached, nymphal_ticks_attached, adult_ticks_attached),
               names_to = "tick_life_stage",
               values_to = "tick_presence") %>%
  mutate(tick_life_stage = str_remove_all(string = tick_life_stage, pattern = "_ticks_attached"))

## just checking the breakdown of tick presence by life stage and site
ggplot(mammal_captures_long, aes(x = tick_life_stage, fill = tick_presence)) +
  geom_bar(position = "fill") +
  facet_wrap(~site_id)

## let's filter out the rows where tick presence is either missing ("NA") or unknown ("U")
mammal_captures_tick <- filter(mammal_captures_long, tick_presence %in% c("Y", "N"))

## now checking the breakdown of tick presence with only the rows where tick presence is known
ggplot(mammal_captures_tick, aes(x = tick_life_stage, fill = tick_presence)) +
  geom_bar(position = "fill") +
  facet_wrap(~site_id)

## let's do a similar plot, but with a facet grid including mammal species and tick life stage
ggplot(mammal_captures_tick, aes(x = taxon_id, fill = tick_presence)) +
  geom_bar(position = "fill") +
  coord_flip() +
  facet_grid(cols = vars(site_id), rows = vars(tick_life_stage))

## and with abundacne of observations shown
ggplot(mammal_captures_tick, aes(x = taxon_id, fill = tick_presence)) +
  geom_bar() +
  coord_flip() +
  facet_grid(cols = vars(site_id), rows = vars(tick_life_stage))

## things to note: BART and HARV's "PELEPEMA" taxon_id could be lumped in with "PESP"
## however, there aren't many PESP with tick data, so maybe unnecessary?
## short-tailed shrews ("BLBR") don't seem to have any ticks recorded (likely because they aren't fully processed)
## Most common species with tick data include PEMA, PELE, PESP, TAST, ZAHU, NAIN, MYGA
## all species can be used for diversity metrics, focus on PELE and TAST for tick presence
## tick presence on MYGA, NAIN, and ZAHU could be interesting as evidence of dilution
## tick presence on PEMA could represent a potential reservoir in D05

## next step, let's characterize the overall small mammal communities by site
## first, let's filter out the rows where the taxon_id is missing,
## and classify all species that are not PEMA, PELE, PESP, TAST, ZAHU, NAIN, or MYGA as "other"
## in a new column called "taxon_group"
mammal_captures_taxon <- mammal_captures_tick %>%
  mutate(taxon_group = case_when(taxon_id %in% c("PEMA", "PELE", "PESP", "PELEPEMA", "TAST", "ZAHU", "NAIN", "MYGA") ~ taxon_id,
                                 TRUE ~ "Other"))

## now let's make a bar plot of the number of captures by site and taxon id
ggplot(mammal_captures_taxon, aes(x = site_id, fill = taxon_group)) +
  geom_bar() +
  coord_flip()

## next steps: proportion of community represented by PELE as x axis, 
## proportion of PELE with ticks as y axis, expect a positive relationship
