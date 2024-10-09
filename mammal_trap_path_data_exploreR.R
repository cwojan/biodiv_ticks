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
## scratch that, instead let's first look at how the relative proportions of taxons vary among
## trapping sessions with different total mammal community sizes

## one key part is identifying trapping sessions among the capture data
## first let's look at the capture dates across sites

## let's make a year column to easily visualize
mammal_captures_taxon <- mammal_captures_taxon %>%
  mutate(collect_date = as.Date(collect_date),
         year = year(collect_date),
         yday = yday(collect_date))

## visualizing trapping sessions
ggplot(mammal_captures_taxon, aes(x = yday, y = site_id)) +
  geom_point() +
  facet_wrap(vars(year))

## trapping session are stratified well, what we need to do is categorize sets of capture dates as trapping sessions
## by plot_id, such that collect_dates within a given number of days of each other are in the same trapping session
## let's start by calculating the number of days between each collect_date and the first collect_date in the plot
mammal_captures_taxon <- mammal_captures_taxon %>%
  group_by(plot_id) %>%
  mutate(trapping_day = as.numeric(collect_date - min(collect_date)))

ggplot(mammal_captures_taxon, aes(x = trapping_day, y = plot_id)) +
  geom_point() +
  coord_cartesian(xlim = c(1, 150))

## now let's identify "clusters" of values in the trapping session column
## so that we can group nearby numbers together
## first let's just separate out the plot_id and trapping_session columns into a new data frame
trapping_sessions <- select(mammal_captures_taxon, plot_id, trapping_day) %>%
  distinct() %>%
  arrange(plot_id, trapping_day) %>%
  mutate(plot_day = str_c(plot_id, "_", trapping_day),
         day_diff = c(0, diff(trapping_day)),
         session = cumsum(day_diff > 10),
         plot_session = str_c(plot_id, "_", session))

ggplot(trapping_sessions, aes(x = trapping_day, y = plot_id)) +
  geom_point(aes(color = factor(session))) +
  coord_cartesian(xlim = c(1, 150))

## now let's join the trapping session data back to the capture data
mammal_captures_taxon <- left_join(mammal_captures_taxon, trapping_sessions, by = c("plot_id", "trapping_day"))

## let's make a data frame of only trapping sessions at TREE
tree_captures <- filter(mammal_captures_taxon, site_id == "TREE")

## now let's look at the captures with missing values for tag_id (AKA captures that escaped or are species that aren't tagged)
View(filter(tree_captures, is.na(tag_id)))

## there are none because we only have data for captures with tick info
## thus wee need to create a session column for our earlier mammal_trap_df data frame
## but first we need to create a trapping_day column in the data frame to join things
mammal_trap_df <- mammal_trap_df %>%
  group_by(plot_id) %>%
  mutate(trapping_date = as.Date(collect_date),
         trapping_day = as.numeric(trapping_date - min(trapping_date)))
## and recreate sessions with full mammal data
trapping_sessions <- select(mammal_trap_df, plot_id, trapping_day) %>%
  distinct() %>%
  arrange(plot_id, trapping_day) %>%
  mutate(plot_day = str_c(plot_id, "_", trapping_day),
         day_diff = c(0, diff(trapping_day)),
         session = cumsum(day_diff > 10),
         plot_session = str_c(plot_id, "_", session))

## now join in sessions
mammal_trap_df <- left_join(mammal_trap_df, trapping_sessions, by = c("plot_id", "trapping_day"))

## now let's filter down to just TREE data
tree_trap_df <- filter(mammal_trap_df, site_id == "TREE")

## now we want to summarize for every session the minimum number alive (total unique captures in a session)
## of each species, including species that are not present but were observed in other sessions
## we will also want to include the number of unique species in each session
## perhaps we can start by making a data frame with every combination of session and taxon_id

## step 1, change sessions and taxon_id to factors with mutate
tree_trap_df <- tree_trap_df %>%
  mutate(plot_session = factor(plot_session),
         taxon_id = factor(taxon_id))

## step 2, create a data frame with every combination of session and taxon_id
tree_mna_by_session <- tibble(plot_session = 
                                rep(levels(tree_trap_df$plot_session), each = length(levels(tree_trap_df$taxon_id))),
                              taxon_id = 
                                rep(levels(tree_trap_df$taxon_id), length(levels(tree_trap_df$plot_session)))
                              )

## step 3, calculate minimum number alive for each combination
## first make data frame of only captures that have tag_id values
tree_trap_df_tagged <- filter(tree_trap_df, !is.na(tag_id))
## then find only unique tag_ids by session
tree_unique <- tree_trap_df_tagged %>%
  distinct(tag_id, plot_session, taxon_id) %>%
  arrange(plot_session, taxon_id)

## now summarize the mna for each session and taxon_id
tree_mna_summary <- tree_unique %>%
  group_by(plot_session, taxon_id) %>%
  summarize(mna = n())

## now join this data into the tree_mna_by_session data frame, inputting zeros for when a taxon_id was not captured
tree_mna_by_session <- left_join(tree_mna_by_session, tree_mna_summary, by = c("plot_session", "taxon_id")) %>%
  mutate(mna = replace_na(mna, 0)) %>%
  group_by(plot_session) %>%
  mutate(total_mna = sum(mna),
         richness = sum(mna > 0),
         prop_mna = if_else(total_mna > 0, mna / total_mna, 0))

## let's do some quick visulaizations
## first, let's visualize the frequency of mna values for each of some select taxons
ggplot(filter(tree_mna_by_session, taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = mna, fill = taxon_id)) +
  geom_histogram() +
  facet_wrap(~taxon_id)

## now let's see how their prop_mna values vary with total_mna
ggplot(filter(tree_mna_by_session, taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = total_mna, y = prop_mna, color = taxon_id)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~taxon_id)

## including all the session with zero captures of a taxon throws things off...

## but also let's look at the same relationship with richness
ggplot(filter(tree_mna_by_session, taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = richness, y = prop_mna, color = taxon_id)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~taxon_id)

## similar problem
## let's try filtering out only the sessions where total_mna = 0
tree_mna_by_session_nonzero <- filter(tree_mna_by_session, total_mna > 0)

## now revisit the previous plots
ggplot(filter(tree_mna_by_session_nonzero, taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = mna, fill = taxon_id)) +
  geom_histogram() +
  facet_wrap(~taxon_id)

ggplot(filter(tree_mna_by_session_nonzero, taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = total_mna, y = prop_mna, color = taxon_id)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~taxon_id)

ggplot(filter(tree_mna_by_session_nonzero, taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = richness, y = prop_mna, color = taxon_id)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~taxon_id) +
  scale_x_continuous(breaks = seq(0, 13, 1)) +
  theme_bw()

## now let's modify the tree_mna_by_session data to evaluate plot differences
tree_mna_by_plot <- tree_mna_by_session %>%
  distinct(plot_session, total_mna, richness) %>%
  separate_wider_delim(cols = plot_session, names = c("site_id", "plot_num", "session"), 
                       delim = "_", cols_remove = FALSE) %>%
  mutate(session = as.numeric(session),
         plot_id = as.factor(str_c(site_id, "_", plot_num)))

## histogram of total_mna by plot
ggplot(tree_mna_by_plot, aes(x = total_mna, fill = plot_id)) +
  geom_histogram() +
  facet_wrap(~plot_id)

ggplot(tree_mna_by_plot, aes(x = richness, fill = plot_id)) +
  geom_histogram() +
  facet_wrap(~plot_id)

## let's look specifically at TREE_028 with all the mna data
tree_mna_by_session <- tree_mna_by_session %>%
  separate_wider_delim(cols = plot_session, names = c("site_id", "plot_num", "session"), 
                       delim = "_", cols_remove = FALSE) %>%
  mutate(session = as.numeric(session),
         plot_id = as.factor(str_c(site_id, "_", plot_num)))

tree_mna_by_session_nonzero <- tree_mna_by_session_nonzero %>%
  separate_wider_delim(cols = plot_session, names = c("site_id", "plot_num", "session"), 
                       delim = "_", cols_remove = FALSE) %>%
  mutate(session = as.numeric(session),
         plot_id = as.factor(str_c(site_id, "_", plot_num)))

ggplot(filter(tree_mna_by_session_nonzero, plot_id == "TREE_028", 
              taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = richness, y = prop_mna, color = taxon_id)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~taxon_id)

ggplot(filter(tree_mna_by_session_nonzero, plot_id == "TREE_027", 
              taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = richness, y = prop_mna, color = taxon_id)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~taxon_id)

ggplot(filter(tree_mna_by_session_nonzero, plot_id == "TREE_016", 
              taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = total_mna, y = prop_mna, color = taxon_id)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~taxon_id)

ggplot(filter(tree_mna_by_session_nonzero, plot_id == "TREE_019", 
              taxon_id %in% c("PELE","PEMA","TAST","NAIN","ZAHU","BLBR","SOCI","MYGA")), 
       aes(x = total_mna, y = prop_mna, color = taxon_id)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~taxon_id)




###
## old and outdated code below

## now filter down to just captures from TREE
tree_captures <- filter(mammal_trap_df, trap_status_code %in% c("4","5"), site_id == "TREE")

## checking captures with no tag_id
View(filter(tree_captures, is.na(tag_id)))

## let's filter out fate = "escaped", as it is not possible to determine if they are recaptures from the same session
tree_captures <- filter(tree_captures, fate != "escaped")

## actually, let's just filter out all of the captures with missing tag_id
tree_captures <- filter(tree_captures, !is.na(tag_id))

## now let's make a new df called "tree_mna" that feature only unique tag_id values for each trapping session
tree_mna <- tree_captures %>%
  distinct(tag_id, plot_session)

## now let's join in the rest of the mammal data into this new data frame
tree_mna <- left_join(tree_mna, tree_captures, by = c("tag_id", "plot_session"))

## adding taxon grouping to tree_mna
tree_mna <- tree_mna %>%
  mutate(taxon_group = case_when(taxon_id %in% c("PEMA", "PELE", "PESP", "PELEPEMA", "TAST", "ZAHU", "NAIN", "MYGA") ~ taxon_id,
                                 TRUE ~ "Other"))

## let's visualize the number of captures of each taxon group by session
ggplot(tree_mna, aes(x = plot_session, fill = taxon_group)) +
  geom_bar() +
  coord_flip()

## let's summarize this into a summary df
tree_mna_summary <-  tree_mna %>%
  group_by(plot_session, taxon_group) %>%
  summarize(mna = n()) %>%
  mutate(total_mna = sum(mna),
         richness = n(),
         prop_mna = mna / total_mna) %>%
  separate(col = plot_session, into = c("site_id", "plot_num", "session"), sep = "_") %>%
  mutate(session = as.numeric(session),
         plot_id = str_c(site_id, "_", plot_num),
         session_id = str_c(plot_id, "_", session)) %>% 
  arrange(total_mna)

tree_ranks <- distinct(tree_mna_summary, total_mna, session_id) %>%
  mutate(mna_rank = row_number())

tree_mna_summary <- left_join(tree_mna_summary, tree_ranks, by = c("total_mna", "session_id"))

## let's visualize the mna of each taxon group by session
ggplot(tree_mna_summary, aes(x = mna_rank, y = mna, fill = taxon_group)) +
  geom_bar(stat = "identity", position = "stack")

## visualize how the number of taxon groups varies with total mna
ggplot(tree_mna_summary, aes(x = total_mna, y = mna, color = taxon_group)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

ggplot(tree_mna_summary, aes(x = richness, y = prop_mna, color = taxon_group)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

ggplot(filter(tree_mna_summary, taxon_group == "PELE"), aes(x = total_mna, y = prop_mna)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

ggplot(filter(tree_mna_summary, taxon_group == "PELE"), aes(x = richness, y = prop_mna)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)





