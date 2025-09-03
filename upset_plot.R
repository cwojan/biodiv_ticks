##
# upset plot of mammal communities
##

library(tidyverse)
library(UpSetR)
library(ComplexUpset)
library(DescTools)

# Load the data
data_date <- "2025-09-03"
mammals <- read_csv(str_c("processed_data/mammal_community_effects_", data_date,".csv"))

# select only site, mammal and presence columns
mammals <- mammals %>%
  select(site_id, plot_id, nlcd_class, plot_session, session,
         taxon_id, year, mean_month, mean_yday, presence)

## check most common mammals
mamm_sum <- mammals %>%
  group_by(taxon_id) %>%
  summarise(sum = sum(presence)) %>%
  arrange(desc(sum))

mamm_filt <- mammals %>%
  filter(taxon_id %in% c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN"))

# create a wide format data frame for the upset plot
mamm_wide <- mamm_filt %>%
  pivot_wider(names_from = taxon_id, values_from = presence) %>%
  mutate(total = PEMA + PELE + MYGA + BLBR + TAST + NAIN,
         pele_pres = if_else(PELE == 1, TRUE, FALSE)) %>%
  filter(total > 0) %>%
  data.frame() %>%
  mutate(region = case_when(site_id == "BART" ~ "New Hampshire",
                            site_id == "HARV" ~ "Massachusetts",
                            site_id == "BLAN" ~ "Virginia A",
                            site_id == "SCBI" ~ "Virginia B",
                            site_id == "SERC" ~ "Maryland",
                            site_id == "STEI" ~ "Northwest WI",
                            site_id == "TREE" ~ "Northeast WI",
                            site_id == "UNDE" ~ "Upper MI"),
         region = factor(region, levels = c("Massachusetts", "New Hampshire",
                                         "Virginia A", "Virginia B", "Maryland",
                                         "Upper MI", "Northeast WI", "Northwest WI"))) 

mamm_labels <- c(
  PELE = "White-footed\nMouse",
  PEMA = "Deer Mouse",
  BLBR = "Short-tailed\nShrew",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)

mamms <- c("PELE", "PEMA", "BLBR", "MYGA", "TAST", "NAIN")
mamm_combs <- CombSet(mamms, m = 1:6)
comb_lists <- map(.x = mamm_combs, .f = function(x){
  split(x, rep(1:nrow(x)))
}) %>%
  unlist(recursive = FALSE)

pele_combs <- comb_lists[map_lgl(comb_lists, ~ "PELE" %in% .x)]
non_pele_combs <- comb_lists[!map_lgl(comb_lists, ~ "PELE" %in% .x)]
sort_list <- c(pele_combs, non_pele_combs)
sort_list <- sort_list[-c(15,16,26)]

query_list <- map(.x = comb_lists, .f = function(x){
  if ("PELE" %in% x){
    upset_query(intersect = as.vector(x), color = "#E66101",
                only_components = "intersections_matrix")
  } else {
    upset_query(intersect = as.vector(x), color = "#008080",
                only_components = "intersections_matrix")
  }
})
query_list

query_list[[1]]

summary(mamm_wide$NAIN)
summary(c(1,2,3,NA))

# create upset plots

upset(
  data = mamm_wide, intersect = c("PEMA", "PELE", "MYGA", "BLBR", "TAST", "NAIN"),
  name = "Community Composition (by trap session)",
  base_annotations = list(
    "Intersection size" = intersection_size(counts = FALSE,
                                            mapping = aes(fill = pele_pres)) +
      scale_fill_manual(values = c("#008080", "#E66101"),
                        name = "White-footed Mouse\nPresence",
                        labels = c("Absent", "Present"))
  ),
  set_sizes = (
    upset_set_size(geom = geom_bar(aes(fill = pele_pres))) +
      scale_fill_manual(values = c("#008080", "#E66101"),
                        name = "White-footed Mouse\nPresence",
                        labels = c("Absent", "Present"))
  ),
  labeller = as_labeller(mamm_labels),
  themes = upset_modify_themes(
    list(
      "Intersection size" = theme(
        legend.position = c(-0.15,0.2),
        axis.title = element_blank()
      ),
      "overall_sizes" = theme(
        axis.title = element_blank(),
        legend.position = "none"
      )
    )
  ),
  queries = query_list,
) +
  theme(
    text = element_text(size = 16),
    axis.text.y = element_text(size = 12)
  )



upset(
  data = mamm_wide, intersect = c("PEMA", "PELE", "MYGA", "BLBR", "TAST", "NAIN"),
  name = "Community Composition (by trap session)",
  base_annotations = list(),
  set_sizes = FALSE,
  annotations = list(
    "Day of Year" = (
      ggplot(mapping = aes(y = mean_yday)) +
        geom_jitter(aes(color = pele_pres), width = 0.2, alpha = 0.5) +
        geom_boxplot(aes(fill = pele_pres), alpha = 0.2) +
        scale_color_manual(values = c("#008080", "#E66101"),
                           name = "White-footed Mouse\nPresence",
                           labels = c("Absent", "Present")) +
        scale_fill_manual(values = c("#008080", "#E66101"),
                          name = "White-footed Mouse\nPresence",
                          labels = c("Absent", "Present")) +
        ylab("Day of Year") +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
        )
    )
  ),
  labeller = as_labeller(mamm_labels),
  themes = upset_modify_themes(
    list(
      "Intersection size" = theme(
        axis.title = element_blank()
      ),
      "overall_sizes" = theme(
        axis.title = element_blank(),
        legend.position = "none"
      )
    )
  ),
  queries = query_list,
) +
  theme(
    text = element_text(size = 16),
    axis.text.y = element_text(size = 12)
  )


upset(
  data = mamm_wide, intersect = c("PEMA", "PELE", "MYGA", "BLBR", "TAST", "NAIN"),
  name = "Community Composition (by trap session)",
  base_annotations = list(),
  set_sizes = FALSE,
  annotations = list(
    "NLCD Class" = (
      ggplot() +
        geom_bar(aes(fill = nlcd_class), position = "fill") +
        geom_vline(xintercept = 29.5) +
        scale_fill_viridis_d(name = "NLCD Class",
                             labels = c("Cultivated Crops", "Deciduous Forest",
                                        "Evergreen Forest", "Mixed Forest",
                                        "Pasture Hay", "Woody Wetlands")) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
        )
    )
  ),
  labeller = as_labeller(mamm_labels),
  themes = upset_modify_themes(
    list(
      "Intersection size" = theme(
        axis.title = element_blank()
      ),
      "overall_sizes" = theme(
        axis.title = element_blank(),
        legend.position = "none"
      )
    )
  ),
  queries = query_list,
  sort_intersections = FALSE,
  intersections = sort_list
) +
  geom_vline(xintercept = 29.5) +
  theme(
    text = element_text(size = 16),
    axis.text.y = element_text(size = 12)
  )

upset(
  data = mamm_wide, intersect = c("PEMA", "PELE", "MYGA", "BLBR", "TAST", "NAIN"),
  name = "Community Composition (by trap session)",
  base_annotations = list(),
  set_sizes = FALSE,
  annotations = list(
    "Region" = (
      ggplot() +
        geom_bar(aes(fill = region), position = "fill") +
        geom_vline(xintercept = 29.5) +
        scale_fill_viridis_d(name = "Region", option = "turbo") +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
        )
    )
  ),
  labeller = as_labeller(mamm_labels),
  themes = upset_modify_themes(
    list(
      "Intersection size" = theme(
        axis.title = element_blank()
      ),
      "overall_sizes" = theme(
        axis.title = element_blank(),
        legend.position = "none"
      )
    )
  ),
  queries = query_list,
  sort_intersections = FALSE,
  intersections = sort_list
) +
  geom_vline(xintercept = 29.5) +
  theme(
    text = element_text(size = 16),
    axis.text.y = element_text(size = 12)
  )





