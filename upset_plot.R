##
# upset plot of mammal communities
##

library(tidyverse)
library(UpSetR)
library(ComplexUpset)

# Load the data
mammals <- read_csv("processed_data/mammal_community_effects.csv")

# select only site, mammal and presence columns
mammals <- mammals %>%
  select(site_id, plot_id, nlcd_class, plot_session, session,
         taxon_id, presence)

## check most common mammals
mamm_sum <- mammals %>%
  group_by(taxon_id) %>%
  summarise(sum = sum(presence)) %>%
  arrange(desc(sum))

mamm_filt <- mammals %>%
  filter(taxon_id %in% c("PEMA", "PELE", "MYGA", "TAST", "NAIN"))

# create a wide format data frame for the upset plot
mamm_wide <- mamm_filt %>%
  pivot_wider(names_from = taxon_id, values_from = presence) %>%
  data.frame()

mamm_labels <- c(
  PEMA = "Deer Mouse",
  PELE = "White-footed\nMouse",
  MYGA = "Red-backed\nVole",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)

# create upset plot
upset(
  data = mamm_wide, intersect = c("PEMA", "PELE", "MYGA", "TAST", "NAIN"),
  name = "Community Composition (by trap session)",
  group_by = "sets",
  base_annotations = list(
    "Intersection size" = intersection_size(counts = FALSE)
  ),
  labeller = as_labeller(mamm_labels),
  queries = list(
    upset_query(
      set = "PELE",
      fill = "blue",
      name = "Highest Quality Disease Host"
    ),
    upset_query(
      group = "PELE",
      fill = "blue", color = "blue",
      name = "Highest Quality Disease Host"
    )
  )
)
