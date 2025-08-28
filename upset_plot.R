##
# upset plot of mammal communities
##

library(tidyverse)
library(UpSetR)
library(ComplexUpset)
library(DescTools)

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
  filter(taxon_id %in% c("PEMA", "PELE", "MYGA", "BLBR", "TAST", "NAIN"))

# create a wide format data frame for the upset plot
mamm_wide <- mamm_filt %>%
  pivot_wider(names_from = taxon_id, values_from = presence) %>%
  mutate(total = PEMA + PELE + MYGA + BLBR + TAST + NAIN,
         pele_pres = if_else(PELE == 1, TRUE, FALSE)) %>%
  filter(total > 0) %>%
  data.frame()

mamm_labels <- c(
  PEMA = "Deer Mouse",
  PELE = "White-footed\nMouse",
  MYGA = "Red-backed\nVole",
  BLBR = "Short-tailed\nShrew",
  TAST = "Eastern\nChipmunk",
  NAIN = "W. Jumping\nMouse"
)

mamms <- c("PEMA", "PELE", "MYGA", "BLBR", "TAST", "NAIN")
mamm_combs <- CombSet(mamms, m = 1:6)
comb_lists <- map(.x = mamm_combs, .f = function(x){
  split(x, rep(1:nrow(x)))
}) %>%
  unlist(recursive = FALSE)

query_list <- map(.x = comb_lists, .f = function(x){
  if ("PELE" %in% x){
    upset_query(intersect = as.vector(x), color = "blue",
                only_components = "intersections_matrix")
  } else {
    upset_query(intersect = as.vector(x), color = "black",
                only_components = "intersections_matrix")
  }
})
query_list

query_list[[1]]

# create upset plot
upset(
  data = mamm_wide, intersect = c("PEMA", "PELE", "MYGA", "BLBR", "TAST", "NAIN"),
  name = "Community Composition (by trap session)",
  base_annotations = list(
    "Intersection size" = intersection_size(counts = FALSE,
                                            mapping = aes(fill = pele_pres))
  ),
  set_sizes = (
    upset_set_size(geom = geom_bar(aes(fill = pele_pres)))
  ),
  labeller = as_labeller(mamm_labels),
  themes = upset_modify_themes(
    list(
      "Intersection size" = theme(
        axis.title = element_blank()
      ),
      "overall_sizes" = theme(
        axis.title = element_blank()
      )
    )
  ),
  queries = query_list,
  guides = "collect"
) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 16),
    axis.text.y = element_text(size = 12)
  )

