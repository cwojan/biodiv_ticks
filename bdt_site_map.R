##
#' file: bdt_site_map.R
#' author: chris wojan
#' description: generate map of the sites used for data analysis
##

## load libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggspatial)

## load NEON spatial data
sites <- read_csv("spatial_data/neon_sites.csv")

## download Natural Earth data for states at 1:110m scale
states110 <- ne_download(scale = 110, type = "admin_1_states_provinces_lakes", 
                         category = "cultural", returnclass = "sf")

## download Natural Earth data for countries at 1:110m scale
countries <- ne_download(type = "admin_0_countries_lakes", 
                         category = "cultural", scale = 110, returnclass = "sf")

## give the sites a region label based on their domain code
sites_sf <- sites %>%
  mutate(region = case_when(domainCode %in% c("D01", "D02") ~ "Northeast",
                            domainCode == "D05" ~ "Upper Midwest",
                            .default = "Other")) %>%
  select(region, siteCode, latitude, longitude) %>%
  #c convert to sf object
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(states110)) 


## draw map with state outlines, canada, and the sites colored by region
ggplot() +
  geom_sf(data = countries, fill = "gray90", color = "black") +
  geom_sf(data = states110, fill = "gray80", color = "black") +
  geom_sf(data = sites_sf, color = "black", aes(fill = region),
          shape = 21, stroke = 1, size = 4, alpha = 0.7) +
  coord_sf(crs = 3857, xlim = c(-10500000, -7550000), ylim = c(4300000, 6300000)) +
  scale_fill_viridis_d(name = "Region") +
  annotation_scale(location = "br", width_hint = 0.2) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  labs(caption = "Site points from National Ecological Observatory Network (NEON) \n
       Basemap made with Natural Earth. Free vector and raster map data @ naturalearthdata.com.") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.15, 0.15),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.background = element_rect(fill = "lightblue", color = NA))

  
    



