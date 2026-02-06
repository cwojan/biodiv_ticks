## 
# site map
##

## load libraries
library(tidyverse)
library(sf)
library(ggspatial)
library(ggmap)


## load data
domains <- read_sf("spatial_data/NEONDomains_2024", layer = "NEON_Domains")
sites <- read_csv("spatial_data/neon_sites.csv")


domain_subset <- domains %>%
  filter(domainID %in% c("D01", "D02", "D05")) %>%
  mutate(domainID = factor(domainID, levels = c("D01", "D02", "D05")),
         domainName = recode(domainID,
                             "D01" = "Domain 01 (Northeast)",
                             "D02" = "Domain 02 (Northeast)",
                             "D05" = "Domain 05 (Upper Midwest)"))


sites_sf <- sites %>%
  select(siteCode, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

domain_proj <- st_transform(domain_subset, crs = 4326)

## ggplot
ggplot() +
  annotation_map_tile(type = "osm", zoom = 6) +
  geom_sf(data = domain_proj, fill = NA, linewidth = 0.5) +
  geom_sf_label(data = domain_proj, aes(label = domainName), nudge_y = c(50000, 1, -0.3)) +
  geom_sf(data = sites_sf, color = "black", fill = "blue", shape = 21, stroke = 1, size = 4, alpha = 0.7) +
  coord_sf(crs = 3857, xlim = c(-10500000, -7550000), ylim = c(4300000, 6300000)) +
  labs(caption = "Domain outlines and site points from National Ecological Observatory Network (NEON) \n
       Basemap © OpenStreetMap contributors") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank())

ggplot() +
  get_stadiamap() +
  geom_sf(data = domain_proj, fill = NA, linewidth = 0.5)
  geom_sf_label(data = domain_proj, aes(label = domainName), nudge_y = c(-0.8, 1, -0.3)) +
  geom_sf(data = sites_sf, color = "black", fill = "blue", shape = 21, 
          stroke = 1, size = 4, alpha = 0.7) +
  coord_sf(crs = 4326) +
  theme_minimal()
  
    



