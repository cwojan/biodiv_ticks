# Evaluating Small Mammal Community Assembly and Tick Attachment with NEON Data

## Broad Overview

This repository contains code for importing and analyzing small mammal trapping data from the National Ecological Observatory Network (NEON). It accompanies the article \_\_\_\_\_.

The repository is composed of two data processing scripts, three analysis scripts, and two supplementary scripts.

1.  bdt_neon_data_import.R uses the neonUtilities packages to download small mammal trapping data from the NEON data portal and save it by site into folders inside a "raw_data" directory.
2.  bdt_mamm_trap_data_process.R wrangles the small mammal data into distinct trapping sessions, adds useful variables, and calculates species richness, population sizes, and community composition.
3.  bdt_presence_analyses.R uses bootstrapping to assess how often species are present in trapping sessions of varying species richness compares to the expectation given neutral community assembly.
4.  bdt_abundance_analyses.R uses bootstrapping to assess how the relative abundance of species in trapping sessions of various species richness compares to the expectation given neutral community assembly.
5.  bdt_tick_attach_analyses.R uses bootstrapping to assess how ticks are distributed among mammals compared to what one would expect if ticks attached to mammals randomly.
6.  bdt_site_map.R generates a map of the NEON sites from which data was analyzed.
7.  bdt_supp_mna_richness.R draw a supplementary figure showing the relationship between mammal community size and species richness.

The "old_misc_files" folder contains developmental / exploratory scripts. The "logistics_data" folder contains a pair of small .csv files with site and data info that are used for data importing and a few other things.

## Data Import and Processing

### Import

The bdt_neon_data_import.R script will first load the logistic data, i.e. the names of the sites of interest in the Upper Midwest and Northeast, and the ID values fro the data products of interest. The script downloads mammal trapping data, mammal pathogen data, mammal DNA barcode data, tick colleciton data, and tick pathogen data. However, only the small mammal trapping data is used in analyses.

The script creates a directory in a "raw_data" folder for each site of interest, named after its four-letter site code. With those directories, the script creates and iterates a function that checks for data for each site and data product ID and downloads it with neonUtilities::zipsByProduct. The function then combines data from different time periods with neonUtilities::stackByTable, and organizes the data into meaningfully named folders.

The script ends by recording the date of download into a .txt file.

### Processing

The bdt_mamm_trap_data_process.R script iteratively reads in the downloaded mammal trapping data across the sites of interest and creates three processed data frames for analysis. The mammal trapping data is set up such that each row represents a single trap on a single trapping night.

First, it reformats existing columns and adds some helpful columns to the trapping data. Notably, it creates a year, month, and day of year column out of the trapping date. It also identifies captures where the animal was identified to species, and whether it was tagged.

Next, it separates trapping nights into distinct sessions, i.e., consecutive or semi-consecutive nights. It does this by calculating the difference in date for each plot's trapping nights, and then cumulatively summing the number of "breaks", or date differences \> 10 in each plot's sequence. The session IDs are merged into the main trapping data and saved as an .rds file - mammal_session_df.

Next, we reconfigure the data into a form where each plot and session has a row for each species. We calculate the minimum number alive (MNA) for each species and each session with 1.) the maximum number observed in a night for species that are not given permanent tag; and 2.) the total unique tags observed across nights of a session for species that are given permanent tags. These MNA values are then used to determine species presence, species richness, total MNA (all species), and relative abundance. The data frame of each species for each plot/session is saved as an .rds file - mammal_community_df.

Finally, for only the 6 most common mammals, tick attachment data is filtered out. We take each animal with complete tick attachment data in a session and consider it to have ticks if any tick life stage was attached to it during any night in a session. The resulting data frame features a row for each identifiable mammal with tick info, and columns on species, plot, session, tick attachment etc. It is saved as an .rds file - mammal_tick_captures_df.

## Data Analysis

The analyses are all performed separately for sites in the Upper Midwest (NEON Domain 05) and those in the Northeast (NEON Domains 01 and 02).

### Presence Analyses

The first analysis script, bdt_presence_analyses.R, uses the mammal_community_df .rds file to examine the presence of the 6 most common mammal species across sessions.

It starts by creating species pools - lists of species that are observed at a given site and in a given National Land Cover Database class. There are multiple plots in a site with the same NLCD class, so this can pull from multiple plots in a site. It aims to assess which species would be reasonably found in an area.

For a given plot and session, the sample_species function draws species from the relevant species pool until the species richness value for that session is reached (without replacement). For example, if a plot's relevant species pool has 10 species, and a given session features 5 species, the function would randomly draw 5 of those 10 species. This process is replicated a default of 1000 times to generate many "simulated" communities.

The sampling function is run for each plot and session, and then filtered down to only the 6 most common species. The observed data is filtered such that it only includes rows where the given species is present in the relevant species pool.

For both the observed and simulated communities, binomial regressions of presence as a function of species richness are run for each species and region (Upper Midwest sites and Northeast sites). The coefficients and significance are not of interest, just the comparison of predictions between observed and simulated models are. The non-independence due to site, year, etc. is accounted for because the observed and simulated data have the same non-independence structure.

The observed and simulated model prediction are all plotted together for each species and region.

### Abundance Analyses

The "bdt_abundance_analyses.R" is very similar to the presence analysis script. However, the bootstrapping function simulates full communities instead of just species assemblages. The "sample_community" function creates communities for plots and sessions that have the same number of total animals, and of the same species that were there. Each individual in the simulated community is drawn from the species that were there, with replacement. Thus, a simulated community will have different proportions of species, generally close to equal proportions among species but with stochasticity.

For each of the 6 most common species, binomial models are run for the proportion of the community the species makes up as a function of species richness, with the total number of animals as the weight for each observation. E.g., if there are 10 animals, and 8 of them were deer mice, that translates to 8 out of 10 "successes" for the binomial model. Predictions from models based on observed proportions and simulated proportions are compared visually in a ggplot.

### Tick Attachment Analyses

The "bdt_tick_attach_analyses.R" script examines tick attachment data from multiple angles. But first it summarizes the mammal_tick_captures_df (only mammals with tick data) by finding the number and proportion of each species that have ticks in each session. This is what can be modeled and compared.

The first approach is to simply model the number of animal with tick by taxon and region as a function of species richness, using a negative binomial regression with site, year, and month as random effects (NLCD class as a fixed effect). This is summarized in a table.

The next approaches require simulation / bootstrapping. The presence of ticks on an animal is shuffled by individual 1000 times, to create 100 simulated replicates where ticks have assembled randomly on the mammals present in a session.

Then, we find the 95% quantiles of the proportions of species with ticks and compare the observed proportions to those quantiles. We summarize how many are outside the 95% range, and visualize them with ggplot by species and region, across species richness.

Finally, we calculate the mean and standard deviations of the proportions of mammals with ticks from the simulated data. We then use those values to calculate how many standard deviations each observed proportion is from the simulated mean (i.e., standardized differences), and visualize by species and region with ggplot.

## Supplementary Code

### Site Map

The "bdt_site_map.R" script is very straightforward. It loads the site locations from the site info .csv file in the logistics_data folder. It uses the rnaturalearth::ne_download() to load shape files of the US state boundaries and the boundaries of countries, in "sf" package format. The site coordinates are transformed into sf format with the projection of the rnaturalearth shape files. The three files are mapped with ggplot, and a north arrow and scale bar are added from the "ggspatial" package.

The site coordinates can be found here:

<https://www.neonscience.org/field-site-map-and-info>

### Supplementary Figure

The "bdt_supp_mna_rich.R" script is also straightforward. It loads the mammal_sommunity_df that was created in the processing script, and pares it down to just one row for each plot and trapping session, with columns for total mammals observed (sums of the minimum number alive of each species) and species richness. Then ggplot is used to plot a figure showing how the total mammal count increases with species richness, in separate panels for the Upper Midwest and Northeast sites.
