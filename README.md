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

## Data Import and Processing

## Data Analysis

## Supplementary Code
