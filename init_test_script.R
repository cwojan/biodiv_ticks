## init script for biodiversity tick project

## load libraries
library(tidyverse)
library(neonUtilities)
library(neonOS)


## load by product for quick import, redownloads with each run
ticks_test <- loadByProduct(dpID = "DP1.10093.001", site = "TREE",
                            startdate = "2017-06", enddate = "2017-07")

tck_fielddata <- ticks_test$tck_fielddata
