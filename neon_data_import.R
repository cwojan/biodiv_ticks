##
#' file: neon_data_import
#' description: download neon data for use in small mammal biodiversity and tick analysis
#' using neonUtilities package
#' Note: this can't overwrite data currently
##

## load packages
library(neonUtilities)
library(tidyverse)
library(janitor)

## load data product ids
data_products <- read_csv("neon_data_ids.csv")
## reformat data product names
data_products$data_name <- make_clean_names(data_products$data_category)
## load sites of interest
sites <- read_csv("neon_sites.csv")

## create folder for each site code with map
map(sites$siteCode, ~dir.create(paste0("raw_data/", .x)))

## write function to download data
download_data <- function(site, dpID, startdate, enddate){
  ## download data
  tryCatch({
    zipsByProduct(site = site, dpID = dpID, savepath = paste0("raw_data/", site),
                  package = "basic", startdate = startdate, enddate = enddate, check.size = FALSE)
  }, error = function(e) {
    message("No data found for ", data_products$data_name[data_products$dpID == dpID])
  })
  ## check if data was downloaded
  if(file.exists(paste0("raw_data/", site, "/filesToStack", str_remove_all(dpID, "DP1\\.|\\.\\d+")))){
    ## rename folder
    file.rename(paste0("raw_data/", site, "/filesToStack", str_remove_all(dpID, "DP1\\.|\\.\\d+")), 
                paste0("raw_data/", site, "/", data_products$data_name[data_products$dpID == dpID], "_", site, "_", startdate, "_", enddate))
    print(paste("Data downloaded for", data_products$data_name[data_products$dpID == dpID], "at", site))
  }
}

## the above function works, but let's add stacking the downloaded data with stackByTable, and renaming the folder of stacked data
stack_data <- function(site, dpID, startdate, enddate){
  ## download data
  tryCatch({
    zipsByProduct(site = site, dpID = dpID, savepath = paste0("raw_data/", site),
                  package = "basic", startdate = startdate, enddate = enddate, check.size = FALSE)
  }, error = function(e) {
    message("No data found for ", data_products$data_name[data_products$dpID == dpID])
  })
  ## check if data was downloaded
  if(file.exists(paste0("raw_data/", site, "/filesToStack", str_remove_all(dpID, "DP1\\.|\\.\\d+")))){
    ## stack data
    stackByTable(filepath = paste0("raw_data/", site, "/filesToStack", str_remove_all(dpID, "DP1\\.|\\.\\d+")))
    ## copy all files from "stackedFiles" subdirectory to the parent folder
    file.copy(from = 
                list.files(paste0("raw_data/", site, "/filesToStack", str_remove_all(dpID, "DP1\\.|\\.\\d+"), "/stackedFiles"), 
                           full.names = TRUE), 
              to = paste0("raw_data/", site, "/filesToStack", str_remove_all(dpID, "DP1\\.|\\.\\d+")))
    ## remove "stackedFiles" subdirectory
    unlink(paste0("raw_data/", site, "/filesToStack", str_remove_all(dpID, "DP1\\.|\\.\\d+"), "/stackedFiles"), recursive = TRUE)
    ## rename folder
    file.rename(paste0("raw_data/", site, "/filesToStack", str_remove_all(dpID, "DP1\\.|\\.\\d+")), 
                paste0("raw_data/", site, "/", site, "_", data_products$data_name[data_products$dpID == dpID],  
                       "_", startdate, "_", enddate))
    print(paste("Data downloaded for", data_products$data_name[data_products$dpID == dpID], "at", site))
  }
}

## download data for each site
map2(.x = rep(data_products$dpID, nrow(sites)), .y = rep(sites$siteCode, each = nrow(data_products)), 
     ~stack_data(site = .y, dpID = .x, startdate = "2012-01-01", enddate = "2024-12-31"))

## create txt file in raw_data folder recording date of download
write(paste("Data downloaded on", Sys.Date()), file = "raw_data/download_date.txt")
