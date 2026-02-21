##
#' file: neon_data_import
#' description: download neon data for use in small mammal biodiversity and tick analysis
#' using neonUtilities package
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

## iteratively download data for each data product and site (testing on small date range)
map2(.x = rep(data_products$dpID, nrow(sites)), .y = rep(sites$siteCode, each = nrow(data_products)), 
     ~zipsByProduct(site = .y, dpID = .x, savepath = paste0("raw_data/", .y),
                    package = "basic", startdate = "2021-01-01", enddate = "2021-12-31", 
                    check.size = FALSE))

## the above ends when no data is found, and also saves files in non-descriptive folders
## need to write a function that checks if the files are there, downloads them if so, 
## renames the "filesToStack#####" folder to the data product name, and moves to the next site/data product if no data is found

## function to download data for a site and data product, does not end if no data is found, 
## and renames filesToStack folder to data product name and site code
download_data <- function(site, dpID, startdate, enddate){
  ## download data
  zipsByProduct(site = site, dpID = dpID, savepath = paste0("raw_data/", site),
                package = "basic", startdate = startdate, enddate = enddate, check.size = FALSE)
  ## check if data was downloaded
  if(file.exists(paste0("raw_data/", site, "/filesToStack", dpID))){
    ## rename folder
    file.rename(paste0("raw_data/", site, "/filesToStack", dpID), 
                paste0("raw_data/", site, "/", data_products$data_name[data_products$dpID == dpID]))
  }
}
## the above is pretty close, but I think "zipsByProduct" will still error out when no data is found
## so I need to wrap it in a tryCatch statement
## now the same function as above but with zipByProduct wrapped in tryCatch
download_data <- function(site, dpID, startdate, enddate){
  ## download data
  tryCatch({
    zipsByProduct(site = site, dpID = dpID, savepath = paste0("raw_data/", site),
                  package = "basic", startdate = startdate, enddate = enddate, check.size = FALSE)
  }, error = function(e) {
    message("No data found for ", data_products$data_name[data_products$dpID == dpID])
  })
  ## check if data was downloaded
  if(file.exists(paste0("raw_data/", site, "/filesToStack", dpID))){
    ## rename folder
    file.rename(paste0("raw_data/", site, "/filesToStack", dpID), 
                paste0("raw_data/", site, "/", data_products$data_name[data_products$dpID == dpID]))
  }
}
## this should work, but I need to test it
## I will test it on a single site and data product
download_data(site = "BART", dpID = data_products$dpID[1], startdate = "2022-01-01", enddate = "2022-12-31")
## this did not work because the dpID features a "DP1." prefix and a ".##" suffix,
## while the "filesToStack" folder does not, thus I need to use stringr to cut out the prefix and suffix
## I will also add a print statement to see if the function is working
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
                paste0("raw_data/", site, "/", data_products$data_name[data_products$dpID == dpID]))
    print(paste("Data downloaded for", data_products$data_name[data_products$dpID == dpID], "at", site))
  }
}
## test function
download_data(site = "BART", dpID = data_products$dpID[1], startdate = "2022-01-01", enddate = "2022-12-31")
## now let's try for data that does not exist
download_data(site = "BART", dpID = data_products$dpID[3], startdate = "2022-01-01", enddate = "2022-12-31")
## seems to work, but I want to add site and date range to the renamed folder name
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
## test function
download_data(site = "BART", dpID = data_products$dpID[1], startdate = "2022-01-01", enddate = "2022-12-31")
download_data(site = "BART", dpID = data_products$dpID[3], startdate = "2022-01-01", enddate = "2022-12-31")

## now I will run the function for all sites and data products
map2(.x = rep(data_products$dpID, nrow(sites)), .y = rep(sites$siteCode, each = nrow(data_products)), 
     ~download_data(site = .y, dpID = .x, startdate = "2022-01-01", enddate = "2022-12-31"))
## seemed to work great, now I can consider if the folder hierarchy and naming scheme are appropriate
## next steps: clean up this script to be more streamlined
## and then write a script to stack the data