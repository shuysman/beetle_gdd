library(terra)
library(lubridate)
library(foreach)
library(doParallel)
library(parallel)
library(iterators)
library(glue)
library(tidyverse)

beetle_year <- function(x) {
    index <- ymd(x)
    index <- as.character(year(add_with_rollback(index, months(start_month))))
    prename <- "y_"
    return(paste0(prename, index))
}

get_gdd <- function(t_k, tbase = tbase_k) {
    return(pmax(t_k - tbase, 0))
}

cores <- detectCores() - 2

terraOptions(verbose = TRUE,
             memfrac = 0.9)

tbase_k <- 5.5 + 273.15

models <- c(
    "bcc-csm1-1",
    "BNU-ESM",
    "CanESM2",
    "CNRM-CM5",
    "CSIRO-Mk3-6-0",
    "GFDL-ESM2G",
    "GFDL-ESM2M",
    "HadGEM2-CC365",
    "HadGEM2-ES365",
    "inmcm4",
    "IPSL-CM5A-LR",
    "IPSL-CM5A-MR",
    "IPSL-CM5B-LR",
    "MIROC-ESM",
    "MIROC-ESM-CHEM",
    "MIROC5",
    "MRI-CGCM3"
 ##   "NorESM1-M"  ## Missing rcp85 scenario
)

scenarios <- c("rcp45", "rcp85")

data_dir <- file.path("~/data/MACA/gye/forecasts/daily")
out_dir <- file.path("~/out")

## t_file <- file.path(data_dir, "tavg_inmcm4_rcp45_2006-2099_daily_gye.nc")

## r <- rast(t_file)

start_month <- 5

mclapply(models,
    function(model) {
        for (scenario in scenarios) {
            print(model)
            print(scenario)
            in_filename <- glue("tavg_{model}_{scenario}_2006-2099_daily_gye.nc")
            r <- terra::rast(file.path(data_dir, in_filename))
            
            gdd <- terra::app(r, fun = get_gdd, tbase = tbase_k)

            time(gdd) <- time(r)

            accumgdd <- terra::tapp(gdd, index = beetle_year, fun = "sum")

            ## cut off first and last layers in spatraster, they are incomplete because of "beetle year" indexing
            accumgdd2 <- subset(accumgdd, seq(2, nlyr(accumgdd) - 1))
            
            time(accumgdd2, tstep = "years") <- 2007:2099
            
            out_filename <- glue("Beetle_GDD_{model}_{scenario}_2007-2099.nc")
            
            writeCDF(accumgdd2,
                        filename = file.path(out_dir, out_filename),
                        varname = "Beetle_GDD",
                        longname = "Beetle year growing degree days above 5.5C",
                        unit = "deg_K")
            
        }
    },
    mc.cores = cores
)
