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

cores <- 34

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

start_month <- 7

cl <- parallel::makeCluster(cores)
registerDoParallel(cl)

foreach(
    model = iter(models),
    scenario = iter(scenarios),
    .export = c("start_month", "tbase_k", "data_dir", "out_dir", "beetle_year", "get_gdd"),
    .packages = c("terra", "lubridate", "glue")
) %dopar% {
    in_filename <- glue("tavg_{model}_{scenario}_2006-2099_daily_gye.nc")
    r <- terra::rast(file.path(data_dir, in_filename))
    
    gdd <- terra::app(r, fun = get_gdd, tbase = tbase_k)

    time(gdd) <- time(r)

    accumgdd <- terra::tapp(gdd, index = beetle_year, fun = "sum")

    accumgdd2 <- subset(accumgdd, seq(2, nlyr(accumgdd) - 1))

    out_filename <- glue("Beetle_GDD_{model}_{scenario}_2007-2099.nc")
    writeRaster(accumgdd2, file.path(out_dir, out_filename))
}

stopCluster(cl)
