SDM Benchmark Study Part 1: Data Preparation
================

## Overview

This is just one of several posts that I will be making regarding a
study that I am working on, focused on benchmarking and comparing
different species distribution models (SDMs) using presence-only data. A
review of the literature that I wrote on this subject can be found
[here](https://raw.githubusercontent.com/benton-tripp/benton-tripp.github.io/main/_docs/lit_review.pdf).
My GitHub repo for this project can be accessed
[here](https://github.com/benton-tripp/presence-only-sdm).

This post is focused on gathering and preparing the data that will be
used throughout my project. More details on the data can be found in the
links in the next section.

## Data Sources

- State Boundary Data
  - [ArcGIS
    Hub](https://hub.arcgis.com/datasets/1612d351695b467eba75fdf82c10884f/explore?filters=eyJTVEFURV9BQkJSIjpbIkNPIiwiVlQiLCJOQyIsIk9SIl19&location=48.814319%2C163.610769%2C2.35)
    (Shapefile download with filters set to the 4 states in question)
- eBird Observation Data
  - [eBird Data Access](https://ebird.org/data/download)
- Raster Data
  - [DEM \| Source
    Page](https://www.sciencebase.gov/catalog/item/5540e111e4b0a658d79395d9)
    - Download by regions (Great Plains, Northeast, Northwest,
      Southeast, Southwest, and Upper Midwest)
  - [Urban Imperviousness \| Source Page](https://www.mrlc.gov/data)
    - *NLCD 2019 Percent Developed Imperviousness (CONUS), NLCD 2019
      Developed Imperviousness Descriptor (CONUS)*
  - [Land Cover](https://www.mrlc.gov/data)
    - *NLCD 2019 Land Cover (CONUS)*
  - [Canopy](https://www.mrlc.gov/data)
    - *NLCD 2016 USFS Tree Canopy Cover (CONUS)*
  - [Weather (min/max temperature, avg
    precipitation)](https://www.nacse.org/prism/)
    - Download weather raster data for “ppt”, “tmax”, “tmin”, 2017-2019
      at a 4km resolution and 30-year monthly normals at an 800m
      resolution
    - URL to download 4km data is:
      *<https://services.nacse.org/prism/data/public/4km/>*<VARIABLE>/<YEAR>
    - URL to download 800m data is
      *<https://services.nacse.org/prism/data/public/normals/800m/>*<VARIABLE>/<MONTH>
  - [Hydrography (Water Bodies & Coast) \| Source
    Page](https://apps.nationalmap.gov/downloader/#/)
    - [Download](https://prd-tnm.s3.amazonaws.com/StagedProducts/Small-scale/data/Hydrography/hydrusm010g.gdb_nt00897.tar.gz)
  - Vegetation Index
    - [USGS Earth Explorer](https://earthexplorer.usgs.gov/) (eVIIRS
      NDVI, 02/23/21-03/08/21 1km; 05/04/21-05/17/21 1km;
      09/07/21-09/20/21 1km; 11/30/21-12/13/21 1km)

## Basic Project Setup

``` r
# Load necessary libraries
library(sf)
library(terra)
library(stars)
library(auk)
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(purrr)
library(stringr)
library(fs)

# Hard coded path for downloaded data (in this case, saved in external hard drive)
ext.data.path <- "D:/AvianAnalyticsData"
```

## Boundary Data

There are four regions explored in this analysis from 2016 through 2019
(corresponding to the 4 boundary datasets):

- Colorado, 2016-2019
- North Carolina, 2016-2019
- Oregon, 2016-2019
- Vermont, 2016-2019

For simplicity, only data for the four states that are being used as
observation area
(<https://hub.arcgis.com/datasets/1612d351695b467eba75fdf82c10884f/explore?filters=eyJTVEFURV9BQkJSIjpbIkNPIiwiVlQiLCJOQyIsIk9SIl19&location=48.814319%2C163.610769%2C2.35>)
needs to be downloaded (as a .shp file). By default, the zipped files
should be saved in *US_State_Boundaries.zip*. Extract the files into a
folder of the same name within the data directory, and use the following
code to split the states into separate files:

``` r
# Define the state abbreviations
states <- c("CO", "NC", "OR", "VT")

# Set the output directory
state.output.dir <- "data/US_State_Boundaries"

output.paths <- file.path(state.output.dir, 
                         paste0(states, "_State_Boundaries.shp"))
if (!all(file.exists(output.paths))) {
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(state.output.dir)) {
    dir.create(state.output.dir)
  }
  
  # Set the path for the shapefile
  shapefile.path <- file.path(ext.data.path,
                              "US_State_Boundaries/US_State_Boundaries.shp")
  
  # Load the shapefile into an sf object
  gdf <- st_read(shapefile.path)
  
  # Transform the coordinate reference system to EPSG:5070
  gdf <- st_transform(gdf, crs = 5070)
  
  
  # Loop through state abbreviations and save individual shapefiles
  for (sa in states) {
    tryCatch({
      state.gdf <- gdf[gdf$STATE_ABBR == sa, ]
      output.path <- file.path(state.output.dir, paste0(sa, "_State_Boundaries.shp"))
      if (!file.exists(output.path)) st_write(state.gdf, output.path, quiet = T)
    }, error = function(e) {
      message(paste0('Failed to save shapefile for state ', sa, ': ', e$message))
    })
  }
}
```

## eBird Data

There are four regions explored in this analysis from 2016 through 2019
(corresponding to the 4 boundary datasets):

- Colorado, 2016-2019
- North Carolina, 2016-2019
- Oregon, 2016-2019
- Vermont, 2016-2019

Species included are:

- Belted Kingfisher
- Cedar Waxwing
- Downy Woodpecker
- Ruddy Duck
- Sanderling
- Sandhill Crane
- Sharp-shinned Hawk
- Wild Turkey

## Observation Data Pre-Processing

1.  Download each of the regions from the [eBird Download
    Page](https://ebird.org/data/download), filtered by date range (you
    will need to request access annually). They will initially be
    downloaded as compressed folders, so the contents will need to be
    extracted. There is also an eBird API, but the use of the API is
    beyond the scope of this document.
2.  Using the R [`auk`
    package](https://cornelllabofornithology.github.io/auk/), the
    contents can be filtered and saved.
3.  Do some basic pre-processing (i.e., select/rename relevant fields,
    filter by date, filter by approval, and filter by invalid
    observations).
4.  Save results; View data summary.

``` r
# Specify where your eBird datasets were downloaded to;
ebird.download.dirs <- list.dirs(file.path(ext.data.path, "ebird_downloads"))[-1]

# Specify where the outputs should be saved
ebird.output.dir <- file.path(ext.data.path, "ebird")

if (!dir.exists(ebird.output.dir)) {
  dir.create(ebird.output.dir)
}

# Define species
species <- c("Sandhill Crane", "Sharp-shinned Hawk",
             "Wild Turkey", "Downy Woodpecker",
             "Sanderling", "Cedar Waxwing", 
             "Belted Kingfisher", "Ruddy Duck")

# Parse eBird downloads
for (dir in ebird.download.dirs) {
  # Set the path to EBD text files
  # auk_set_ebd_path(dir, overwrite = T)
  cat("Processing species observation data at", dir, "\n")
  # List .txt files that start with "ebd_"
  in.files <- list.files(path = dir, pattern = "^ebd_.*\\.txt$", full.names = T)

  .sampl <- in.files[grepl("*sampling\\.txt$", in.files)]
  in.file <- in.files[!grepl("*sampling\\.txt$", in.files)]
  
  # Use regular expression to extract state abbreviation
  state.abbreviation <- sub(".*_US-([A-Z]{2})_.*", "\\1", in.file)
  out.file <- file.path(ebird.output.dir, paste0(state.abbreviation, ".txt"))
  
  out.sampling.file <- file.path(ebird.output.dir, 
                                 paste0("sampling_", state.abbreviation, ".txt"))
  if (!file.exists(out.file) & !file.exists(out.sampling.file)) {
    # Read in the filtered data using the `auk` library, saving to `out.file`
    auk_ebd(in.file, .sampl) %>%
      auk_species(species = species) %>% 
      auk_complete() %>% # Add this to keep only complete checklists
      auk_filter(file = out.file, file_sampling = out.sampling.file, 
                 overwrite=T, execute=T)
    
    # Remove "Problem" records
    df <- readr::read_delim(out.file, delim = "\t", 
                            show_col_types = F) %>% 
      suppressWarnings()
    df.samp <- readr::read_delim(out.sampling.file, delim = "\t", 
                                 show_col_types = F) %>%
      suppressWarnings()
    p <- problems(df)
    p.s <- problems(df.samp)
    if (nrow(p) > 0) {
      df <- df %>% slice(-p$row)
      readr::write_delim(df, out.file, delim="\t")
    }
    if (nrow(p.s) > 0) {
      df.samp <- df.samp %>% slice(-p$row)
      readr::write_delim(df.samp, out.sampling.file, delim="\t")
    }
  }
}

# Some basic pre-processing
preprocess.obs <- function(data.path) {
  data <- readr::read_delim(data.path, delim = "\t", show_col_types = F) %>%
    suppressMessages()
  names(data) <- gsub(" ", "\\.", tolower(names(data)))
  data <- data %>%
    filter(observation.date >= as.Date("2016-01-01") & 
             observation.date < as.Date("2020-01-01") &
             approved == 1 & observation.count != "X") %>%
    dplyr::select(common.name, observation.count, latitude, longitude) %>% #observation.date
    group_by(common.name, latitude, longitude) %>% 
    summarize(observation.count = sum(as.numeric(observation.count), na.rm=T),
              .groups="keep") %>%
    ungroup() %>%
    as.data.table()
  return(data)
}

ebird.path <- "data/ebird"
if (!dir.exists(ebird.path)) dir.create(ebird.path)

obs <- purrr::map(states, function(.x) {
  cat("Getting preprocessed observation data in", .x, "\n")
  out.file <- file.path(ebird.path, paste0(.x, ".csv"))
  if (!file.exists(out.file)) {
    out <- preprocess.obs(file.path(ebird.output.dir, paste0(.x, ".txt")))
    fwrite(out, out.file)
  } else {
    out <- fread(out.file)
  }
  out
})

names(obs) <- states
```

``` r
# View summary of data
map_df(states, ~obs[[.x]][, .(.N, state=.x), by=.(common.name)]) %>%
  setorderv(c("common.name", "state")) %>%
  rename(Species="common.name", `Total Obs.`="N", State="state") %>%
  knitr::kable( 
    caption = "Summary of Observation Data", 
    align = c("l", "l", "l"),
  )
```

| Species            | Total Obs. | State |
|:-------------------|:-----------|:------|
| Belted Kingfisher  | 4552       | CO    |
| Belted Kingfisher  | 4854       | NC    |
| Belted Kingfisher  | 6367       | OR    |
| Belted Kingfisher  | 2089       | VT    |
| Cedar Waxwing      | 3447       | CO    |
| Cedar Waxwing      | 4443       | NC    |
| Cedar Waxwing      | 8896       | OR    |
| Cedar Waxwing      | 4009       | VT    |
| Downy Woodpecker   | 7513       | CO    |
| Downy Woodpecker   | 11280      | NC    |
| Downy Woodpecker   | 8848       | OR    |
| Downy Woodpecker   | 4699       | VT    |
| Ruddy Duck         | 1728       | CO    |
| Ruddy Duck         | 1345       | NC    |
| Ruddy Duck         | 2097       | OR    |
| Ruddy Duck         | 51         | VT    |
| Sanderling         | 131        | CO    |
| Sanderling         | 1495       | NC    |
| Sanderling         | 665        | OR    |
| Sanderling         | 39         | VT    |
| Sandhill Crane     | 1534       | CO    |
| Sandhill Crane     | 125        | NC    |
| Sandhill Crane     | 2472       | OR    |
| Sandhill Crane     | 77         | VT    |
| Sharp-shinned Hawk | 2258       | CO    |
| Sharp-shinned Hawk | 1408       | NC    |
| Sharp-shinned Hawk | 2807       | OR    |
| Sharp-shinned Hawk | 758        | VT    |
| Wild Turkey        | 2617       | CO    |
| Wild Turkey        | 2420       | NC    |
| Wild Turkey        | 2460       | OR    |
| Wild Turkey        | 2211       | VT    |

Summary of Observation Data

Below is a brief snapshot of the data, filtered to
`common.name == "Sandhill Crane"` in Oregon:

``` r
# OR Sandhill Crane Map
or.shc <- obs$OR[common.name == "Sandhill Crane"]

ggplot(data = ggplot2::map_data("state") %>% filter(region == "oregon")) +
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill = "#ffffff", color = "black") +
    geom_point(data = or.shc, 
               aes(x = longitude, y = latitude), 
               size=2, alpha=1, fill = "red", shape=21) +
    coord_map() +
    labs(title = "Sandhill Crane observations in OR") + 
    theme_minimal() +
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
```

![](https://raw.githubusercontent.com/benton-tripp/benton-tripp.github.io/main/_posts/2023-09-10-sdm-benchmark-study-part-1-data-preparation_files/figure-gfm/nc-obs-example-1.png)

## Rasters

Raster data is a type of digital image represented by a grid of pixels,
where each pixel has an associated value that represents information
about a particular geographic area. It is one of the two primary ways to
represent geospatial data (the other being vector data). In
environmental studies and geographic analysis, raster data is
extensively used to represent continuous data such as elevation,
temperature, or land cover. In this project, raster data serves as the
basis for various explanatory variables that influence bird species
distributions, including digital elevation, urban imperviousness, land
cover, canopy, weather, hydrography, and vegetation index. Each of these
variables represents environmental conditions that can affect bird
habitat and distribution.

## Raster Pre-Processing

Before using raster data for analysis, it is crucial to preprocess them
to ensure compatibility, accuracy, and relevancy. Raster pre-processing
includes a series of steps that prepare and standardize the raster
datasets for subsequent analyses. Along with some general pre-processing
steps (and corresponding code), the following sections outline specific
pre-processing steps tailored for each type of explanatory variable
being used in the analysis.

### General Raster Pre-Processing Steps

1.  Load the environmental raster datasets & State Boundary data
2.  Reproject to CRS EPSG:5070
3.  Resample each to a common resolution (e.g., 5000 x 5000 meters)
4.  Mask each raster using State Boundary data

Following the pre-processing of all rasters, it will also be necessary
to ensure that the extents of each of the rasters (by state) match, and
to crop them if necessary to produce uniform extents.

``` r
# `terra` equivalent of base R `any()` function
terra.any <- function(r) {
  freqs <- freq(r)
  any(freqs[,1] == 1)
}

# Function to fix NA values in a raster (e.g., NULL values are equal to 999, 
# but should be NA)
basic.na.fix <- function(r, na.val) {
  msk <- r == na.val
  cat("Checking for improperly formatted NA values...\n")
  if (terra.any(msk)) {
    r[msk] <- NA
  }
  r
}


# Recursive function to list directories containing .adf files
# Usage:
# raster.directories <- list.raster.dirs.recursive("your_start_directory_path")
list.adf.rasters.recursive <- function(directory.path) {
  
  # List all directories recursively
  all.dirs <- list.dirs(directory.path, recursive = T, full.names = T)
  
  # Filter directories that contain .adf files
  raster.dirs <- all.dirs[sapply(all.dirs, function(dir.path) {
    any(grepl("\\.adf$", list.files(dir.path, full.names = F)))
  })]
  
  return(raster.dirs)
}


general.raster.preprocessing <- function(
    raster.name,
    out.raster.name,
    data.path = "data",
    state.boundary.path = "data/US_State_Boundaries",
    states = c("NC", "CO", "OR", "VT"),
    state.file.suffix = "_State_Boundaries.shp",
    out.path = "../gis630/data",
    crs = 5070,
    wildcard = "*.tif",
    resolution = 5000,
    agg="bilinear",
    recursive.adf.path=F,
    crop.by.state=T,
    na.val=NULL
) {
  # Example usage:
  # general.raster.preprocessing(
  #   raster.name = "canopy/nlcd_2016_treecanopy_2019_08_31",
  #   out.raster.name = "canopy",
  #   out.path = "data/canopy"
  # )
  
  # Make sure output directory exists
  if (!dir.exists(out.path)) dir.create(out.path)
  
  if (!recursive.adf.path) {
    # List all the raster files in the workspace directory
    rasters <- list.files(path = file.path(data.path, raster.name), 
                          pattern = wildcard, 
                          full.names = T)
  } else {
    # List directories containing .adf files
    rasters <- list.adf.rasters.recursive(data.path)
  }
  
  cat("Rasters:", paste(rasters, collapse=", "), "\n")

  # Loop through each raster file in the list
  for (i in 1:length(rasters)) {
    
    raster.path <- rasters[[i]]
    
    if (crop.by.state) {
      
      # Read raster
      cat(paste0("[", i, "/", length(rasters), "]:"), "Reading", raster.path, "\n")
      raster <- rast(raster.path)
      
      # Loop through each state
      for (state in states) {
        
        state.out.file <- file.path(out.path, 
                                    paste0(out.raster.name, "_", 
                                           state, ".tif"))
        if (!file.exists(state.out.file)) {
          # Construct the file path to the state boundary shapefile
          state.path <- file.path(state.boundary.path, 
                                  paste0(state, state.file.suffix))
          
          # Load state shape and reproject it to match raster's CRS
          state.shape <- vect(state.path)
          reprojected.shape <- project(state.shape, crs(raster))
          
          # Crop raster by reprojected state shape
          cat("\tCropping raster for", state, "\n")
          # masked.raster <- terra::mask(raster, reprojected.shape)
          masked.raster <- terra::crop(raster, reprojected.shape, mask=T)
          
          # Reproject
          cat("\tReprojecting raster for", state, "\n")
          reprojected.raster <- project(masked.raster, 
                                        crs(paste0("EPSG:", crs)))
          
          # Clean up
          rm(masked.raster)
          gc()
          
          # Fixing NA values
          if (!is.null(na.val)) {
            reprojected.raster <- basic.na.fix(reprojected.raster, na.val)
          }
          
          template.raster <- ext(reprojected.raster) %>% 
            rast(res=rep(resolution, 2), crs=crs(reprojected.raster))
          
          
          
          # Resample
          current.res <- terra::res(reprojected.raster)
          cat("\tCurrent Resolution:", current.res, "\n")
          cat("\tTarget Factor:", resolution/terra::res(reprojected.raster)[1], "\n")
          if (any(current.res != resolution)) {
            cat("\tResampling raster for", state, "\n")
            # resampled.raster <- terra::aggregate(
            #   reprojected.raster, 
            #   fact=c(resolution/current.res[1], resolution/current.res[2]), 
            #   fun = agg, 
            #   expand = T) %>% suppressWarnings()
            resampled.raster <- terra::resample(reprojected.raster, 
                                                template.raster,
                                                method=agg)
          } else {
            resampled.raster <- reprojected.raster
          }
          
          # Clean up
          rm(reprojected.raster)
          gc()
          
          # Rename final raster
          names(resampled.raster) <- paste0(out.raster.name, "_", state)
          
          # Save raster
          cat("\tSaving to", state.out.file, "\n")
          writeRaster(resampled.raster, state.out.file)
          
          # Clean up
          rm(resampled.raster)
          gc()
        }
      }
    } else {
      path.parts <- stringr::str_split(raster.path, "\\.")[[1]]
      if (length(path.parts) > 1) path.parts <- path.parts[[1]]
      path.parts <- stringr::str_split(path.parts, "/")[[1]]
      r.name <- path.parts[[length(path.parts)]]
      out.file <- file.path(out.path, 
                            paste0(out.raster.name, "_", 
                                   r.name, ".tif"))
      if (!file.exists(out.file)) {
        
        # Read raster
        cat(paste0("[", i, "/", length(rasters), "]:"), "Reading", raster.path, "\n")
        raster <- rast(raster.path)
        
        # Reproject
        cat("\tReprojecting raster...\n")
        reprojected.raster <- project(raster, 
                                      crs(paste0("EPSG:", crs)))
        
        # Clean up
        rm(raster)
        gc()
        
        # Resample
        current.res <- terra::res(reprojected.raster)
        cat("\tCurrent Resolution:", current.res, "\n")
        cat("\tTarget Factor:", resolution/terra::res(reprojected.raster)[1], "\n")
        if (any(current.res != resolution)) {
          cat("\tResampling raster...\n")
          resampled.raster <- terra::aggregate(
            reprojected.raster, 
            fact=c(resolution/current.res[1], resolution/current.res[2]), 
            fun = agg, 
            expand = T) %>% suppressWarnings()
        } else {
          resampled.raster <- reprojected.raster
        }
        
        # Clean up
        rm(reprojected.raster)
        gc()
        
        # Save raster
        cat("\tSaving to", out.file, "\n")
        writeRaster(resampled.raster, out.file)
        
        # Clean up
        rm(resampled.raster)
        gc()
      }
    }
    
    # Clean up
    rm(raster)
    gc()
  }
  cat("Finished general pre-processing.\n")
}
```

### Digital Elevation Model (DEM) Pre-Processing

The following steps are taken when merging and Pre-processing the DEM
data downloaded from its
[source](https://www.sciencebase.gov/catalog/item/5540e111e4b0a658d79395d9):

1.  Merge multiple DEM parts (the raw data download is divided into
    Great Plains, Northeast, Northwest, Southeast, Southwest, and Upper
    Midwest regions).
2.  Apply general pre-processing of DEM raster (see section on general
    pre-processing).

``` r
# Combine DEM rasters (Source data has the US split into parts)

# By default, `terra` does most of its work on the disk (as 
# opposed to RAM). However, the temp data for this particular
# process was using ~70GB of storage. This function is set up 
# to reduce this overhead as much as possible.
combine.dem.parts <- function(dem.dir = "data/dem/raw_dem", 
                              output.dir = "data/dem", 
                              output.raster.name = "mosaic_dem") {
  
  if (!dir.exists(output.dir)) dir.create(output.dir)
  
  out.name <- file.path(output.dir, paste0(output.raster.name, ".tif"))
  
  # Get list of rasters recursively
  rasters <- list.adf.rasters.recursive(dem.dir)
  
  # Read the first raster
  cat("1: Reading first raster from", rasters[1], "out of", 
      length(rasters), "rasters. \nWriting to", out.name, ".\n\n")
  
  # Save the combined raster
  writeRaster(rast(rasters[1]), 
              out.name,
              gdal=c("COMPRESS=DEFLATE", "TFW=YES"),
              verbose=T,
              overwrite=T)
  cat("----------------------------\n")
  
  # If there are more rasters, mosaic them iteratively
  if (length(rasters) > 1) {
    for (i in 2:length(rasters)) {
      cat(paste0(i, ":"), "Reading raster from", rasters[i], "out of", 
          length(rasters), "rasters.\nAdding [", i, "/", length(rasters), 
          "] rasters to mosaic.\n\n")
      
      # Save the combined raster
      writeRaster(mosaic(rast(out.name), rast(rasters[i]), fun=mean), 
                  out.name,
                  gdal=c("COMPRESS=DEFLATE", "TFW=YES"),
                  verbose=T,
                  overwrite=T)
      cat("----------------------------\n")
      
      # Clean up
      gc()
    }
  }
  
  cat("Mosaic completed successfully.\n")
}

dem.path <- file.path(ext.data.path, "dem")

if (!file.exists(file.path(dem.path, "mosaic_dem.tif"))) {
  # Combine rasters 
  combine.dem.parts(dem.dir=file.path(dem.path, "raw_dem"),
                    output.dir=dem.path)
}

if (!all(file.exists(paste0("data/dem/dem_", states, ".tif")))) {
  # General pre-processing on each of the DEM parts
  general.raster.preprocessing(
    data.path = ext.data.path,
    raster.name = "dem", 
    out.raster.name = "dem",
    out.path = "data/dem",
    resolution = 5000,
    wildcard="\\.tif$",
    agg="bilinear"
  )
}
```

### Weather Rasters

#### Download Weather Raster Data

In this section, the weather raster data is fetched from its source
programmatically. The weather data variables considered are
precipitation (ppt), maximum temperature (tmax), and minimum temperature
(tmin). The annual weather data at a 4km resolution for the years 2017,
2018, and 2019. Additionally, monthly 30-year normal data is retrieved
at an 800m resolution. After downloading, the data is stored in zipped
format and then extracted.

``` r
# Downloads and processes weather raster data for specified variables and 
# years at a 4km resolution and 30-year monthly normals at an 800m resolution.
get.weather.data <- function(data.path,
                             out.dir="weather") {
  
  out.path <- file.path(data.path, out.dir)
  if (!dir.exists(out.path)) dir.create(out.path)
  
  vars <- c("ppt", "tmax", "tmin")
  
  # Setup for yearly 4km resolution 
  yrs <- c(2017, 2018, 2019)
  pairs <- expand.grid(vars, yrs)
  
  # Setup for 30 year normal monthly 800m resolution
  mnths <- sprintf("%02d", 1:12)
  norm.pairs <- expand.grid(vars, mnths)
  
  ### Get Raster Data ######
  cat("Getting explanatory Weather Rasters...\n")
  
  # Data documentation:
  # https://www.prism.oregonstate.edu/documents  /PRISM_downloads_web_service.pdf
  
  # 4km yearly data (for 2017-2019)
  for (i in 1:nrow(pairs)) {
    v <- pairs[i, 1]
    y <- pairs[i, 2]
    
    dwnld.out <- file.path(out.path, paste0(v, "_", y, ".zip"))
    dwnld.path <- file.path(out.path, paste0(v, "_", y))
    
    if (!dir.exists(dwnld.path)) {
      dir.create(dwnld.path)
      url <- paste0("https://services.nacse.org/prism/data/public/4km/", v, "/", y)
      cat(paste("Downloading weather data from", url, "...\n"))
      download.file(url, dwnld.out)
      cat(paste("Saved", v, "/", y, "to", dwnld.out, "\n"))
      unzip(dwnld.out, exdir = dwnld.path)
      cat(paste("Extracted", v, "/", y, "from", dwnld.out, "to", dwnld.path, "\n"))
      # file.remove(dwnld.out)
    }
  }
  
  # 800m monthly data (30 year normals)
  for (i in 1:nrow(norm.pairs)) {
    v <- norm.pairs[i, 1]
    m <- norm.pairs[i, 2]
    
    dwnld.out <- file.path(out.path, paste0(v, "_", m, ".zip"))
    dwnld.path <- file.path(out.path, paste0(v, "_", m))
    
    if (!dir.exists(dwnld.path)) {
      url <- paste0("https://services.nacse.org/prism/data/public/normals/800m/", 
                    v, "/", m)
      cat(paste("Downloading weather data from", url, "...\n"))
      download.file(url, dwnld.out)
      cat(paste("Saved", v, "/", m, "to", dwnld.out, "\n"))
      unzip(dwnld.out, exdir = dwnld.path)
      cat(paste("Extracted", v, "/", m, "from", dwnld.out, "to", dwnld.path, "\n"))
      # file.remove(dwnld.out)
    }
  }
  
  cat("Finished getting weather data.\n")
}

get.weather.data(ext.data.path)
```

#### Aggregate Weather Raster Data

The downloaded weather rasters are aggregated, considering both the
yearly data and monthly 30-year normals. The aggregation involves
determining weights for the different resolutions (4km and 800m) and
then combining them accordingly.

``` r
aggregate.weather <- function(data.path,
                              raster.dir="weather",
                              out.dir="weather/aggregated",
                              states=c("CO", "NC", "OR", "VT"),
                              vars=c("ppt", "tmax", "tmin"),
                              yrs=2017:2019) {
  out.path <- file.path(data.path, out.dir)
  raster.path <- file.path(data.path, raster.dir)
  if (!dir.exists(out.path)) dir.create(out.path)
  
  max.temp.data <- "max_temp.tif"
  min.temp.data <- "min_temp.tif"
  avg.prcp.data <- "avg_prcp.tif"
  
  pairs <- expand.grid(vars, yrs)
  mnths <- sprintf("%02d", 1:12)
  norm.pairs <- expand.grid(vars, mnths)
  
  for(v in vars) {
    raster.name <- ifelse(v == "tmax", max.temp.data, 
                          ifelse(v == "tmin", min.temp.data, 
                                 avg.prcp.data))
    v.agg.file <-  file.path(out.path, paste0("aggregated_", raster.name))
    if (!file.exists(v.agg.file)) {
      
      cat(paste0("Checking ", v, "...\n"))
      v.path <- file.path(out.path, v)
      agg.func <- switch(v,
                         tmax = max,
                         tmin = min,
                         ppt = mean,
                         function(x) x) # default to just return the value
      
      # Yearly rasters
      rasters <- list.files(raster.path, full.names=T, 
                            recursive=T, pattern="ppt_stable_4km.*\\.bil$")
      r <- rast(rasters)
      agg.r <- terra::app(r, fun=agg.func)
      
      # Monthly rasters
      norm.rasters <- list.files(raster.path, full.names=T, 
                                 recursive=T, pattern="ppt_30yr.*\\.bil$")
      
      r.norm <- rast(norm.rasters)
      agg.r.norm <- terra::app(r.norm, fun=agg.func)
      
      # Aggregate with weights
      initial.weight.4km <- 3.0
      initial.weight.800m <- 3.0 / 30.0
      total.weight <- initial.weight.4km + initial.weight.800m
      normalized.weight.4km <- initial.weight.4km / total.weight
      normalized.weight.800m <- initial.weight.800m / total.weight
      
      # Resample to match res
      agg.r.resampled <- resample(agg.r, agg.r.norm, method="bilinear")

      combined.raster <- (agg.r.resampled * normalized.weight.4km) + 
        (agg.r.norm * normalized.weight.800m)
      
      # Write result
      writeRaster(combined.raster, 
                  file = v.agg.file, 
                  overwrite = T)
    }
  }
  
  cat("Finished weather data pre-processing.\n")
}

aggregate.weather(ext.data.path)
```

#### Average Precipitation, Min/Max Temperature Pre-Processing

``` r
# Average Precipitation

# Apply General Raster Pre-Processing to precipitation
if (!all(file.exists(paste0("data/prcp/avg_prcp_", states, ".tif")))) {
  general.raster.preprocessing(
    data.path=ext.data.path,
    raster.name="weather/aggregated", 
    out.path="data/prcp",
    out.raster.name="avg_prcp",
    resolution=5000,
    wildcard="avg_prcp\\.tif$",
    agg="bilinear"
  )
}

# Minimum Temperature 

# Apply General Raster Pre-Processing to minimum temp
if (!all(file.exists(paste0("data/tmin/tmin_", states, ".tif")))) {
  general.raster.preprocessing(
    data.path=ext.data.path,
    raster.name="weather/aggregated", 
    out.path="data/tmin",
    out.raster.name="tmin",
    resolution=5000,
    wildcard="min_temp\\.tif$",
    agg="bilinear"
  )
}

# Maximum Temperature 

# Apply General Raster Pre-Processing to maximum temp
if (!all(file.exists(paste0("data/tmax/tmax_", states, ".tif")))) {
  general.raster.preprocessing(
    data.path=ext.data.path,
    raster.name="weather/aggregated", 
    out.path="data/tmax",
    out.raster.name="tmax",
    resolution=5000,
    wildcard="max_temp\\.tif$",
    agg="bilinear"
  )
}
```

### Hydrography Raster Data Conversion

#### Extract from ESRI Geodatabase

This section handles the extraction of hydrographic data from an ESRI
Geodatabase. Two of the data sets, Coastline and Waterbody, within the
geodatabase are converted into individual shapefile formats.

``` r
convert.hydro.gdb <- function(data.path = "data/hydrography/",
                              gdb.name = "hydrusm010g.gdb_nt00897/hydrusm010g.gdb",
                              data.sets = c("Coastline", "Waterbody"),
                              output.path = "data/hydrography/") {
  if (!dir.exists(output.path)) dir.create(output.path)
  gdb.path <- file.path(data.path, gdb.name)
  # Iterate over each dataset and convert
  for (ds in data.sets) {
    out.file <- file.path(output.path, paste0(ds, ".shp"))
    if (!file.exists(out.file)) {
      # Read the geodatabase layer using st_read
      data <- st_read(dsn = file.path(data.path, gdb.name), layer = ds, quiet = T)
      # Write the shapefile
      st_write(obj = data, dsn = out.file, quiet = T)
      cat(sprintf("Converted %s from geodatabase to shapefile.\n", ds))
    }
  } 
}

convert.hydro.gdb(data.path=file.path(ext.data.path, "hydrography"),
                  output.path=file.path(ext.data.path, "hydrography")) %>%
  suppressWarnings()
```

#### Convert Waterbody Shapefile to Raster & Pre-Process

The waterbody data previously extracted as a shapefile is then
transformed into raster format. The steps entail:

1.  Cleaning the data by filtering out specific undesirable features.
2.  Reprojecting the data into a standard coordinate reference system
    (CRS).
3.  Rasterizing the shapefile data with a default resolution.
4.  Using dilation to expand the waterbody boundaries in the raster.
    Dilation is a process to enlarge or expand a feature. In this
    context, the waterbody raster is dilated in multiple stages to
    represent different proximities to the water.
5.  Applying the general raster pre-processing steps.

``` r
update.wb.array <- function(r) {
  r.updated <- copy(r)
  
  w <- matrix(1,3,3)
  # Perform the dilation operation multiple times
  dilated1 <- terra::focal(r, w = w, fun = max)
  dilated2 <- terra::focal(dilated1, w =w, fun = max)
  dilated3 <- terra::focal(dilated2, w = w, fun = max)
  dilated4 <- terra::focal(dilated3, w = w, fun = max)
  
  r.updated[dilated4 == 1] <- 0.2
  r.updated[dilated3 == 1] <- 0.4
  r.updated[dilated2 == 1] <- 0.6
  r.updated[dilated1 == 1] <- 0.8
  r.updated[r == 1] <- 1
  
  return(r.updated)
}

# Convert waterbody shapefile to raster
process.waterbody.shp <- function(out.path,
                                  data.path="data/hydrography",
                                  states=c("CO", "VT", "NC", "OR")) {
  out.file <- file.path(out.path, 'Waterbody.tif')
  if (!file.exists(out.file)) {
    if (!dir.exists(out.path)) dir.create(out.path)
    
    cat("Cleaning waterbody data...\n")
    gdf <- st_read(file.path(data.path, 'Waterbody.shp'))
    
    # Filter out records where Feature == "Lake Dry"
    gdf <- gdf %>% filter(Feature != "Lake Dry")
    
    gdf$waterbody <- 1
    
    cat("Reprojecting...\n")
    gdf <- st_transform(gdf, crs = 5070)
    
    cat("Converting to raster using template...\n")
    
    # Initialize raster resolution at 1km prior to dilation
    # (will be updated to desired resolution)
    template.raster <- ext(gdf) %>% 
      rast(res=rep(1e3, 2), crs=crs(gdf))
    
    r <- terra::rasterize(gdf, template.raster, 
                               field="waterbody", values=1, background=0)
    
    # Add dilation
    r.updated <- update.wb.array(r)
    
    # Save raster
    terra::writeRaster(r.updated, out.file, overwrite=T)
  }
}

process.waterbody.shp(out.path=file.path(ext.data.path, "waterbody"),
                      data.path=file.path(ext.data.path, "hydrography"))

# Apply General Raster Pre-Processing to Waterbodies
if (!all(file.exists(paste0("data/waterbody/waterbody_", states, ".tif")))) {
  general.raster.preprocessing(
    data.path=ext.data.path,
    raster.name="waterbody", 
    out.path="data/waterbody",
    out.raster.name="waterbody",
    resolution=5000,
    wildcard="\\.tif$",
    agg="bilinear"
  )
}
```

#### Convert Coastline Shapefile to Raster & Pre-Process

The coastline data, like the waterbody data, undergoes a conversion from
shapefile to raster format. The steps are:

1.  Basic data cleaning and reprojection to a standard CRS.
2.  A buffer is applied to the coastline vectors. This adds a specified
    distance around the coastline, effectively creating a zone around
    the coast.
3.  Rasterization of the buffered coastline data.
4.  Dilation of the coastline raster to indicate proximity zones.
5.  General raster pre-processing.

``` r
# Convert coastline shapefile to raster
process.coastline.shp <- function(out.path,
                                  data.path="data/hydrography",
                                  states=c("CO", "VT", "NC", "OR")) {
  out.file <- file.path(out.path, 'Coastline.tif')
  if (!file.exists(out.file)) {
    if (!dir.exists(out.path)) dir.create(out.path)
    
    cat("Cleaning coastline data...\n")
    
    # Read data
    gdf <- st_read(file.path(data.path, 'Coastline.shp'))
    
    # Reproject
    cat("Reprojecting...\n")
    gdf <- st_transform(gdf, crs = 5070)
    
    # Buffer
    cat("Buffering coastline vectors...\n")
    gdf <- st_buffer(gdf, dist = 500)
    
    gdf$coastline <- 1
    
    cat("Converting to raster using template...\n")
    
    # Initialize raster resolution at 1km prior to dilation
    # (will be updated to desired resolution)
    template.raster <- ext(gdf) %>% 
      rast(res=rep(1e3, 2), crs=crs(gdf))
    
    r <- terra::rasterize(gdf, template.raster, 
                          field="coastline", values=1, background=0)
    
    # Add dilation
    r.updated <- update.wb.array(r)
    
    # Save raster
    terra::writeRaster(r.updated, out.file, overwrite=T)
  }
}


process.coastline.shp(out.path=file.path(ext.data.path, "coastline"),
                      data.path=file.path(ext.data.path, "hydrography"))

# Apply General Raster Pre-Processing to coastline
if (!all(file.exists(paste0("data/coastline/coastline_", states, ".tif")))) {
  general.raster.preprocessing(
    data.path=ext.data.path,
    raster.name="coastline", 
    out.path="data/coastline",
    out.raster.name="coastline",
    resolution=5000,
    wildcard="\\.tif$",
    agg="bilinear"
  )
}
```

### General Pre-Processing of Remaining Rasters

The remaining rasters, upon downloading them from their respective
sources, require no additional pre-processing other than the “general”
pre-processing steps described previously:

- Urban Imperviousness
- Land Cover
- NDVI
- Tree Canopy

For land cover, note the 20 different Land Cover hierarchical categories
(i.e., each of them falls under a “parent” category, see [National Land
Cover Database Class Legend and
Description](https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description)):

- Water:
  - 11: Open Water
  - 12: Perennial Ice/Snow
- Developed
  - 21: Developed, Open Space
  - 22: Developed, Low Intensity
  - 23: Developed, Medium Intensity
  - 24: Developed, High Intensity
- Barren
  - 31: Barren Land (Rock/Sand/Clay)
- Forest
  - 41: Deciduous Forest
  - 42: Evergreen Forest
  - 43: Mixed Forest
- Shrubland
  - 51: Dwarf Shrub
  - 52: Shurb/Scrub
- Herbaceous
  - 71: Grassland/Herbaceous
  - 72: Sedge/Herbaceous
  - 73: Lichens
  - 74: Moss
- Planted/Cultivated
  - 81: Pasture/Hay
  - 82: Cultivated Crops
- Wetlands
  - 90: Woody Wetlands
  - 95: Emergent Herbaceous Wetlands

``` r
# Urban Imperviousness Pre-Processing

if (!all(file.exists(paste0("data/urban_imperviousness/urban_imperviousness_", 
                            states, ".tif")))) {
  # Use the combined raster for the general raster pre-processing
  general.raster.preprocessing(
    data.path = ext.data.path,
    raster.name="urban_imperviousness", 
    out.raster.name="urban_imperviousness",
    out.path="data/urban_imperviousness",
    resolution=5000,
    wildcard="\\.tif$",
    agg="bilinear"
  )
}


# Land Cover Pre-Processing

if (!all(file.exists(paste0("data/land_cover/land_cover_", states, ".tif")))) {
  # Apply "general raster preprocessing" 
  general.raster.preprocessing(
    data.path = ext.data.path,
    raster.name="land_cover/nlcd_2019_land_cover_l48_20210604", 
    out.raster.name="land_cover",
    out.path="data/land_cover",
    resolution=5000,
    wildcard="\\.tif$",
    agg="near"
  )
}


# NDVI Pre-Processing

# Iterate through each season
for (season in c("Spring", "Summer", "Fall", "Winter")) {
  if (!all(file.exists(file.path("data/NDVI", paste0(season, "_NDVI_", 
                                                     states, ".tif"))))) {
    cat(paste0("Applying raster preprocessing for ", season, "...\n"))
    tryCatch({
      # Using the combined raster, apply "general raster preprocessing" 
      # (resample, reproject, mask)
      
      general.raster.preprocessing(
        data.path=ext.data.path,
        raster.name=paste0("NDVI/US_eVSH_NDVI-", season, "-2021"), 
        out.raster.name=paste0(season, "_NDVI"),
        out.path="data/NDVI",
        wildcard="*1KM\\.VI_NDVI.*\\.tif$",
        resolution=5000,
        agg="bilinear")
      cat("-----------------\n")
    }, error = function(e) {
      cat(paste0("An error occurred while processing ", season, ": ", e$message, "\n"))
    })
  }
}

# Canopy Pre-Processing

# Apply "general raster preprocessing" 

if (!all(file.exists(paste0("data/canopy/canopy_", states, ".tif")))) {
  general.raster.preprocessing(
    data.path=ext.data.path,
    raster.name="canopy/nlcd_tcc_CONUS_2016_v2021-4", 
    out.raster.name="canopy",
    out.path="data/canopy",
    resolution=5000,
    wildcard="\\.tif$",
    agg="bilinear"
  )
}
```

## Final Raster Pre-Processing Steps

### Ensuring Shape/Size Conformity

In order to use the rasters to be saved as layers in a single raster,
they must exactly conform in resolution, CRS, and shape/extent. Due to
slight variations in the datasets, the raster extents are not perfect
(although they are nearly so). Now that all of the rasters are
available, a common “intersect” region can be determined for each state,
and all of the rasters for each state can be cropped to fit that region.

Once this is completed, the rasters can be joined into combined,
multi-layer rasters and saved.

``` r
# Function to intersect two extents
intersect.extents <- function(ext1, ext2) {
  xmin <- max(ext1[1], ext2[1])
  xmax <- min(ext1[2], ext2[2])
  ymin <- max(ext1[3], ext2[3])
  ymax <- min(ext1[4], ext2[4])
  
  if (xmin > xmax | ymin > ymax) {
    stop("The extents do not overlap!")
  }
  
  ext(c(xmin, xmax, ymin, ymax))
}

# Get all rasters

all.rasters <- map(states, function(s) {
  files <- list.files("data", pattern=paste0("_", s, "\\.tif$"), 
                      recursive=T, full.names=T)
  
  rasters <- map(files, ~rast(.x))
  
  # Get extent intersect, and update raster extents
  extents <- map(rasters, ~ext(.x) %>% as.vector())
  ext.intersect <- reduce(extents, intersect.extents)
  rasters <- map(rasters, ~crop(.x, ext.intersect))
  
  # Align all rasters to common grid
  rasters <- map(rasters[2:length(rasters)], ~project(.x, rasters[[1]]))
  
  # Check extents, resolutions, and CRS of each raster
  extents <- map(rasters, ~ext(.x) %>% as.vector())
  resolutions <- map(rasters, ~res(.x))
  .crs <- map(rasters, ~crs(.x))
  
  # Set names
  .names <- map_chr(rasters, ~names(.x))
  names(files) <- .names
  names(rasters) <- .names
  names(extents) <- .names
  names(resolutions) <- .names
  names(.crs) <- .names
  # Return list
  list(
    files=files,
    rasters=rasters %>% reduce(c), # Make raster stack
    ext=extents,
    ext.intersect=ext.intersect,
    res=resolutions,
    crs=.crs,
    names=.names
  )
})

names(all.rasters) <- states
if (!all(file.exists(paste0("data/final_rasters/", states, ".tif")))) {
  if (!dir.exists("data/final_rasters")) dir.create("data/final_rasters")
  walk(states, ~writeRaster(all.rasters[[.x]]$rasters, 
                            paste0("data/final_rasters/", .x, ".tif"),
                            overwrite=T))
} 
```

## Combine Raster and Observation Data

When both the observation data and raster data pre-processing has been
completed, points can be extracted from each of the rasters
corresponding to the observation points. These tabular datasets are what
will be used to train/test presence-only models.

``` r
# Combine rasters and observations
for (state in states) {
  
  input.file <- file.path("data", "ebird", paste0(state, ".csv"))
  output.file <- file.path("data", "ebird", paste0(state, ".shp"))
  
  if (!file.exists(output.file)) {
    cat(sprintf("Reading %s data...\n", state))
    
    df <- fread(input.file)
    
    # Convert the bird data to a sf object
    cat(sprintf("Converting %s data to sf...\n", state))
    geo.df <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
    
    
    # Get centroid of the entire sightings
    cat(sprintf("Getting centroid of %s data...\n", state))
    centroid <- st_centroid(st_union(geo.df))
    
    # Assign each point a distance attribute, being the distance from the centroid
    cat(sprintf("Calculating distances from centroid for %s...\n", state))
    geo.df$distance <- st_distance(geo.df, centroid)
    
    # Sort the sf object by the distance attribute
    cat(sprintf("Sorting by distance from centroid in %s data...\n", state))
    geo.df <- geo.df %>% arrange(distance)
    
    # Drop the distance attribute 
    cat(sprintf("Saving %s data to %s...\n", state, output.file))
    geo.df$distance <- NULL
    
    st_write(geo.df, output.file) %>% suppressWarnings()
    cat("--------------\n")
  }
  
  if (!dir.exists(file.path("data", "final"))) dir.create(file.path("data", "final"))
  out.file.all <- file.path("data", "final", paste0("all_data_", state, ".rds"))
  if (!file.exists(out.file.all)) {
    r <- all.rasters[[state]]$rasters
    r.names <- all.rasters[[state]]$names
    
    cat(sprintf("Extracting points to values for %s...\n", state))
    # Load observations shapefile
    geo.df <- st_read(output.file) 
    
    
    # Define target CRS and update
    target.crs <- "EPSG:5070"
    cat(sprintf("Updating CRS for %s...\n", state))
    geo.df <- st_transform(geo.df, target.crs)
    
    # Extract raster values
    for (r.name in r.names) {
      cat("\tExtracting", r.name, "values for", state, "\n")
      x <- terra::extract(r[[r.name]], geo.df)[[r.name]]
      geo.df[[gsub(paste0("_", state), "", r.name)]] <- x
    }
    
    # Update crs back
    geo.df <- st_transform(geo.df, 4326)
    
    # Fix names; Filter NA values
    r.names <- gsub(paste0("_", state), "", r.names)
    names(geo.df) <- c("common.name", "observations", "geometry", r.names)
    geo.df$state <- state
    coords <- st_coordinates(geo.df) %>% as.data.frame() %>% setnames(c("lon", "lat"))
    geo.df <- geo.df %>%
      cbind(coords) %>%
      filter(dplyr::if_all(r.names, ~!is.na(.))) %>%
      select(common.name, state, lon, lat, everything()) %>%
      mutate(urban_imperviousness = urban_imperviousness %>% 
               as.character() %>% 
               as.numeric(),
             
      ) %>%
      suppressWarnings() 
    
    saveRDS(geo.df, out.file.all)
    cat("--------------\n")
  }
} 
```
