---
layout: post
title: "SDM Benchmark Study Part 6: Re-Sampling Pseudo-Absence Points Using Iterative Sampling and BIOCLIM"
permalink: /blog_posts/sdm-benchmark-study-part-6-resampling-pseudoabsence-points
---

## Overview

This part of the study presents a revised approach to selecting pseudo-absence points,
addressing some of the challenges encountered in part 5. An outline of the work done is
presented below:

1. Setup (load data, train/test splits)
2. BIOCLIM (by species, state)
  a. For each raster, output a binary raster that is either (1) for binary 
     rasters, equal to the mode of the raster, or (2) for continuous numeric
     rasters, less than the 10th percentile or greater than the 90th percentile
  b. Get the sum of all of the rasters
  c. Identify suitable pseudo-absence sampling regions as (1) points where the
     BIOCLIM sum is less than the median BIOCLIM sum, and (2) points that are
     at least 5 kilometers away from an observation point
3. Sample Pseudo-Absence Points from suitable regions
  a. Sample the same number of points as observations, except with a minimum of 
     200 points
  b. Resample (without replacement) 10 times if possible. In some cases when
     species are more widely distributed in a smaller region, there aren't 
     enough points available; In these cases, resample the maximum allowed times
  c. Split into training/test sets
4. Re-Fit LASSO Models with Updated Pseudo-Absence Points
  a. Iterate through each pseudo-absence resampling and training set fit LASSO 
     models for each (each with the same presence training set)
  b. Identify variable importance (use all models, sort by importance)
  c. Predict the probabilities of each resampled point; Select those points 
     that are most likely to be absence points and save those as the final
     pseudo-absence dataset for modeling (select from training and test sets,
     corresponding to the number of training/test presence points)
5. Using the updated datasets, complete pre-processing for IPP models (i.e.,
   set up `spatstat` data objects)

## Setup

Setup is pretty standard now for each part of the study; for more in-depth details, 
refer to prior phases of the study.

```

# Source jquery
htmltools::tags$head(
  htmltools::tags$script(
    src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.4/jquery.min.js"
  )
)

# Load libraries
library(sf)
library(terra)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
library(data.table)
library(knitr)
library(purrr)
library(caret)
library(tidyr)

# Set seed for splitting
set.seed(19)

# Load dataset containing observation points, 
# but we can remove the "absence" points
df <- setDT(readRDS("artifacts/final_data/final_data.rds"))[
  observations > 0, .(common.name, state, lon, lat, geometry)]

# Define some global variables that will be referenced throughout the modeling 
states <- sort(unique(df$state))
species <- sort(unique(df$common.name))
spec.state <- expand.grid(species=species, 
                          state=states, 
                          stringsAsFactors=F)
                          
# Load Rasters
r.list <- states %>%
  set_names() %>%
  purrr::map(~rast(file.path("artifacts/final_rasters", 
                             paste0(.x, ".tif"))))

# Pre-processed rasters
# Define min/max scaling function for rasters
min.max.scale.raster <- function(r, x=NULL, na.rm=T) {
  # If training data is provided
  if (!is.null(x)) {
    min.x <- min(x, na.rm=na.rm)
    max.x <- max(x, na.rm=na.rm)
    # Or just scale based on full raster values
  } else {
    min.x <- min(values(r), na.rm=na.rm)
    max.x <- max(values(r), na.rm=na.rm)
  }
  (r - min.x) / (max.x - min.x)
}

# Get model or other object from cache if it has been saved before
get.object <- function(obj, file.name, obj.path, 
                       read.func=readRDS, save.func=saveRDS, ...) {
  f.path <- file.path(obj.path, file.name)
  if (!dir.exists(obj.path)) {
    dir.create(obj.path)
  }
  # Object cache check
  if (file.exists(f.path)) {
    obj <- read.func(f.path)
  } else {
    save.func(obj, f.path, ...)
  }
  obj
}

# Get training/testing spatstat data, and 
# pre-process (scale) rasters 
r.list <- states %>%
  set_names() %>%
  purrr::map(function(st) {
    get.object(obj = {
      # Get raster
      r <- r.list[[st]]
      # Compute filtered & pre-processed rasters
      names(r) %>%
        set_names() %>%
        purrr::map(function(.x) {
          cat("Pre-processing", .x, "data in", st, "\n")
          r.layer <- r[[.x]]
          if (length(unique(values(r.layer))) > 2) {
            # Scale the data
            return(min.max.scale.raster(r.layer))
          } else {
            return(r.layer)
          }
        }) %>%
        terra::rast()
    },
    file.name=paste0(st, ".tif"),
    obj.path="artifacts/preprocessed_rasters_updated",
    read.func=terra::rast,
    save.func=terra::writeRaster)
  })

```

### Split into Training/Test Sets

As was done in Part 2 of the study, the data is split into a training set and a test set. The data is stratified based on latitude, longitude, species, and state to ensure that the distribution of these variables is consistent between the train and test sets. This stratification helps maintain representative samples and mitigates potential biases.

```

stratified.split.idx <- function(df, p=0.7, lat.lon.bins=25) {
  # Cut along lat/lon values to create grids (lat.bin & lon.bin)
  # lat.lon.bins is the number of divisions you want
  df$lat.bin <- cut(df$lat, breaks=lat.lon.bins, labels = F)
  df$lon.bin <- cut(df$lon, breaks=lat.lon.bins, labels = F)
  
  # Create a new variable combining the stratification variables
  df %>%
    # stratify on lat/lon bins, species, state
    mutate(strata = paste(lat.bin, lon.bin, common.name, state)) %>%
    pull(strata) %>%
    # Create the data partitions
    createDataPartition(., p = p, list = F) %>%
    suppressWarnings()
}

prepare.data <- function(df, p=.7, lat.lon.bins=25) {
  train.index <- stratified.split.idx(df, p=p, lat.lon.bins = lat.lon.bins)
  df.train <- df[train.index, ]
  df.test <- df[-train.index, ]
  
  list(train = df.train, 
       test = df.test,
       index = train.index)
}

# Get train/test indices
train.test <- prepare.data(df, .7)

# Split datatset
df.train <- df[train.test$index,] 
df.test <- df[-train.test$index,]

```

## Updated Pseudo-Absence Selection Process

### BIOCLIM

This section outlines a method for assessing suitable habitats for the bird species
using the rasters. For each species and state, it extracts values from the raster
layersat locations where the species has been observed. Based on these observed values,
a suitability mask is created using either mode (for binary layers) or percentile-based
thresholds (for continuous layers).

```
                         
bioclim.suitability <- function(sf.subset, r) {
  
  # Remove spatial layers (lat, lon, etc.)
  r <- r[[names(r)[!names(r) %in% 
                     c("lat", "lat.sqrd", "lon", 
                       "lon.lat.interaction", "lon.sqrd")]]]
  
  # Extract values from rasters at presence locations
  values.at.presences <- terra::extract(r, st_transform(sf.subset, crs=crs(r)))
  
  # Determine suitability masks based on unique values
  suitable.masks <- names(r) %>%
    set_names() %>%
    purrr::map(function(name) {
      unique.vals <- unique(values.at.presences[[name]])
      if (length(unique.vals) == 2) {
        # Binary method
        mode.val <- as.numeric(names(sort(table(values.at.presences[[name]]), 
                                          decreasing = T)[1]))
        return(r[[name]] == mode.val)
      } else {
        # 10th and 90th percentile method
        lower.bound <- quantile(values.at.presences[[name]], 0.1, na.rm = T)
        upper.bound <- quantile(values.at.presences[[name]], 0.9, na.rm = T)
        return((r[[name]] >= lower.bound) & (r[[name]] <= upper.bound))
      }
    }) %>% rast()
  
  # Combine individual masks to get areas that have the highest suitablity
  sum(suitable.masks)
}

bcr.dir <- "artifacts/bioclim_suitable_rast"
if (!dir.exists(bcr.dir )) dir.create(bcr.dir )
bioclim.list <- purrr::map(1:nrow(spec.state), function(i) {
  st=spec.state[i,]$state
  spec=spec.state[i,]$species
  sf.subset <- sf::st_as_sf(df[state == st & common.name==spec], 
                            crs=4326, sf_column_name = "geometry")
  f <- file.path(bcr.dir, paste0(st, "_", spec, ".tif"))
  if (!file.exists(f)) {
    suitable.areas <- bioclim.suitability(sf.subset, r.list[[st]])
    terra::writeRaster(suitable.areas, f)
  } else {
    suitable.areas <- rast(f)
  }
  list(
    bioclim=suitable.areas,
    state=st,
    species=spec,
    df=sf.subset
  )
})

```

### Combine BIOCLIM with Masked Sampling Regions

This section is focused on integrating the BIOCLIM suitability layers with pre-defined
sampling masks (used in the prior pseudo-absence selection process). For each species
and state, it applies a 5-kilometer resolution mask to restrict areas from which
pseudo-absence points can be sampled. The combined suitability mask is then further
refined by using the median BIOCLIM value to delineate regions suitable for
pseudo-absence point selection.

```

# Load 5k resolution mask rasters
mask.dir <- "artifacts/masks_5k"
masks <- states %>%
  set_names() %>%
  purrr::map(function(st) {
    species %>%
      set_names() %>%
      purrr::map(function(spec) {
        fname <- file.path(mask.dir, paste0(st, "_", spec, ".tif"))
        rast(fname)
    }) %>% rast()
  })

# Load the original sample regions for pseudo-absence points
pa.regions.dir <- "artifacts/pseudo_absence_regions"
pa.regions <- states %>%
  set_names() %>%
  purrr::map(function(st) {
    species %>%
      set_names(species) %>%
      purrr::map(function(spec) {
        fname <- file.path(
          pa.regions.dir,
          gsub(" |\\-", "_", paste0(paste(st, tolower(spec), sep="_"), ".tif"))
        )
        rast(fname)
      }) %>%
      rast()
  })

purrr::walk(1:nrow(spec.state), function(i) {
  st <- spec.state[i,]$state
  spec <- spec.state[i,]$species
  msk <- masks[[st]][[spec]]
  pa.region <- pa.regions[[st]][[spec]]
  idx <- purrr::detect_index(bioclim.list, ~.x$state == st & 
                        .x$species == spec)
  bioclim.list[[idx]]$mask <<- msk
  bioclim.list[[idx]]$pseudo.absence.region <<- pa.region
  qnt.5 <- quantile(values(bioclim.list[[idx]]$bioclim), .5, na.rm=T)
  bioclim.list[[idx]]$qnt.5 <<- qnt.5
  bioclim.list[[idx]]$final.pseudo.absence.region <<-
    bioclim.list[[idx]]$bioclim <= qnt.5 & terra::lapp(pa.region, function(.x) {
    ifelse(is.na(.x), 0, 1)}) == 1
})
   
```

### Visualize Suitable Sampling Regions

Below, visualizations of the regions deemed suitable for pseudo-absence point sampling
are generated. Multiple plots are generated for each species-state combination, illustrating the original "Sum BIOCLIM" suitability, suitable sampling regions based 
on BIOCLIM regions under the threshold (median), and the final sampling regions (same
as the previous, but with 5 kilometer masks removing areas around actual observation
points).

```

# Function to plot a raster with ggplot
plot.raster.gg <- function(r, fill, scale.c=T) {
  df <- terra::as.data.frame(r, xy = T)
  p <- ggplot(df, aes(x = x, y = y, fill = !!sym(fill))) + 
    geom_raster() +
    coord_cartesian()
  if (scale.c) p <- p + scale_fill_viridis_c()
  p
}

plot.to.svg <- function(p, width, height) {
  # Create a temporary file
  tmpfile <- tempfile(fileext = ".svg")
  # Save the ggplot to this file
  ggsave(filename = tmpfile, p, width=width, height=height)
  # Use img() function from htmltools to display it within a div
  img.src <- paste0("data:image/svg+xml;base64,", base64enc::base64encode(tmpfile))
  file.remove(tmpfile)
  img.src
}

# Create plots for each combination
plots.list <- purrr::map(1:nrow(spec.state), function(i) {
  st <- spec.state[i,]$state
  spec <- spec.state[i,]$species
  r <- bioclim.list[[i]]$bioclim
  qnt.5 <- bioclim.list[[i]]$qnt.5
  fpar <- bioclim.list[[i]]$final.pseudo.absence.region
  width <- case_when(st == "VT"~5,
                     st == "NC"~16,
                     st == "CO"~12,
                     st == "OR"~8,
                     T~8)
  height <- case_when(st == "VT"~7,
                     st == "NC"~8,
                     st == "CO"~9,
                     st == "OR"~8,
                     T~8)
  p1 <- (plot.raster.gg(r, "sum") + 
    labs(title = paste0("BIOCLIM Suitability for\n", spec, " in ", st),
         fill="BIOCLIM Sum")) %>%
    plot.to.svg(p=., width=width, height=height) %>%
    get.object(obj=., 
               paste0(spec, "_", st, "_1.rds"),
               "artifacts/bioclim_suitable_plots")
  
  p2 <- (plot.raster.gg(r < qnt.5, "sum", F) + 
    labs(title = paste0("BIOCLIM Suitable sample regions for\n", 
                        spec, " in ", st),
         fill="BIOCLIM Suitable\nPseudo-absence Region")) %>%
    plot.to.svg(p=., width=width, height=height) %>%
    get.object(obj=.,paste0(spec, "_", st, "_2.rds"),
               "artifacts/bioclim_suitable_plots")
  
  p3 <- (plot.raster.gg(fpar, "sum", F) + 
    labs(title = paste0("Final Suitable sampleregions for\n", 
                        spec, " in ", st),
         fill="Suitable Sampling\nRegions")) %>%
    plot.to.svg(p=., width=width, height=height) %>%
    get.object(obj=., paste0(spec, "_", st, "_3.rds"),
               "artifacts/bioclim_suitable_plots")
  plots <- htmltools::tags$div(
    id=paste0(st, "_", gsub(" ", "-", spec), "_bioclim_plots"),
    style=paste0("padding:5px; overflow:auto; width:100%;", 
                 "min-height:1400px; ", 
                 ifelse(i==1, "", "display:none;")),
    htmltools::tags$div(
      style="margin:5px;",
      htmltools::tags$img(src=p1)
    ),
    htmltools::tags$div(
      style="margin:5px;",
      htmltools::tags$img(src=p2)
    ),
    htmltools::tags$div(
      style="margin:5px;",
      htmltools::tags$img(src=p3)
    )
  )
})

htmltools::div(
  htmltools::tags$script(
    '$(document).ready(function(){
        $("#state_selector").change(function(){
          var selectedState = $(this).val();
          var selectedSpecies = $("#species_selector").val();
          // Hide all raster plots
          $("[id$=_bioclim_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_bioclim_plots").show();
        });
        $("#species_selector").change(function(){
          var selectedState = $("#state_selector").val();
          var selectedSpecies = $(this).val();
          // Hide all raster plots
          $("[id$=_bioclim_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_bioclim_plots").show();
        });
      });'
  ),
  htmltools::tags$div(
    style="display:flex; flex-direction:row;",
    htmltools::tags$div(
      style="margin-right:5px;",
      htmltools::tags$select(id='state_selector',
                             style="font-size:17px;",
                             lapply(states, function(s) {
                               htmltools::tags$option(value=s, s)
                             }))
    ),
    htmltools::tags$div(
      htmltools::tags$select(id='species_selector',
                             style="font-size:17px;",
                             lapply(species, function(s) {
                               htmltools::tags$option(value=gsub(" ", "-", s), s)
                             }))
    )
  ),
  plots.list
)
```

### Sampling (and Re-Sampling) Pseudo-Absences 

This section focuses on the extraction of pseudo-absence points from regions identified as suitable for each bird species in distinct states. A multiple-resampling technique is employed to ensure robustness in the sample. Points are randomly selected without replacement from the non-masked areas of the suitability raster, with a minimum threshold of 200 points and up to 10 resampling iterations. These points are subsequently partitioned into training and test sets to facilitate model evaluation.

```

get.samples <- function(gdf, sample.idx, spec, st) {
  samples <- gdf[sample.idx,] %>%
    mutate(common.name = spec, 
           state = st, 
           lon = NA, 
           lat = NA,
           observations=0)
  
  # Populate lon and lat values:
  coords <- st_coordinates(samples)
  samples$lon <- coords[, "X"]
  samples$lat <- coords[, "Y"]
  samples
}

# Function to sample n points from the non-masked parts
sample.inverse.mask <- function(r.final, n, 
                                spec, st,
                                sample.crs=4326, 
                                min.n=200,
                                resamples=10,
                                train.sample.size=0.7) {
  
  # Convert the raster to SpatialPoints
  gdf <- terra::as.points(r.final) %>%
    st_as_sf() %>%
    st_transform(crs=sample.crs)
  if (nrow(gdf) > 0) {
    gdf <- gdf %>%
      filter(!is.na(layer)) %>%
      select(geometry)
  } else {
    return(gdf)
  }
  
  # Set to min.n size if n < min.n
  if (n < min.n) n <- min.n
  # Last n check
  if (n > nrow(gdf)) n <- nrow(gdf)
  # Make sure there is sufficient available sample points
  resamples <- case_when(floor(nrow(gdf) / n) == 0 ~ 1, 
                         floor(nrow(gdf) / n) > resamples ~ resamples,
                         T ~ floor(nrow(gdf) / n))
  cat("\tGetting", resamples, "sample(s) of", n, "for", spec, "in", st, "\n")
  if (resamples > 1) {
    sample.indices <- sample(1:nrow(gdf), n*resamples, replace = F)
    indices <- cut(sample.indices, breaks=resamples, labels = F)
  }
  # Resample multiple times
  resampled.data <- lapply(1:resamples, function(i) {
    if (resamples > 1) {
      sample.idx <- sample.indices[indices == i]
    } else {
      sample.idx <- 1:n
    }
    train.idx <- 1:floor(train.sample.size*length(sample.idx))
    test.idx <- (max(train.idx)+1):length(sample.idx)
    train.idx <- sample.idx[train.idx]
    test.idx <- sample.idx[test.idx]
    train.samples <- get.samples(gdf, train.idx, spec, st)
    test.samples <- get.samples(gdf, test.idx, spec, st)
    
    list(train.indices = train.idx,
         test.indices = test.idx,
         train.samples = train.samples,
         test.samples = test.samples)
  })
  
  # Separate the sample indices and samples into two lists
  train.indices.list <- lapply(resampled.data, `[[`, "train.indices")
  test.indices.list <- lapply(resampled.data, `[[`, "test.indices")
  train.list <- lapply(resampled.data, `[[`, "train.samples")
  test.list <- lapply(resampled.data, `[[`, "test.samples")
  
  list(train.indices = train.indices.list, 
       test.indices = test.indices.list,
       train.samples = train.list,
       test.samples = test.list)
}

# Get totals by species and state
totals <- df[, .(.N), by=.(state, common.name)]
n.resamples <- 10
train.sample.size <- 0.7
test.sample.size <- 1 - train.sample.size

resampled.pa <- states %>%
  set_names() %>%
  purrr::map(function(st) {
    species %>%
      set_names() %>%
      purrr::map(function(spec) {
        n <- totals[state == st & common.name == spec] %>% pull(N)
        get.object({
          idx <- purrr::detect_index(bioclim.list, ~.x$state == st & 
                                       .x$species == spec)
          r.final <- bioclim.list[[idx]]$final.pseudo.absence.region
          cat(idx, "- Generating ~", n, "pseudo-absence points resampled",
              n.resamples, "times for the", spec, "in", st, "\n")
          sample.inverse.mask(
            r.final, spec, st, n=n, 
            sample.crs=4326, 
            resamples=n.resamples,
            train.sample.size=train.sample.size)
        },
        paste0(st, "_", spec, ".rds"),
        "artifacts/resampled_pseudo_absences_with_bioclim")
      })
  })

```

### Combine Presence/Pseudo-Absence Data

Each sample of pseudo-absence points should be concatenated with the corresponding
presence data (i.e., training/test set, by species and state). The result is up to
20 datasets (10 training, 10 test) for each species/state combination. Note that
the presence data is NOT resampled, and is thus recycled for each set of
corresponding pseudo-absence data.

```

get.pres.df <- function(df, st, spec, r) {
  sf.pres <- sf::st_as_sf(
    df[state == st & common.name==spec, 
       .(geometry, common.name, state)], 
    crs=4326, sf_column_name = "geometry")
  # Extract values from rasters at presence locations
  pres.vals <- sf.pres %>%
    cbind(terra::extract(r, st_transform(sf.pres, crs=crs(r)))) %>%
    select(-ID) %>%
    st_transform(crs(sf.pres)) %>%
    filter(dplyr::if_all(names(.), ~!is.na(.))) %>%
    suppressWarnings()
}

get.abs.df <- function(df.list, r, sf.pres) {
  abs.vals <- df.list %>%
    purrr::map(~{
      .x %>%
        cbind(terra::extract(r, st_transform(.x, crs=crs(r)))) %>%
        select(-ID) %>%
        st_transform(crs(sf.pres)) %>%
        filter(dplyr::if_all(names(.), ~!is.na(.))) %>%
        suppressWarnings()
    })
}

dfs <- states %>%
  set_names() %>%
  purrr::map(function(st) {
    r <- r.list[[st]]
    species %>%
      set_names() %>%
      purrr::map(function(spec) {
        cat("Getting values for", spec, "in", st, "\n")
        get.object({
          sf.pres.train <- get.pres.df(df.train, st, spec, r)
          sf.pres.test <- get.pres.df(df.test, st, spec, r)
          abs.vals.train <- get.abs.df(resampled.pa[[st]][[spec]]$train.samples, 
                                       r, sf.pres.train)
          abs.vals.test <- get.abs.df(resampled.pa[[st]][[spec]]$test.samples, 
                                      r, sf.pres.test)

          list(
            presence.train=sf.pres.train,
            presence.test=sf.pres.test,
            pseudo.absence.train=abs.vals.train,
            pseudo.absence.test=abs.vals.test
          )
        },
        paste0(spec, "_", st, ".rds"),
        "artifacts/resampled_pa_with_bioclim_values")
      })
  })
```

### Post-Processing Checks

To confirm that the number of presence/absence points in each of the created 
datasets is approximately correct, create a table of all of the counts.

```

# Compare the total number of values in each dataset (just spot-checking here)
purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  d <- dfs[[st]][[spec]]
  prs.trn <- nrow(d$presence.train)
  prs.tst <- nrow(d$presence.test)
  p.abs.trn <- purrr::map_dbl(d$pseudo.absence.train, ~nrow(.x))
  p.abs.tst <- purrr::map_dbl(d$pseudo.absence.test, ~nrow(.x))
  tibble(
    state=st,
    species=spec,
    presence.train=prs.trn,
    presence.test=prs.tst,
    min.pseudo.abs.train=min(p.abs.trn),
    max.pseudo.abs.train=max(p.abs.trn),
    min.pseudo.abs.test=min(p.abs.tst),
    max.pseudo.abs.test=max(p.abs.tst)
  )
})

```

## Updated LASSO Models

### Identify "Redundant" Rasters

As previously stated, many of the features are redundant. In particular, the
"engineered" features were meant to improve upon the original features, and reduce
multi-collinearity. Based on the original results of the feature selection process, it
was shown that the new variables weren't always necessarily "better". Below is the list
of raster layers that were determined should be removed from each of the rasters in
order to minimize multi-collinearity and maximize model performance. 

```

# Define "redundant" rasters/features
redund.feat <- c("open_water", "developed_open_space", "developed_low_intensity",
                 "developed_medium_intensity", "developed_high_intensity",
                 "hay_pasture", "cultivated_crops", "wetlands", "forest",
                 "lon.lat.interaction", "waterbody", "urban_imperviousness", "tmax",
                 "tmin", "dem", "Fall_NDVI", "Spring_NDVI", "Summer_NDVI",
                 "Winter_NDVI", "lat.sqrd", "lon.sqrd")

```

### Re-Fit LASSO Models with Updated Pseudo-Absence Points

As was done in part 2, a GLM model is fit for each of the species/state subsets using the training data, but this time with the refined raster subset and repeated for
each iteration of sampled pseudo-absence data.

```

# Load variable importance from fitted LASSO models
lasso.model.path="artifacts/models/lasso_3_fs"

species %>% 
  purrr::walk(function(spec) {
    states %>% 
      purrr::walk(function(st) {
        # Define the control for the train function
        ctrl <- trainControl(method = "cv", number = 5)
        
        cat("Fitting LASSO model for", spec, "in", st, "\n")
        spec.df <- dfs[[st]][[spec]]$presence.train %>%
          setDT()
        spec.df <- spec.df[
          common.name == spec & state == st][
            , `:=` (state=NULL, common.name = NULL, geometry=NULL, 
                    presence=factor(T, levels=c("TRUE", "FALSE")))]
        
        # Remove any columns where all values are the same
        .remove <- c(names(which(sapply(spec.df, function(x) {
          length(unique(x)) <= 1}))),
          redund.feat, "lat", "lon", "observations") %>% 
          unique()
        .remove <- .remove[.remove %in% names(spec.df)]
        .remove <- .remove[.remove != "presence"]
        if (!is_empty(.remove)) {
          spec.df <- spec.df %>% dplyr::select(-.remove)
        }
        train.samples <- length(dfs[[st]][[spec]]$pseudo.absence.train)
        purrr::walk(
          1:train.samples, function(i) {
            fname <- paste0(tolower(gsub(" ", "_", spec)), "_", st,
                            "_", i, "_regr_l1.rds")
            if (!file.exists(fname)) {
              if (i == 1) gc()
              cat("Fitting LASSO model for", i, "/", train.samples,
                  "iterations, for", spec, "observations in", st, "\n")
              abs.df <- dfs[[st]][[spec]]$pseudo.absence.train[[i]] %>% setDT()
              abs.df <- abs.df[
                common.name == spec & state == st][
                  , `:=` (state=NULL, common.name = NULL, geometry=NULL, 
                          presence=factor(F, levels=c("TRUE", "FALSE")))]
              if (!is_empty(.remove)) {
                abs.df <- abs.df %>% dplyr::select(names(spec.df))
              }
              train.df <- data.table::rbindlist(list(spec.df, abs.df))
              
              # LASSO (L1); Elastic Net, where alpha = 1
              fit <- get.object(
                train(presence ~ (.)^2, 
                      data = train.df, 
                      method = "glmnet",
                      family = "binomial",
                      trControl = ctrl,
                      tuneGrid = expand.grid(alpha = 1, 
                                             lambda = seq(0, 1, by = 0.1)),
                      metric = "Accuracy"),
                file.name=fname,
                obj.path=lasso.model.path)
              coef.df <- coef(fit$finalModel, s = fit$bestTune$lambda) %>%
                as.matrix() %>%
                as.data.frame()
              # Remove the intercept
              coef.df <- coef.df[-1, , drop = F]
              if (sum(coef.df$s1, na.rm=T) == 0) {
                # Remove file
                file.remove(file.path(lasso.model.path, fname))
                # Re-train model, this time with variable alpha in 
                # the train grid
                fit <- get.object(
                  train(presence ~ (.)^2, 
                        data = train.df, 
                        method = "glmnet",
                        family = "binomial",
                        trControl = ctrl,
                        tuneGrid = expand.grid(
                          alpha = seq(0, 1, by = 0.1), 
                          lambda = seq(0, 1, by = 0.1)),
                        metric = "Accuracy"),
                  file.name=fname,
                  obj.path=lasso.model.path)
                coef.df <- coef(fit$finalModel, s = fit$bestTune$lambda) %>%
                  as.matrix() %>%
                  as.data.frame()
                # Remove the intercept
                coef.df <- coef.df[-1, , drop = F]
                if (sum(coef.df$s1, na.rm=T) == 0) {
                  file.remove(file.path(lasso.model.path, fname))
                  cat("\tERROR: Try another method for obtaining coefs for",
                      spec, "in", st, "\n")
                }
              }
            }
          })
      })
  })
    
```
    
### Get Updated Variable Importance

Select the most important variables of each model, by species/state (combine the
results of each of the models fit for each species/state if there are multiple).
    
```
    
# Define min/max scaling function for rasters
min.max.scale <- function(x, na.rm=T) {
  min.x <- min(x, na.rm=na.rm)
  max.x <- max(x, na.rm=na.rm)
  (x - min.x) / (max.x - min.x)
}

get.var.imp <- function(st, spec, dir="artifacts/models/lasso_3_fs") {
  files <- list.files(dir) %>%
    .[grepl(paste(tolower(gsub(" ", "_", spec)), st, sep="_"), .)] %>%
    file.path(dir, .)
  if (length(files) == 0) return(
    tibble(
      common.name=character(0),
      state=character(0),
      variable=character(0),
      importance=numeric(0),
      wt=numeric(0),
      weighted.imp=numeric(0)
    )
  )
  var.imp <- purrr::map_df(files, ~{
    fit <- readRDS(.x)
    coef.df <- coef(fit$finalModel, s = fit$bestTune$lambda) %>%
      as.matrix() %>%
      as.data.frame()
    # Remove the intercept
    coef.df <- coef.df[-1, , drop = F]
    # Create a data frame of variable importance
    var.importance <- tibble(
      common.name = spec,
      state = st,
      variable = rownames(coef.df),
      importance = abs(coef.df[,1])
    ) %>%
      # Rank variables by importance
      arrange(state, common.name, -importance, variable) %>%
      # Only look at variables where imp. > 0
      filter(importance > 0)
  }) %>%
    mutate(n=n()) %>%
    group_by(common.name, state, variable, n, .drop=F) %>%
    summarize(importance=median(importance), 
              wt=n()) %>%
    ungroup() %>%
    mutate(wt=wt/n) %>%
    select(-n) %>%
    mutate(weighted.imp=min.max.scale(wt*importance)) %>%
    arrange(-weighted.imp, variable) 
}

var.imp <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  get.var.imp(st, spec)
})

var.imp

```

### Select Final Pseudo-Absence Points Using Model Predictions

Finally, the pseudo-absence points will be finalized by using the model(s) for each 
state/species to estimate the probability each of the original points are actually
absence points. The results are sorted, and the most likely points are selected as the
final pseudo-absence points. This is done for both the training and the test sets.

```

preds.abs.df <- function(df.list, tot, spec, st, files, .train=T, verbose=F) {
  preds <- purrr::map(1:tot, function(i) {
    abs.df <- df.list[[i]] %>%
      setDT() %>%
      .[common.name == spec & state == st] %>%
      .[, `:=` (state=NULL, common.name = NULL, 
                geometry=geometry,
                presence=factor(F, levels=c("TRUE", "FALSE")))]
    geom <- abs.df$geometry
    abs.df[, geometry := NULL]
    if (!is_empty(.remove)) {
      abs.df <- abs.df %>% dplyr::select(-.remove)
    }
    preds.j <- files %>%
      purrr::imap(function(f, j) {
        if (verbose) cat(i, "/", tot, "Getting predictions for model", j, 
                         "for", spec, "in", st, "\n")
        get.object({
          fit <- readRDS(f)
          yhat <- predict(fit, newdata=abs.df, type="prob") %>%
            setDT() %>%
            .[, .(`FALSE`)] %>%
            setnames("FALSE", paste0("abs.prob.", j))},
          paste0(st, "_", spec, "_", i, "_", j, ".rds"),
          ifelse(.train, "artifacts/pa_probs/train", "artifacts/pa_probs/test"))
      }) %>% do.call("cbind", .)
    
    preds.j[, `:=` (med=apply(.SD, 1, median), 
                    id=.I, sample=i, train=.train,
                    geometry=geom, state=st,
                    common.name=spec)] %>% 
      setorderv(x=., cols="med", order=-1)
    preds.j[, .(state, common.name, train, med, geometry)]
  }) %>% rbindlist()
}

if (!dir.exists("artifacts/pa_probs")) dir.create("artifacts/pa_probs")
get.pa.prob <- function(st, spec, dfs,
                        dir="artifacts/models/lasso_3_fs", 
                        verbose=F) {
  files <- list.files(dir) %>%
    .[grepl(paste(tolower(gsub(" ", "_", spec)), st, sep="_"), .)] %>%
    file.path(dir, .)
  if (length(files) == 0) return(NULL)
  
  # Get columns that should be removed (based on original presence data)
  spec.df <- dfs[[st]][[spec]]$presence.train %>%
    setDT() %>%
    .[common.name == spec & state == st] %>%
    .[, `:=` (state=NULL, common.name = NULL, geometry=NULL, 
              presence=factor(T, levels=c("TRUE", "FALSE")))]
  .remove <- c(names(which(sapply(spec.df, function(x) {
    length(unique(x)) <= 1}))),
    redund.feat, "lat", "lon", "observations") %>% 
    unique()
  .remove <- .remove[.remove %in% names(spec.df)]
  .remove <- .remove[.remove != "presence"]
  tot.train <- length(dfs[[st]][[spec]]$pseudo.absence.train)
  tot.test <- length(dfs[[st]][[spec]]$pseudo.absence.test)
  preds.train <- preds.abs.df(dfs[[st]][[spec]]$pseudo.absence.train, 
                              tot.train, spec, st, files, .train=T)
  preds.test <- preds.abs.df(dfs[[st]][[spec]]$pseudo.absence.test, 
                              tot.test, spec, st, files, .train=F)
  rbindlist(list(preds.train, preds.test))
}

purrr::walk(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  cat(spec, st, "\n")
  tryCatch({get.object(
    {get.pa.prob(st, spec, dfs, 
                 verbose=T) %>% 
        setorderv(c("train", "med"), order=c(-1, -1))},
    paste0(st, "_", spec, ".rds"),
    "artifacts/pa_prob_results"
  )}, error=function(e) NULL)
})

```


```

purrr::walk(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  f <- file.path("artifacts/pa_prob_results", paste0(st, "_", spec, ".rds"))
  if (file.exists(f)) {
    get.object({
      dt <- readRDS(f)[, .(common.name, state, train, geometry)]
      pres.train <- df.train[common.name == spec & state == st] %>% 
        sf::st_as_sf(x=., crs=4326) %>%
        mutate(presence=factor(T, levels=c("TRUE", "FALSE")))
      pres.test <- df.test[common.name == spec & state == st] %>% 
        sf::st_as_sf(x=., crs=4326) %>%
        mutate(presence=factor(T, levels=c("TRUE", "FALSE")))
      ntrain <- pres.train %>% 
        nrow() %>%
        ifelse(. < 500, 500, .)
      ntest <- pres.test %>% 
        nrow() %>%
        ifelse(. < 200, 200, .)
      cat(spec, st, ntrain, ntest, "\n")
      trn <- dt[train==T][1:ntrain] %>% 
        sf::st_as_sf(x=., crs=4326) %>%
        cbind(data.frame(sf::st_coordinates(.$geometry[1:nrow(.)])) %>% 
                set_names(c("lon", "lat"))) %>%
        mutate(presence=factor(F, levels=c("TRUE", "FALSE"))) %>%
        select(names(pres.train))
      trn <- do.call(what="rbind", args=list(trn, pres.train))
      tst <- dt[train==F][1:ntest] %>% 
        sf::st_as_sf(x=., crs=4326) %>%
        cbind(data.frame(sf::st_coordinates(.$geometry[1:nrow(.)])) %>% 
                set_names(c("lon", "lat"))) %>%
        mutate(presence=factor(F, levels=c("TRUE", "FALSE"))) %>%
        select(names(pres.test))
      tst <- do.call(what="rbind", args=list(tst, pres.test))
      
      list(train=trn, test=tst)
    }, paste0(st, "_", spec, ".rds"),
    "artifacts/train_test_pseudo_abs_updated")
  }
})

```

## Pre-Processing IPP Model Covariates

Finally, as was done previously, the updated data are preprocessed and saved to
`spatstat.geom::quadscheme` objects so that updated IPP models can be fit using
the resampled data.

```

# Ensure counter-clockwise direction
is.ccw <- function(p) {
  tryCatch({owin(poly=p); T}, error=function(e) F)
}

# Check that all polygons were traversed in the right direction
ensure.ccw <- function(polygon.list) {
  lapply(polygon.list, function(p) {
    # Check the first polygon (outer boundary)
    if (!is.ccw(p)) {
      p$x <- rev(p$x)
      p$y <- rev(p$y)
    }
    p
  })
}

# Function to convert polygon to required list format
convert.to.list.format <- function(sf.object) {
  polygons.list <- lapply(1:nrow(sf.object), function(i) {
    sfc <- st_geometry(sf.object[i,])
    if (class(sfc)[[1]] == "sfc_MULTIPOLYGON") {
      sfc <- st_cast(sfc, "POLYGON")
    }
    lapply(sfc, function(poly) {
      list(x = st_coordinates(poly)[,1], 
           y = st_coordinates(poly)[,2])
    })
  })
  
  # If the object has only one row, we unlist one level to adhere 
  # to the described format for `ppp` objects
  if (length(polygons.list) == 1) {
    polygons.list <- polygons.list[[1]]
  }
  
  # Ensure counter-clockwise
  polygons.list <- ensure.ccw(polygons.list)

  return(polygons.list)
}

prep.ipp.data <- function(data, st, spec, 
                          r.list, just.presence=F) {
  # Filter Raster List
  r <- r.list[[st]]
  
  # filter by state & species
  ipp.df <- data %>% 
    filter(state == st & common.name == spec) %>%
    # select location point
    select(c("lon", "lat", "presence")) %>% 
    # Get the unique points, since we are not accounting 
    # for the temporal nature of the data
    unique() 
    
  # Get the first layer, set it to either NA or TRUE
  r.poly <- terra::project(x=r[[1]], 
                           y=st_crs(4326, parameters=T)$Wkt) %>% 
    lapp(function(z) ifelse(is.na(z), NA, T)) %>%
    terra::as.polygons() %>%
    # Convert to polygon
    st_as_sf()
  
  # Convert polygon to list
  r.poly.list <- convert.to.list.format(r.poly)
  
  # Get indices of points that are within the polygon
  valid.pts <- sapply(st_intersects(ipp.df, r.poly), function(x) length(x) > 0)
  
  # Filter out invalid points
  ipp.df <- filter(ipp.df, valid.pts)
  
  # Subset df by presence
  ipp.presence <- filter(ipp.df, presence == T)
  ipp.dummies <- filter(ipp.df, presence == F)
  
  # Convert the data to a ppp objects
  locations <- spatstat.geom::ppp(ipp.presence$lon, 
                                  ipp.presence$lat, 
                                  poly=r.poly.list) 
  
  if (just.presence) return(locations)
  
  dummy.locs <- spatstat.geom::ppp(ipp.dummies$lon, 
                                   ipp.dummies$lat, 
                                   poly=r.poly.list) 
  
  # Create Quadscheme
  Q <- spatstat.geom::quadscheme(locations, dummy.locs)
}

# Save training/testing spatstat data
purrr::walk(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  cat("Prepping final data for", spec, "in", st, "\n")
  # Get raster
  r <- r.list[[st]]
  data <- file.path("artifacts/train_test_pseudo_abs_updated",
                    paste0(st, "_", spec, ".rds")) %>% readRDS()
  train.data <- data$train
  test.data <- data$test
  cat("\tGetting `spatstat.geom::ppp` object with training points...\n")
  Q <- get.object(
    prep.ipp.data(train.data, st, spec, r.list),
    paste0(st, "_", spec, "_Q.rds"),
    file.path("artifacts", "train_spatstat_Q_2")
  )
  cat("\tGetting `spatstat.geom::ppp` object with testing points...\n")
  Q.test <- get.object(
    prep.ipp.data(test.data, st, spec, r.list),
    paste0(st, "_", spec, "_Q.rds"),
    file.path("artifacts", "test_spatstat_Q_2")
  )
})

```