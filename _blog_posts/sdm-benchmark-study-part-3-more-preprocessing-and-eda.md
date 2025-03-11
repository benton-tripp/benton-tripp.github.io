---
layout: post
title: "SDM Benchmark Study Part 3: More Pre-Processing and EDA"
permalink: /blog_posts/sdm-benchmark-study-part-3-more-preprocessing-and-eda
---

## Overview

The Species Distribution Model (SDM) Benchmark Study Part 3 continues an exploration of the data from from eight bird species across Colorado, Oregon, North Carolina, and Vermont, as well as the corresponding environmental rasters. Along with in-depth testing of the data ensuring that it meets model assumptions, this part of the study further refines the feature selection process and completes the pre-processes of the raster datasets.

## Setup

The data is already available from Parts 1 and 2 of the study. To begin, it will need to be loaded into the session, along with the required libraries. A helper function will also be declared, which will be used to help simplify the caching process of the analysis outputs.

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
library(spatstat)
library(plotly)

# Get model or other object from cache if it has been saved before
get.object <- function(obj, file.name, obj.path, read.func=readRDS, save.func=saveRDS, ...) {
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
```

## Loading the Data

### Load the Tabular Data

```
# Load the dataset saved in part 2 of the study 
df <- readRDS("artifacts/final_data/final_data.rds") %>% setDT()

# Define some global variables that will be referenced throughout the modeling 
states <- sort(unique(df$state))
species <- sort(unique(df$common.name))
spec.state <- expand.grid(species=species, 
                          state=states, 
                          stringsAsFactors=F)

# Get a binary target for presence/absence instead of obs. counts for MaxEnt
set(df, j="presence", value=factor(ifelse(df$observations == 0, 0, 1), levels=c(0,1)))

# Save geometry for reference, and remove from dataset
geometry <- df$geometry

df[, `:=` (geometry=NULL)]

# View output
df %>% as_tibble()
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
  <script src="{{ site.baseurl }}/assets/data/sdm-3-eda-df.js"></script>
  <script data-pagedtable-source data-global="sdm3DF" type="application/json"></script>
</div>

### Split into Training/Test Sets

As was done in Part 2 of the study, the data is split into a training set and a test set. The data is stratified based on latitude, longitude, species, state, and absence to ensure that the distribution of these variables is consistent between the train and test sets. This stratification helps maintain representative samples and mitigates potential biases.

```
# Set seed for splitting and modeling
set.seed(19)

stratified.split.idx <- function(df, p=0.7, lat.lon.bins=25) {
  # Cut along lat/lon values to create grids (lat.bin & lon.bin)
  # lat.lon.bins is the number of divisions you want
  df$lat.bin <- cut(df$lat, breaks=lat.lon.bins, labels = F)
  df$lon.bin <- cut(df$lon, breaks=lat.lon.bins, labels = F)
  
  # Create a new variable combining the stratification variables
  df %>%
    # stratify on lat/lon bins, species, state, and presence/absence
    mutate(strata = paste(lat.bin, lon.bin, common.name, state, presence)) %>%
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

### Load Rasters

All of the rasters saved in parts 1 and 2 of the study are loaded into 4 multi-layer rasters (one for each state).

```
# Load "Feature Engineered" rasters and original rasters into a 
# single multi-layer raster by state
r.list <- set_names(states) %>%
  purrr::map(~rast(c(paste0("data/final_rasters/", .x, ".tif"),
                     file.path("artifacts/feature_engineered_final", 
                               paste0(.x, ".tif")))))
```

### Visualize the Rasters

In the prior sections, not all of the raster visualizations were output to the rendered results. To gain a better understanding of each of the raster layers, an interactive plot can be generated to simplify parsing through all of the layers at each state.

```
plt.ncol <- list(
  NC=4,
  CO=5,
  VT=6,
  OR=5
)
# Get plots as encoded html objects
r.plots <- get.object(
  states %>%
    set_names() %>%
    purrr::map(function(st) {
      cat("Getting raster plots for", st, "\n")
      ncol <- plt.ncol[[st]]
      nrow <- ceiling(length(names(r.list[[st]])) / ncol)
      p <- get.object(
        # Create plots of all rasters, by state; Save plot
        purrr::map(names(r.list[[st]]), function(r.name) {
          r.df <- terra::as.data.frame(r.list[[st]][[r.name]], xy=T)
          p <- ggplot(r.df, aes(x=x, y=y, fill=!!sym(r.name))) +
            geom_raster() +
            coord_cartesian() 
          if (r.name != "NLCD_Land") p <- p + scale_fill_viridis_c()
          p + labs(title=r.name) + theme(legend.position="none")
        }) %>%
          ggarrange(plotlist=., ncol=ncol, nrow=nrow) +
          ggtitle(paste0("All Raster Layers for ", st)),
        paste0(st, "_rasters.rds"),
        "artifacts/plots"
      )
      
      # Create a temporary file
      tmpfile <- tempfile(fileext = ".svg")
      # Save the ggplot to this file
      ggsave(filename = tmpfile, p, width = 24, height = 24)
      # Use img() function from htmltools to display it within a div
      img.src <- paste0("data:image/svg+xml;base64,", base64enc::base64encode(tmpfile))
      file.remove(tmpfile)
      out <- htmltools::div(id=paste0(st, "_raster_plts"), 
                            style=ifelse(st == states[[1]], "", "display:none;"), 
                            htmltools::tags$img(src=img.src))
      gc()
      out
    }),
  "all_rasters_html.rds",
  "artifacts/plots"
)
```

```
htmltools::div(
  htmltools::tags$script(
    '$(document).ready(function(){
        $("#state_selector").change(function(){
          var selectedState = $(this).val();
          // Hide all raster plots
          $("[id$=_raster_plts]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_raster_plts").show();
        });
      });'
  ),
  htmltools::tags$select(id='state_selector',
                                       lapply(states, function(st) {
                                         htmltools::tags$option(value=st, st)
                                       })),
  r.plots
)
```

<div>
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.4/jquery.min.js"></script>
  <script>
  $(document).ready(function(){
    $("#state_selector").change(function() {
      var selectedState = $(this).val();
      // Hide all raster plots
      $("[id$=_raster_plts]").hide();
      // Show the selected raster plot
      $("#" + selectedState + "_raster_plts").show();
    });
  });
  </script>
    <select id="state_selector">
    <option value="CO">CO</option>
    <option value="NC">NC</option>
    <option value="OR">OR</option>
    <option value="VT">VT</option>
  </select>
<div id="raster-select-plts">
    <div id="CO_raster_plts" style="">
      <img src="{{ site.baseurl }}/assets/plots/sdm-eda-co.svg">
    </div>
    <div id="NC_raster_plts" style="display:none;">
      <img src="{{ site.baseurl }}/assets/plots/sdm-eda-nc.svg">
    </div>
    <div id="OR_raster_plts" style="display:none;">
      <img src="{{ site.baseurl }}/assets/plots/sdm-eda-or.svg">
    </div>
    <div id="VT_raster_plts" style="display:none;">
      <img src="{{ site.baseurl }}/assets/plots/sdm-eda-vt.svg">
    </div>
  </div>
</div>
