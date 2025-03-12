---
layout: post
title: "SDM Benchmark Study Part 5: Fitting and Testing Inhomogeneous Poisson Process Models"
permalink: /blog_posts/sdm-benchmark-study-part-5-fitting-and-testing-ipp-models
---


## Overview

This part of the study shows my first iteration of fitting the IPP models using the data
pre-processed in prior parts of the study. Because it is mostly just large bodies of code for the 
main portion, I will just give a rough outline here of the steps that are taken. Then, I will go 
into greater detail discussing the results.

### IPP Model Process Outline

- Load libraries
- Read in observation and raster data from [part 2]({{ site.baseurl }}/blog/blog_posts/sdm-benchmark-study-part-2-exploratory-analysis.html) and [part 3]({{ site.baseurl }}/blog/blog_posts/sdm-benchmark-study-part-3-more-preprocessing-and-eda.html) of the study
- Establish utility functions and variables
- Load variable importance (obtained from fitted LASSO models in [part 3]({{ site.baseurl }}/blog/blog_posts/sdm-benchmark-study-part-3-more-preprocessing-and-eda.html) of the study)
- Iterate over each state/species pair
  * Define the base number of maximum covariates and covariate interactions to keep (50)
  * Load pre-processed train/test observation data (`spatstat.geom::quadscheme` objects cached in 
    [part 3]({{ site.baseurl }}/blog/blog_posts/sdm-benchmark-study-part-3-more-preprocessing-and-eda.html) of the study)
  * Load the corresponding pixel images (`spatstat.geom::im` objects, also cached in 
    [part 3]({{ site.baseurl }}/blog/blog_posts/sdm-benchmark-study-part-3-more-preprocessing-and-eda.html)) for each of the covariates up to the maximum allowed number
  * Using the training data and pixel images, fit an IPP (GLM) Model using the Method of 
    Maximum Pseudo-Likelihood (`spatstat.model::ppm`)
  * Check for issues (model does not converge, prediction error, etc.);
    - If there is an issue, reduce the number of maximum allowed covariates
      and re-fit the model
    - Repeat this process until there are either no problems or else there
      are no more available covariates (in which case fit a model with no trend)
  * Compute the Receiver Operating Characteristic curve for the fitted point process 
    model (using `spatstat.model::roc.ppm`)
  * Generate predictions with standard errors and a confidence interval for the 
    test data
  * Using the estimates and C.I. values, estimate probabilities 
  * Generate predictions for the entire region
  * Estimate prediction intervals with standard errors
  * Find the optimal threshold using Youden's J statistic (i.e., find the threshold 
    that maximizes the J statistic)
  * Generate model accuracy metrics using the test predictions/actuals

## Setup

```
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
library(tidyr)

# Set seed for splitting and modeling
set.seed(19)

# Load the dataset saved in part 2 of the study 
df <- readRDS("artifacts/final_data/final_data.rds") %>% setDT()

# Define some global variables that will be referenced throughout the modeling 
states <- sort(unique(df$state))
species <- sort(unique(df$common.name))
spec.state <- expand.grid(species=species, 
                          state=states, 
                          stringsAsFactors=F)


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

### Load Variable Importance

```
# Load variable importance from fitted LASSO models
lasso.model.path="artifacts/models/lasso_2_fs"

var.imp <- species %>% purrr::map_df(function(spec) {
    spec.state.fit <- states %>% purrr::map_df(function(st) {
      fname <- paste0(tolower(gsub(" ", "_", spec)), "_", st, "_regr_l1.rds")
      fit <- readRDS(file.path(lasso.model.path, fname))
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
    })
  })

```

## Fit IPP Models

```
optimize.threshold <- function(roc.obj) {
  # Calculate Youden's J statistic for each threshold
  j.values <- roc.obj$fobs + (1 - roc.obj$ftheo) - 1
  # Find the threshold that maximizes the J statistic
  optimal.threshold <- roc.obj$p[which.max(j.values)]
  return(optimal.threshold)
}

# Iterate over each state/species pair
purrr::walk(1:nrow(spec.state), function(i) {
  # Define the base number of maximum covariates to keep
  covariates.keep <- 50
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  glm.results.path <- file.path("artifacts/test_results/ipp_glm_mpl",
                                paste0(spec, "_", st, "_ipp_glm_mpl.rds"))
  if (!file.exists(glm.results.path)) {
    cat("Getting IPP model for", spec, "in", st, "\n")
    
    # Load observation data
    cat("\tGetting pre-processed quadrature obj. (train)\n")
    Q <- readRDS(file.path("artifacts", "train_spatstat_Q",
                           paste(st, "_", spec, "_Q.rds")))
    cat("\tGetting quadrature obj. (test)\n")
    Q.test <- readRDS(file.path("artifacts", "test_spatstat_Q",
                                paste(st, "_", spec, "_Q.rds")))
    
    glm.base.path <- "artifacts/models/ipp_glm_mpl"
    # Regular expression pattern to match the filename
    pattern <- paste0("^", spec, "_", st, "_([0-9]{1,2})_ipp_glm_mpl\\.rds$")
    glm.path <- list.files(path = glm.base.path, pattern = pattern, full.names = T)
    if (file.exists(glm.path)) {
      fit.glm <- readRDS(glm.path)
      roc.obj <- spatstat.explore::roc(fit.glm)
      covariates.keep <- as.numeric(regmatches(glm.path, 
                                               regexpr("[0-9]{1,2}", 
                                                       glm.path)))
      if (covariates.keep > 0) {
        covariates <- coef(fit.glm) %>% names() %>% .[-1] %>%
        stringr::str_split(., pattern="\\:") %>% unlist() %>%
        unique() %>% sort() %>% 
        set_names() %>%
        purrr::map(function(.x) {
          file <- file.path("artifacts/spatstat_imgs", 
                            paste0(st, "_", .x, ".rds"))
          readRDS(file)
        })
      } else {
        covariates <- c()
      }
      
    } else {
      converged <- F
      no.pred.error <- F
      roc.obj <- NULL
      while (!converged | !no.pred.error | is.null(roc.obj)) {
        # Set path
        glm.path <- file.path("artifacts/models/ipp_glm_mpl", 
                              paste0(spec, "_", st, "_", covariates.keep,
                                     "_ipp_glm_mpl.rds"))
        if (covariates.keep > 0) {
          # Select covariates based on feature importance
          cat("\tFetching variable importance with `covariates.keep` set to", 
              covariates.keep, "\n")
          fs.df <- var.imp %>% 
            filter(state == st & common.name == spec) %>%
            mutate(var1 = purrr::map_chr(variable, 
                                         ~ stringr::str_split(.x, "\\:")[[1]][1]),
                   var2 = purrr::map_chr(variable, ~ {
                     split_result <- stringr::str_split(.x, "\\:")[[1]]
                     if(length(split_result) > 1) split_result[2] else NA_character_
                   })) %>%
            mutate(variable = ifelse(is.na(var2), 
                                     var1, 
                                     paste(var1, var2, sep = ":"))) %>%
            # Keep only pre-determined # of variables/interactions
            head(covariates.keep)
          covariates.keep <- nrow(fs.df) # in case there are fewer available
          
          if (nrow(fs.df) > 0) {
            covariates <- c(fs.df$var1, fs.df$var2) %>% 
              unique() %>% 
              sort()
          } else {
            cat("\tThere are no specified covariates for", spec, st, "\n")
          }
          
          # Load/compute filtered & pre-processed rasters
          covariates <- covariates %>%
            set_names() %>%
            purrr::map(function(.x) {
              file <- file.path("artifacts/spatstat_imgs", 
                                paste0(st, "_", .x, ".rds"))
              readRDS(file)
            })
          
          # Create formula
          .f <- paste(fs.df$variable, collapse=" + ") %>% 
            paste("~", .) %>% 
            as.formula()
          # Fit the IPP model, using the Method of Maximum PseudoLikelihood (MPL)
          # * gcontrol=list() for additional glm/gam parameters
          
          # GLM Model
          cat("\tFitting GLM Model using the Method of Maximum",
              "PseudoLikelihood...\n")
          fit.glm <- tryCatch({
            ppm(Q=Q, trend=.f, covariates=covariates, 
                rbord=.05, method="mpl") %>%
              get.object(
                obj=.,
                file.name=paste0(spec, "_", st, "_", covariates.keep,
                                 "_ipp_glm_mpl.rds"), 
                obj.path="artifacts/models/ipp_glm_mpl")}, 
            error=function(e) NULL)
          if (is.null(fit.glm)) {
            if (covariates.keep > 15) {
              covariates.keep <- covariates.keep - 5
            } else {
              covariates.keep <- covariates.keep - 1
            }
            if (covariates.keep < 0) stop("\tNo trend identified.\n")
            warning("\tThere was an error fitting the model. Trying a",
                    "different number of covariates\n")
            next
          }
          converged <- fit.glm$internal$glmfit$converged
        } else {
          cat("\tFitting GLM Model using the Method of Maximum",
              "PseudoLikelihood with no trend...\n")
          fit.glm <- ppm(Q=Q, rbord=.05, method="mpl") %>%
            get.object(
              obj=.,
              file.name=paste0(spec, "_", st, "_", covariates.keep,
                               "_ipp_glm_mpl.rds"), 
              obj.path="artifacts/models/ipp_glm_mpl")
          converged <- T
        }
        # Check for errors when predicting
        options(show.error.messages = F)
        no.pred.error <- suppressWarnings(
          tryCatch(
            expr={predict.ppm(fit.glm, se=T); T}, 
            error=function(e) {F}
          )
        )
        options(show.error.messages = T)
        # Compute the ROC curve, check for errors
        cat("\tGetting ROC...\n")
        roc.obj <- tryCatch(expr={spatstat.explore::roc(fit.glm)}, 
                            error=function(e) {NULL})
        if (!converged | !no.pred.error | is.null(roc.obj)) {
          warning("\tThe model converged or there was an error",
                  "for covariates.keep ==", covariates.keep, "\n")
          
          file.remove(glm.path)
          if (covariates.keep > 15) {
            covariates.keep <- covariates.keep - 5
          } else {
            covariates.keep <- covariates.keep - 1
          }
          if (covariates.keep < 0) stop("\tNo trend identified.\n")
        }
      }
    }
    
    glm.results <- get.object(
      obj={
        locations.test <- data.table::rbindlist(
          list(
            data.table(x=Q.test$data$x, y=Q.test$data$y, obs=T), 
            data.table(x=Q.test$dummy$x, y=Q.test$dummy$y, obs=F)
          )
        ) 
        if (covariates.keep > 0) {
          covariates.rasters <- names(covariates) %>%
            set_names() %>%
            purrr::map(function(x) {
              r <- terra::rast(covariates[[x]])
              set.crs(r, st_crs(4326, parameters=T)$Wkt)
              names(r) <- x
              r
            })
          purrr::walk2(names(covariates), covariates.rasters, function(n, r) {
            locations.test[, (n) := 
                             terra::extract(r, cbind(locations.test$x, 
                                                     locations.test$y))]
          })
        }
        
        glm.pred <- spatstat.model::predict.ppm(
          fit.glm, 
          locations=locations.test[, .(x,y)],
          type="trend", se=T)
        glm.ci <- spatstat.model::predict.ppm(
          fit.glm, 
          locations=locations.test[, .(x,y)],
          type="trend", interval="c")
        # Intensity
        inten <- predict(fit.glm)
        pixarea <- with(inten, xstep * ystep)
        glm.test <- cbind(
          locations.test,
          data.table(
            est=glm.pred$estimate,
            se=glm.pred$se,
            lo=glm.ci[1:length(glm.pred$estimate)],
            hi=glm.ci[(length(glm.pred$estimate) + 1):length(glm.ci)])
        )
        glm.test[, `:=` (
          p.obs = est * pixarea, 
          p.obs.lo = lo * pixarea,
          p.obs.hi = hi * pixarea
        )]
        
        all.predictions <- response.ppm(fit.glm)$window %>% as.matrix()
        glm.count.pred <- tryCatch(spatstat.model::predict.ppm(fit.glm, 
                                                               type="count", se=T),
                                   error=function(e) NULL)
        glm.pi <- tryCatch(spatstat.model::predict.ppm(fit.glm, type="count", 
                                                       se=T, interval="p"),
                           error=function(e) NULL)
        
        optimal.threshold <- optimize.threshold(roc.obj)
        list(
          test=glm.test,
          all.preds=all.predictions,
          count=glm.count.pred,
          pi=glm.pi,
          roc=roc.obj,
          thresh=optimal.threshold
        )
      },
      file.name=paste0(spec, "_", st, "_ipp_glm_mpl.rds"),
      obj.path="artifacts/test_results/ipp_glm_mpl"
    )
    cat("\tFinished IPP model for", spec, "in", st, "\n")
  }
  gc()
})
```

### Generate Accuracy Metrics using Test Data

```
get.acc <- function(results.path, thresh) {
  test <- readRDS(results.path)$test
  df <- test[, .(pred=as.factor(ifelse(p.obs > thresh, T, F)), obs=as.factor(obs))]
  
  cm <- confusionMatrix(df$pred, df$obs, positive = "TRUE", mode="everything")
}

ipp.models <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  glm.base.path <- "artifacts/models/ipp_glm_mpl"
  # Regular expression pattern to match the filename
  pattern <- paste0("^", spec, "_", st, "_([0-9]{1,2})_ipp_glm_mpl\\.rds$")
  glm.path <- list.files(path = glm.base.path, pattern = pattern, full.names = T)
  num.of.covariates <- as.numeric(regmatches(glm.path, regexpr("[0-9]{1,2}", glm.path)))
  results.path <- file.path("artifacts/test_results/ipp_glm_mpl",
                            paste0(spec, "_", st, "_ipp_glm_mpl.rds"))
  thresh <- 1-readRDS(results.path)$thresh
  cm <- get.acc(results.path, thresh)
  tibble(
    common.name=spec,
    state=st,
    covariate.count=num.of.covariates,
    optimal.threshold=thresh # ,
    # Fpred.Fref=cm$table[1,1],
    # Fpred.Tref=cm$table[1,2],
    # Tpred.Fref=cm$table[2,1],
    # Tpred.Tref=cm$table[2,2],
    # glm.path=glm.path,
    # results=results.path
  ) %>%
    cbind(as.list(c(cm$overall, cm$byClass)) %>% as_tibble())
}) %>%
  select(common.name:Accuracy,Sensitivity,Specificity,F1)
```

## Results

```
DT::datatable(
  ipp.models,
  filter='none',
  selection='none',
  rownames=F,
  options=list(
    scrollX=T,
    scrollY=T,
    paging=F,
    searching=F,
    orderMulti=T,
    info=F,
    lengthChange = F
  )) %>%
  DT::formatStyle(columns=names(ipp.models), 
                  `font-size`="13px") %>%
  DT::formatSignif(4:ncol(ipp.models), digits=2)

```

<div style="min-width: 320px; overflow-x: auto; border: 1px solid #fff;">
<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
<script src="https://code.jquery.com/jquery-3.7.0.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>

<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item html-widget-static-bound" 
     id="htmlwidget-fce9b5828451b8c22cde" style="min-width: 873px; height:auto; margin: 2px;">
  <table id="ippResultsDT" class="display dataTable no-footer table table-condensed" style="width: 873px;">
    <thead>
      <tr>
        <th>common.name</th>
        <th>state</th>
        <th>covariate.count</th>
        <th>optimal.threshold</th>
        <th>Accuracy</th>
        <th>Sensitivity</th>
        <th>Specificity</th>
        <th>F1</th>
      </tr>
    </thead>
  </table>
</div>

<script>
  // Define the JSON dataset as a JavaScript variable.
  var ippResults = [
    {"common.name": "Belted Kingfisher", "state": "CO", "covariate.count": 50, "optimal.threshold": 0.97, "Accuracy": 0.74, "Sensitivity": 0.39, "Specificity": 1.0, "F1": 0.57},
    {"common.name": "Cedar Waxwing",       "state": "CO", "covariate.count": 50, "optimal.threshold": 0.96, "Accuracy": 0.68, "Sensitivity": 0.30, "Specificity": 1.0, "F1": 0.46},
    {"common.name": "Downy Woodpecker",     "state": "CO", "covariate.count": 50, "optimal.threshold": 0.93, "Accuracy": 0.77, "Sensitivity": 0.46, "Specificity": 1.0, "F1": 0.63},
    {"common.name": "Ruddy Duck",           "state": "CO", "covariate.count": 15, "optimal.threshold": 0.68, "Accuracy": 0.57, "Sensitivity": 0.027, "Specificity": 1.0, "F1": 0.053},
    {"common.name": "Sanderling",           "state": "CO", "covariate.count": 0,  "optimal.threshold": 0.50, "Accuracy": 0.82, "Sensitivity": 0.0,   "Specificity": 1.0, "F1": ""},
    {"common.name": "Sandhill Crane",       "state": "CO", "covariate.count": 9,  "optimal.threshold": 0.44, "Accuracy": 0.53, "Sensitivity": 0.017, "Specificity": 1.0, "F1": 0.034},
    {"common.name": "Sharp-shinned Hawk",   "state": "CO", "covariate.count": 50, "optimal.threshold": 0.94, "Accuracy": 0.57, "Sensitivity": 0.12,  "Specificity": 1.0, "F1": 0.22},
    {"common.name": "Wild Turkey",          "state": "CO", "covariate.count": 50, "optimal.threshold": 0.84, "Accuracy": 0.55, "Sensitivity": 0.047, "Specificity": 1.0, "F1": 0.089},
    {"common.name": "Belted Kingfisher", "state": "NC", "covariate.count": 50, "optimal.threshold": 0.47, "Accuracy": 0.71, "Sensitivity": 0.41, "Specificity": 0.97, "F1": 0.57},
    {"common.name": "Cedar Waxwing",    "state": "NC", "covariate.count": 15, "optimal.threshold": 0.61, "Accuracy": 0.57, "Sensitivity": 0.10, "Specificity": 0.99, "F1": 0.18},
    {"common.name": "Downy Woodpecker",  "state": "NC", "covariate.count": 8,  "optimal.threshold": 0.45, "Accuracy": 0.47, "Sensitivity": 0.98, "Specificity": 0.042,"F1": 0.63},
    {"common.name": "Ruddy Duck",        "state": "NC", "covariate.count": 14, "optimal.threshold": 0.68, "Accuracy": 0.55, "Sensitivity": 0.0,  "Specificity": 1.0, "F1": ""},
    {"common.name": "Sanderling",        "state": "NC", "covariate.count": 3,  "optimal.threshold": 1.0,  "Accuracy": 0.71, "Sensitivity": 0.0,  "Specificity": 1.0, "F1": ""},
    {"common.name": "Sandhill Crane",    "state": "NC", "covariate.count": 35, "optimal.threshold": 1.0,  "Accuracy": 0.89, "Sensitivity": 0.0,  "Specificity": 1.0, "F1": ""},
    {"common.name": "Sharp-shinned Hawk","state": "NC", "covariate.count": 15, "optimal.threshold": 0.48, "Accuracy": 0.52, "Sensitivity": 0.021,"Specificity": 1.0, "F1": 0.041},
    {"common.name": "Wild Turkey",       "state": "NC", "covariate.count": 50, "optimal.threshold": 0.27, "Accuracy": 0.66, "Sensitivity": 0.39, "Specificity": 0.92, "F1": 0.53},
    {"common.name": "Belted Kingfisher", "state": "OR", "covariate.count": 4,  "optimal.threshold": 0.48, "Accuracy": 0.60, "Sensitivity": 0.051,"Specificity": 1.0, "F1": 0.097},
    {"common.name": "Cedar Waxwing",    "state": "OR", "covariate.count": 50, "optimal.threshold": 1.0,  "Accuracy": 0.81, "Sensitivity": 0.53, "Specificity": 1.0, "F1": 0.69},
    {"common.name": "Downy Woodpecker",  "state": "OR", "covariate.count": 50, "optimal.threshold": 1.0,  "Accuracy": 0.83, "Sensitivity": 0.60, "Specificity": 1.0, "F1": 0.74},
    {"common.name": "Ruddy Duck",        "state": "OR", "covariate.count": 50, "optimal.threshold": 0.49, "Accuracy": 0.65, "Sensitivity": 0.17, "Specificity": 1.0, "F1": 0.29},
    {"common.name": "Sanderling",        "state": "OR", "covariate.count": 7,  "optimal.threshold": 0.96, "Accuracy": 0.72, "Sensitivity": 0.0,  "Specificity": 1.0, "F1": ""},
    {"common.name": "Sandhill Crane",    "state": "OR", "covariate.count": 25, "optimal.threshold": 0.62, "Accuracy": 0.60, "Sensitivity": 0.032,"Specificity": 1.0, "F1": 0.062},
    {"common.name": "Sharp-shinned Hawk","state": "OR", "covariate.count": 20, "optimal.threshold": 0.44, "Accuracy": 0.66, "Sensitivity": 0.27, "Specificity": 1.0, "F1": 0.43},
    {"common.name": "Wild Turkey",       "state": "OR", "covariate.count": 9,  "optimal.threshold": 0.48, "Accuracy": 0.56, "Sensitivity": 0.035,"Specificity": 1.0, "F1": 0.066}, 
    {"common.name": "Belted Kingfisher", "state": "VT", "covariate.count": 50, "optimal.threshold": 0.68, "Accuracy": 0.54, "Sensitivity": 0.0,   "Specificity": 1.0, "F1": ""},
    {"common.name": "Cedar Waxwing",    "state": "VT", "covariate.count": 50, "optimal.threshold": 0.52, "Accuracy": 0.28, "Sensitivity": 0.028,"Specificity": 1.0, "F1": 0.055},
    {"common.name": "Downy Woodpecker",  "state": "VT", "covariate.count": 50, "optimal.threshold": 0.73, "Accuracy": 0.33, "Sensitivity": 0.029,"Specificity": 1.0, "F1": 0.056},
    {"common.name": "Ruddy Duck",        "state": "VT", "covariate.count": 4,  "optimal.threshold": 0.48, "Accuracy": 0.93, "Sensitivity": 0.0,   "Specificity": 1.0, "F1": ""},
    {"common.name": "Sanderling",        "state": "VT", "covariate.count": 0,  "optimal.threshold": 0.51, "Accuracy": 0.94, "Sensitivity": 0.0,   "Specificity": 1.0, "F1": ""},
    {"common.name": "Sandhill Crane",    "state": "VT", "covariate.count": 50, "optimal.threshold": 1.0,  "Accuracy": 0.89, "Sensitivity": 0.0,   "Specificity": 0.99,"F1": ""},
    {"common.name": "Sharp-shinned Hawk","state": "VT", "covariate.count": 50, "optimal.threshold": 0.48, "Accuracy": 0.51, "Sensitivity": 0.0048,"Specificity": 1.0, "F1": 0.0095},
    {"common.name": "Wild Turkey",       "state": "VT", "covariate.count": 50, "optimal.threshold": 0.26, "Accuracy": 0.55, "Sensitivity": 0.073,"Specificity": 0.99,"F1": 0.13}
  ];

  // Initialize the DataTable using the JSON dataset.
  jQuery(document).ready(function($){
    $('#ippResultsDT').DataTable({
      data: ippResults,
      columns: [
        { data: function(row) { return row["common.name"]; } },
        { data: "state" },
        { data: function(row) { return row["covariate.count"]; } },
        { data: function(row) { return row["optimal.threshold"]; } },
        { data: "Accuracy" },
        { data: "Sensitivity" },
        { data: "Specificity" },
        { data: "F1" }
      ],
      scrollX: false,
      searching: false
    });
  });
</script>


</div>

### Overall Metric Summaries

```
select(ipp.models, Accuracy:F1) %>% summary()
```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
##     Accuracy       Sensitivity       Specificity            F1          
##  Min.   :0.2766   Min.   :0.00000   Min.   :0.04218   Min.   :0.009479  
##  1st Qu.:0.5487   1st Qu.:0.00000   1st Qu.:0.99603   1st Qu.:0.058968  
##  Median :0.6258   Median :0.03321   Median :0.99932   Median :0.183740  
##  Mean   :0.6462   Mean   :0.15896   Mean   :0.96510   Mean   :0.288013  
##  3rd Qu.:0.7455   3rd Qu.:0.27887   3rd Qu.:1.00000   3rd Qu.:0.549047  
##  Max.   :0.9412   Max.   :0.98264   Max.   :1.00000   Max.   :0.744135  
##               
</code></pre></div><br>

```
# Reshape the data to long format
long.data <- ipp.models %>%
  select(common.name, Accuracy:F1) %>%
  gather(key = "Metric", value = "Value", -common.name)

# Overall box-plot representation
ggplot(long.data, aes(x = Metric, y = Value)) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, color = "red") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "lightgray"),
        plot.background = element_rect(fill = "white"),
        legend.position = "none") +
  coord_flip()

```

<img src="{{ site.baseurl }}/assets/plots/sdm-5-plt-1.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 650px; height: auto;">

### Summaries by Species

```
# Box plot for each species
ggplot(long.data, aes(x = Metric, y = Value)) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, color = "red") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "lightgray"),
        plot.background = element_rect(fill = "white"),
        legend.position = "none") +
  coord_flip() +
  facet_wrap(~ common.name, scales = "free_y", ncol = 2)

```

<img src="{{ site.baseurl }}/assets/plots/sdm-5-plt-2.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 650px; height: auto;">

### Summaries by State

```
# Reshape the data to long format
long.data.state <- ipp.models %>%
  select(state, Accuracy:F1) %>%
  gather(key = "Metric", value = "Value", -state)

# Box plot for each species
ggplot(long.data.state, aes(x = Metric, y = Value)) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, color = "red") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "lightgray"),
        plot.background = element_rect(fill = "white"),
        legend.position = "none") +
  coord_flip() +
  facet_wrap(~ state, scales = "free_y", ncol = 2)

```

<img src="{{ site.baseurl }}/assets/plots/sdm-5-plt-3.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 650px; height: auto;">

### Discussion

To put it frankly, the apparent results of the models given the test data
are not good. On average, accuracy is ~65% and the F1 score is ~0.29. 
Perhaps most noticeable in the results is that the Sensitivity and Specificity are consistently very polarized, averaging ~0.97 and ~0.16 respectively. This means that the models are consistently worse at
identifying the positive class (bird observations) correctly.

Given this information when typically evaluating a predictive model,
sensitivity and specificity metrics like these might suggest that 
the training data is imbalanced. I.e., there are many more examples 
of one class than the other. Models trained on imbalanced data can 
often show high sensitivity but low specificity because they learn 
to predict the majority class more frequently. However, given that 
the data used for the negative values in these models is pseudo-absence 
data that was specifically sampled in balanced sample sizes, this 
cannot be the case. But the problem still most-likely lies with
the pseudo-absence data. If a model is not correctly distinguishing
between the observation/absence points, it might mean that the 
pseudo-absence points do not accurately represent the 
regions where a particular bird is less likely to be observed.

As a next step, I will be re-evaluating the pseudo-absence selection
process and making some improvements. From there, I will re-fit
the IPP models, and proceed with the study using a (hopefully) more
accurate representation of the observation/absence data.
