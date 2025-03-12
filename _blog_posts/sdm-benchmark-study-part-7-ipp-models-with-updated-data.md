---
layout: post
title: "SDM Benchmark Study Part 7: Fitting and Testing Inhomogeneous Poisson Process Models with Updated Data"
permalink: /blog_posts/sdm-benchmark-study-part-7-ipp-models-with-updated-data
---


## Overview

In [part 5](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-5-fitting-and-testing-ipp-models) 
of the study, a first attempt was made to fit IPP models to the data. Although the models are meant to serve as a 
baseline, the models had very poor performance. In an attempt to improve the results, several changes were made
to the pseudo-absence selection process (see 
[part 6](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-6-resampling-pseudoabsence-points)
of the study). In this part of the study, the updated data will be used to re-fit the
models. In addition, a few changes have been made to the model-fitting process. Those
changes are highlighted in the outline below.

As in some of the prior parts, there will not be a detailed description for section
and block of code. Rather, the outline will lay things out at a high-level, and at the
end of this part of the study there will be a summarization and discussion of the
results.

### IPP Model Process Outline (With Updates)

*Steps refined from [part 5](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-5-fitting-and-testing-ipp-models) are highlighted in bold.*

- Load libraries
- Read in observation and raster data from [part 2](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-2-exploratory-analysis) and [part 3](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-3-more-preprocessing-and-eda) of the study
- Establish utility functions and variables
- **Load variable importance (obtained from fitted LASSO models in [part 6](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-6-resampling-pseudoabsence-points) of the study)**
- Iterate over each state/species pair
  * Define the base number of maximum covariates and covariate interactions to keep (50)
  * **Load pre-processed train/test observation data (`spatstat.geom::quadscheme` objects cached in [part 6](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-6-resampling-pseudoabsence-points) of the study)**
  * Load the corresponding pixel images (`spatstat.geom::im` objects, cached in 
    [part 3](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-3-more-preprocessing-and-eda)) for each of the covariates that align with those
    in the variable importance data, up to the maximum allowed number
  * Using the training data and pixel images, fit an IPP (GLM) Model using the Method of
    Maximum Pseudo-Likelihood (`spatstat.model::ppm`) 
  * Check for issues: Model does not converge, Prediction error, **or Specificity and Sensitivity are a 0/1 pair (i.e., Specificity is 0 and Sensitivity is 1, or vice versa) when probabilities are generated from the training data and the optimal threshold (see below)**
    - If there is an issue, reduce the number of maximum allowed covariates
      and re-fit the model
    - Repeat this process until there are either no problems or else there
      are no more available covariates
    - **If there are still problems, restart the process using the same covariates but with no interactions. If the issues still persist, use a model with no covariates**
  * Generate predictions with standard errors and a confidence interval for the 
    test data
  * Using the estimates and C.I. values, estimate probabilities 
  * Generate predictions for the entire region
  * Estimate prediction intervals with standard errors
  * **Find the optimal threshold by maximizing the Specificity and Sensitivity, and ensuring neither is equal to zero**
  * Generate model accuracy metrics using the test predictions/actuals
- Visualize results (using the test data), and **compare with the metrics of the prior models that were fit in [part 5](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-5-fitting-and-testing-ipp-models)**

## Setup

```
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

## IPP Models

### Load Variable Importance

```
# Load variable importance from fitted LASSO models
lasso.model.path="artifacts/models/lasso_3_fs"

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

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script src="{{ site.baseurl }}/assets/data/sdm-6-vi.js"></script>
<script data-pagedtable-source data-global="sdm6VI" type="application/json"></script>
</div>

```
get.f1 <- function(locations.df, thresh=0.5, method="F1") {
  metrics <- locations.df[, .(pred=factor(ifelse(p.obs > thresh, T, F), 
                               levels=c("FALSE", "TRUE")), 
          obs=factor(obs, levels=c("FALSE", "TRUE")))] %>% 
    confusionMatrix(reference=.$pred, data=.$obs, 
                    positive = "TRUE", 
                    mode="everything") %>%
    as.list() %>%
    .[["byClass"]] %>%
    as.list()
  if (method == "F1") {
    return(metrics$F1)
  } else if (method == "SS") {
    # Penalize heavily if either is near zero or NA
    if (is.na(metrics$Sensitivity) || is.na(metrics$Specificity)) {
      return(-1e6) 
    } else if (metrics$Sensitivity < .05 || metrics$Specificity < .05) {
      return(-1e6) 
    }
    
    # Combine Specificity and Sensitivity (e.g., geometric mean)
    combined.metric = sqrt(metrics$Sensitivity * metrics$Specificity)
    return(combined.metric)
  }
}

optimize.f1 <- function(locations.df) {
  objective.fn <- function(thresh) {
    # Negative because we want to maximize
    -get.f1(locations.df=locations.df, thresh=thresh, method="SS")  
  }
  opt.result <- optimize(objective.fn, lower=0.05, upper=0.95)
  return(opt.result$minimum)  # Return the optimal threshold
}

get.acc <- function(test, thresh) {
  df <- test[, .(pred=as.factor(ifelse(p.obs > thresh, T, F)), obs=as.factor(obs))]
  cm <- confusionMatrix(df$pred, df$obs, positive = "TRUE", mode="everything")
}

# 9, 26, 27, 32
purrr::walk(1:nrow(spec.state), function(i) {
  covariates.keep <- 50
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  glm.results.path <- file.path("artifacts/test_results/ipp_glm_mpl_2",
                                paste0(spec, "_", st, "_ipp_glm_mpl_2.rds"))
  if (!file.exists(glm.results.path)) {
    cat("Fitting IPP model for", spec, "in", st, "\n")
    
    cat("\tGetting pre-processed `spatstat.geom::ppp` object (train)...\n")
    Q <- readRDS(file.path("artifacts", "train_spatstat_Q_2",
                           paste0(st, "_", spec, "_Q.rds")))
    cat("\tGetting `spatstat.geom::ppp` object (test)...\n")
    Q.test <- readRDS(file.path("artifacts", "test_spatstat_Q_2",
                                paste0(st, "_", spec, "_Q.rds")))
    
    converged <- F
    no.pred.error <- F
    # roc.obj <- NULL
    spec.sens.check <- F
    alternate.formula <- F
    while (!converged | !no.pred.error | !spec.sens.check) {
      # Set path
      glm.path <- file.path("artifacts/models/ipp_glm_mpl_2", 
                            paste0(spec, "_", st, "_", covariates.keep,
                                   "_ipp_glm_mpl_2.rds"))
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
        
        if (nrow(fs.df) > 0) {
          covariates <- c(fs.df$var1, fs.df$var2) %>% 
            unique() %>% 
            sort()
          
          # Load/compute filtered & pre-processed rasters
          covariates <- covariates %>%
            set_names() %>%
            purrr::map(function(.x) {
              file <- file.path("artifacts/spatstat_imgs", 
                                paste0(st, "_", .x, ".rds"))
              readRDS(file)
            })
          
          # Create formula
          if (alternate.formula) {
            if (length(covariates) < covariates.keep) {
              covariates.keep <- length(covariates) # in case there are fewer available
            } else if (length(covariates) > covariates.keep) {
              covariates <- covariates[1:covariates.keep]
            }
            
            # No interactions
            .f <- paste(names(covariates), collapse = " + ") %>% 
              paste("~", .) %>% 
              as.formula()
          } else {
            if (nrow(fs.df) < covariates.keep) {
              covariates.keep <- nrow(fs.df) # in case there are fewer available
            }
            # With interactions
            .f <- paste(fs.df$variable, collapse=" + ") %>% 
              paste("~", .) %>% 
              as.formula()
          }
          # Make sure glm.path is right
          glm.path <- file.path("artifacts/models/ipp_glm_mpl_2", 
                            paste0(spec, "_", st, "_", covariates.keep,
                                   "_ipp_glm_mpl_2.rds"))
          
        } else {
          cat("\tThere are no specified covariates for", spec, st, "\n")
        }
        
        # Fit the IPP model, using the Method of Maximum PseudoLikelihood (MPL)
        
        # GLM Model
        cat("\tFitting GLM Model using the Method of Maximum",
            "PseudoLikelihood...\n")
        fit.glm <- tryCatch({
          ppm(Q=Q, trend=.f, covariates=covariates, 
              rbord=.05, method="mpl", emend=T) %>%
            get.object(
              obj=.,
              file.name=paste0(spec, "_", st, "_", covariates.keep,
                               "_ipp_glm_mpl_2.rds"), 
              obj.path="artifacts/models/ipp_glm_mpl_2")}, 
          error=function(e) NULL)
        if (is.null(fit.glm)) {
          if (covariates.keep > 15) {
            covariates.keep <- covariates.keep - 5
          } else {
            covariates.keep <- covariates.keep - 1
          }
          if (covariates.keep < 3 & !alternate.formula) {
            cat("\tNo trend identified. Removing interactions...\n")
            alternate.formula <- T
            covariates.keep <- 50
          }  else if (covariates.keep < 1 & alternate.formula) {
            cat("\tNo trend identified. Removing covariates...\n")
          } else if (covariates.keep < 0 & alternate.formula) {
            stop("\tUnable to successfully fit a model given the data.\n")
          } 
          next
        }
        converged <- fit.glm$internal$glmfit$converged
      } else {
        cat("\tFitting GLM Model using the Method of Maximum",
            "PseudoLikelihood with no trend...\n")
        fit.glm <- ppm(Q=Q, rbord=.05, method="mpl", emend=T) %>%
          get.object(
            obj=.,
            file.name=paste0(spec, "_", st, "_", covariates.keep,
                             "_ipp_glm_mpl_2.rds"), 
            obj.path="artifacts/models/ipp_glm_mpl_2")
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
      if (!converged | !no.pred.error) {
        cat("\tThe model converged or there was an error",
                "for covariates.keep ==", covariates.keep, "\n")
        
        file.remove(glm.path)
        if (covariates.keep > 15) {
          covariates.keep <- covariates.keep - 5
        } else {
          covariates.keep <- covariates.keep - 1
        }
        if (covariates.keep < 3 & !alternate.formula) {
          cat("\tNo trend identified. Removing interactions...\n")
          alternate.formula <- T
          covariates.keep <- 50
        }  else if (covariates.keep < 1 & alternate.formula) {
          cat("\tNo trend identified. Removing covariates...\n")
        } else if (covariates.keep < 0 & alternate.formula) {
          stop("\tUnable to successfully fit a model given the data.\n")
        } 
        next
      } else {
        locations.train <- data.table::rbindlist(
          list(
            data.table(x=Q$data$x, y=Q$data$y, obs=T), 
            data.table(x=Q$dummy$x, y=Q$dummy$y, obs=F)
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
            locations.train[, (n) := 
                             terra::extract(r, cbind(locations.train$x, 
                                                     locations.train$y))]
          })
        }
        glm.pred <- spatstat.model::predict.ppm(
          fit.glm, 
          locations=locations.train[, .(x,y)],
          type="trend", se=T)
        # Intensity
        inten <- predict(fit.glm)
        pixarea <- with(inten, xstep * ystep)
        glm.train <- cbind(
          locations.train,
          data.table(
            est=glm.pred$estimate,
            se=glm.pred$se,
            p.obs = glm.pred$estimate*pixarea)
        )
        optimal.threshold <- optimize.f1(glm.train)
        cm <- get.acc(glm.train, optimal.threshold)
        acc <- tibble(
          common.name=spec,
          state=st,
          covariate.count=covariates.keep,
          optimal.threshold=optimal.threshold 
        ) %>%
          cbind(as.list(c(cm$overall, cm$byClass)) %>% 
                  as_tibble()) %>%
        select(common.name:Accuracy, Sensitivity, Specificity, F1)
        cat("Train Results:\n Species:", spec, "\n",
            "State:", st, "\n",
            "Covariates:", covariates.keep, "\n",
            "Optimal Threshold:", optimal.threshold, "\n",
            "Accuracy:", acc$Accuracy, "\n",
            "F1:", acc$F1, "\n",
            "Sensitivity (TP Rate):", acc$Sensitivity, "\n",
            "Specificity (TN Rate):", acc$Specificity, "\n")
        if ((acc$Sensitivity == 1 & acc$Specificity == 0) | 
             (acc$Sensitivity == 0 & acc$Specificity == 1)) {
          cat("\tThe sensitivity/specificity is a 0/1 pair",
              "for covariates.keep ==", covariates.keep, "\n")
          spec.sens.check <- F
          file.remove(glm.path)
          if (covariates.keep > 15) {
            covariates.keep <- covariates.keep - 5
          } else {
            covariates.keep <- covariates.keep - 1
          }
          if (covariates.keep < 3 & !alternate.formula) {
            cat("\tNo trend identified. Removing interactions...\n")
            alternate.formula <- T
            covariates.keep <- 50
          }  else if (covariates.keep == 0 & alternate.formula) {
            cat("\tNo trend identified. Removing covariates...\n")
          } else if (covariates.keep < 0 & alternate.formula) {
            cat("\tNo other meaningful options to check, 
                   assuming the model has no trend;\n")
            spec.sens.check <- T
          } 
          next
        } else {
          spec.sens.check <- T
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
        glm.test <- cbind(
          locations.test,
          data.table(
            est=glm.pred$estimate,
            se=glm.pred$se,
            lo=glm.ci[1:length(glm.pred$estimate)],
            hi=glm.ci[(length(glm.pred$estimate) + 1):length(glm.ci)])
        )
        
        glm.test[, `:=` (
          p.obs = est * pixarea, # probability > 0
          p.obs.lo = est * lo,
          p.obs.hi = est * hi
        )]
        cm <- get.acc(glm.test, optimal.threshold)
        test.acc <- tibble(
          common.name=spec,
          state=st,
          covariate.count=covariates.keep,
          optimal.threshold=optimal.threshold 
        ) %>%
          cbind(as.list(c(cm$overall, cm$byClass)) %>% 
                  as_tibble()) %>%
          select(common.name:Accuracy, Sensitivity, Specificity, F1)
        all.predictions <- response.ppm(fit.glm)$window %>% as.matrix()
        glm.count.pred <- tryCatch(spatstat.model::predict.ppm(fit.glm, 
                                                               type="count", se=T),
                           error=function(e) NULL)
        glm.pi <- tryCatch(spatstat.model::predict.ppm(fit.glm, type="count", 
                                                       se=T, interval="p"),
                           error=function(e) NULL)
        list(
          test=glm.test,
          train=glm.train,
          all.preds=all.predictions,
          count=glm.count.pred,
          pi=glm.pi,
          thresh=optimal.threshold,
          train.accuracy=acc,
          test.accuracy=test.acc
        )
      },
      file.name=paste0(spec, "_", st, "_ipp_glm_mpl_2.rds"),
      obj.path="artifacts/test_results/ipp_glm_mpl_2"
    )
    cat("\tFinished IPP model for", spec, "in", st, "\n")
  }
  gc()
})

```


## Results

```
ipp.models <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  results.path <- file.path("artifacts/test_results/ipp_glm_mpl_2",
                        paste0(spec, "_", st, "_ipp_glm_mpl_2.rds"))
  readRDS(results.path)$test.accuracy
})

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
       id="htmlwidget-fce9b5828451b8c22cde" style="min-width: 865px; height:auto; margin: 2px;">
    <table id="ippResultsDT" class="display dataTable no-footer table table-condensed" style="width: 865px;">
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
      {
        "common.name": "Belted Kingfisher",
        "state": "CO",
        "covariate.count": 50,
        "optimal.threshold": 0.0854754125117056,
        "Accuracy": 0.979985174203113,
        "Sensitivity": 0.968124536693847,
        "Specificity": 0.99184581171238,
        "F1": 0.979744936234059
      },
      {
        "common.name": "Cedar Waxwing",
        "state": "CO",
        "covariate.count": 50,
        "optimal.threshold": 0.0980961138985813,
        "Accuracy": 0.952755905511811,
        "Sensitivity": 0.923228346456693,
        "Specificity": 0.982283464566929,
        "F1": 0.95131845841785
      },
      {
        "common.name": "Downy Woodpecker",
        "state": "CO",
        "covariate.count": 50,
        "optimal.threshold": 0.166043311031947,
        "Accuracy": 0.974284436493739,
        "Sensitivity": 0.954850245864998,
        "Specificity": 0.993736017897092,
        "F1": 0.973786186459995
      },
      {
        "common.name": "Ruddy Duck",
        "state": "CO",
        "covariate.count": 50,
        "optimal.threshold": 0.0536288446520682,
        "Accuracy": 0.969939879759519,
        "Sensitivity": 0.949899799599198,
        "Specificity": 0.98997995991984,
        "F1": 0.969325153374233
      },
      {
        "common.name": "Sanderling",
        "state": "CO",
        "covariate.count": -1,
        "optimal.threshold": 0.949944437279511,
        "Accuracy": 0.869565217391304,
        "Sensitivity": 0,
        "Specificity": 1,
        "F1": null
      },
      {
        "common.name": "Sandhill Crane",
        "state": "CO",
        "covariate.count": 50,
        "optimal.threshold": 0.0516029012081977,
        "Accuracy": 0.959550561797753,
        "Sensitivity": 0.925842696629213,
        "Specificity": 0.993258426966292,
        "F1": 0.958139534883721
      },
      {
        "common.name": "Sharp-shinned Hawk",
        "state": "CO",
        "covariate.count": 30,
        "optimal.threshold": 0.131607996769709,
        "Accuracy": 0.88212927756654,
        "Sensitivity": 0.776595744680851,
        "Specificity": 0.987823439878234,
        "F1": 0.868309260832625
      },
      {
        "common.name": "Wild Turkey",
        "state": "CO",
        "covariate.count": 50,
        "optimal.threshold": 0.0559273100829301,
        "Accuracy": 0.963446475195822,
        "Sensitivity": 0.975195822454308,
        "Specificity": 0.951697127937337,
        "F1": 0.963870967741935
      },
      {
        "common.name": "Belted Kingfisher",
        "state": "NC",
        "covariate.count": 17,
        "optimal.threshold": 0.395158765057634,
        "Accuracy": 0.758103727714749,
        "Sensitivity": 0.548281505728314,
        "Specificity": 0.963884430176565,
        "F1": 0.691791430046464
      },
      {
        "common.name": "Cedar Waxwing",
        "state": "NC",
        "covariate.count": 13,
        "optimal.threshold": 0.30712336120542,
        "Accuracy": 0.882020202020202,
        "Sensitivity": 0.792038992688871,
        "Specificity": 0.971061093247588,
        "F1": 0.869759143621766
      },
      {
        "common.name": "Downy Woodpecker",
        "state": "NC",
        "covariate.count": 50,
        "optimal.threshold": 0.695478595267957,
        "Accuracy": 0.864607345935147,
        "Sensitivity": 0.794761171032357,
        "Specificity": 0.93408951563458,
        "F1": 0.85411491968869
      },
      {
        "common.name": "Ruddy Duck",
        "state": "NC",
        "covariate.count": 30,
        "optimal.threshold": 0.0656676044654924,
        "Accuracy": 0.861370716510903,
        "Sensitivity": 0.817610062893082,
        "Specificity": 0.904320987654321,
        "F1": 0.853858784893268
      },
      {
        "common.name": "Sanderling",
        "state": "NC",
        "covariate.count": 30,
        "optimal.threshold": 0.949944437279511,
        "Accuracy": 0.901140684410646,
        "Sensitivity": 0.609375,
        "Specificity": 0.994974874371859,
        "F1": 0.75
      },
      {
        "common.name": "Sandhill Crane",
        "state": "NC",
        "covariate.count": 35,
        "optimal.threshold": 0.949944437279511,
        "Accuracy": 0.865217391304348,
        "Sensitivity": 0,
        "Specificity": 0.995,
        "F1": null
      },
      {
        "common.name": "Sharp-shinned Hawk",
        "state": "NC",
        "covariate.count": 50,
        "optimal.threshold": 0.0972626576317735,
        "Accuracy": 0.873103448275862,
        "Sensitivity": 0.797222222222222,
        "Specificity": 0.947945205479452,
        "F1": 0.861861861861862
      },
      {
        "common.name": "Wild Turkey",
        "state": "NC",
        "covariate.count": 50,
        "optimal.threshold": 0.156677070020394,
        "Accuracy": 0.865329512893983,
        "Sensitivity": 0.771223021582734,
        "Specificity": 0.958630527817404,
        "F1": 0.850793650793651
      },
      {
        "common.name": "Belted Kingfisher",
        "state": "OR",
        "covariate.count": 50,
        "optimal.threshold": 0.0829636963690363,
        "Accuracy": 0.98136459272618,
        "Sensitivity": 0.975487115021999,
        "Specificity": 0.986751152073733,
        "F1": 0.980416929879975
      },
      {
        "common.name": "Cedar Waxwing",
        "state": "OR",
        "covariate.count": 13,
        "optimal.threshold": 0.382898196341097,
        "Accuracy": 0.947765762089369,
        "Sensitivity": 0.904442581726739,
        "Specificity": 0.988866799204771,
        "F1": 0.944006999125109
      },
      {
        "common.name": "Downy Woodpecker",
        "state": "OR",
        "covariate.count": 50,
        "optimal.threshold": 0.107111419353424,
        "Accuracy": 0.987367154601965,
        "Sensitivity": 0.979046836483155,
        "Specificity": 0.99529964747356,
        "F1": 0.986953820666805
      },
      {
        "common.name": "Ruddy Duck",
        "state": "OR",
        "covariate.count": 15,
        "optimal.threshold": 0.0756447011048129,
        "Accuracy": 0.957559681697613,
        "Sensitivity": 0.917543859649123,
        "Specificity": 0.998217468805704,
        "F1": 0.956124314442413
      },
      {
        "common.name": "Sanderling",
        "state": "OR",
        "covariate.count": 13,
        "optimal.threshold": 0.131230830996797,
        "Accuracy": 0.872427983539095,
        "Sensitivity": 0.369565217391304,
        "Specificity": 0.989847715736041,
        "F1": 0.523076923076923
      },
      {
        "common.name": "Sandhill Crane",
        "state": "OR",
        "covariate.count": 15,
        "optimal.threshold": 0.139665443000992,
        "Accuracy": 0.972305389221557,
        "Sensitivity": 0.968926553672316,
        "Specificity": 0.976114649681529,
        "F1": 0.973740241305891
      },
      {
        "common.name": "Sharp-shinned Hawk",
        "state": "OR",
        "covariate.count": 50,
        "optimal.threshold": 0.14203152064344,
        "Accuracy": 0.859445519019987,
        "Sensitivity": 0.806324110671937,
        "Specificity": 0.910353535353535,
        "F1": 0.848821081830791
      },
      {
        "common.name": "Wild Turkey",
        "state": "OR",
        "covariate.count": 50,
        "optimal.threshold": 0.130835019727487,
        "Accuracy": 0.91220556745182,
        "Sensitivity": 0.837762237762238,
        "Specificity": 0.989795918367347,
        "F1": 0.906888720666162
      },
      {
        "common.name": "Belted Kingfisher",
        "state": "VT",
        "covariate.count": 50,
        "optimal.threshold": 0.148488691116605,
        "Accuracy": 0.803156146179402,
        "Sensitivity": 0.698492462311558,
        "Specificity": 0.906095551894563,
        "F1": 0.778711484593838
      },
      {
        "common.name": "Cedar Waxwing",
        "state": "VT",
        "covariate.count": 50,
        "optimal.threshold": 0.105635695388682,
        "Accuracy": 0.797333333333333,
        "Sensitivity": 0.981212638770282,
        "Specificity": 0.142857142857143,
        "F1": 0.883166794773251
      },
      {
        "common.name": "Downy Woodpecker",
        "state": "VT",
        "covariate.count": 50,
        "optimal.threshold": 0.0701497281086996,
        "Accuracy": 0.752427184466019,
        "Sensitivity": 0.994181818181818,
        "Specificity": 0.058455114822547,
        "F1": 0.856248042593173
      },
      {
        "common.name": "Ruddy Duck",
        "state": "VT",
        "covariate.count": 11,
        "optimal.threshold": 0.262539061372081,
        "Accuracy": 0.943396226415094,
        "Sensitivity": 0,
        "Specificity": 1,
        "F1": null
      },
      {
        "common.name": "Sanderling",
        "state": "VT",
        "covariate.count": 6,
        "optimal.threshold": 0.262539061372081,
        "Accuracy": 0.947619047619048,
        "Sensitivity": 0,
        "Specificity": 0.995,
        "F1": null
      },
      {
        "common.name": "Sandhill Crane",
        "state": "VT",
        "covariate.count": 35,
        "optimal.threshold": 0.0511460594100544,
        "Accuracy": 0.91324200913242,
        "Sensitivity": 0.105263157894737,
        "Specificity": 0.99,
        "F1": 0.173913043478261
      },
      {
        "common.name": "Sharp-shinned Hawk",
        "state": "VT",
        "covariate.count": 50,
        "optimal.threshold": 0.199957992503922,
        "Accuracy": 0.870748299319728,
        "Sensitivity": 0.771689497716895,
        "Specificity": 0.968468468468468,
        "F1": 0.855696202531646
      },
      {
        "common.name": "Wild Turkey",
        "state": "VT",
        "covariate.count": 50,
        "optimal.threshold": 0.199957992503922,
        "Accuracy": 0.709006928406467,
        "Sensitivity": 0.460587326120556,
        "Specificity": 0.955521472392638,
        "F1": 0.611909650924025
      }
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
          { data: function(row) { return row["F1"] === null ? "" : row["F1"]; } }
        ],
        scrollX: true,
        scrollY: true,
        paging: false,
        searching: false
      });
    });
  </script>
</div>


### Load Original Model Results (from Part 5)

```
ipp.models.old <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  glm.base.path <- "artifacts/models/ipp_glm_mpl"
  # Regular expression pattern to match the filename
  pattern <- paste0("^", spec, "_", st, "_([0-9]{1,2})_ipp_glm_mpl\\.rds$")
  glm.path <- list.files(path = glm.base.path, pattern = pattern, full.names = T)
  num.of.covariates <- as.numeric(regmatches(glm.path, regexpr("[0-9]{1,2}", glm.path)))
  results.path <- file.path("artifacts/test_results/ipp_glm_mpl",
                            paste0(spec, "_", st, "_ipp_glm_mpl.rds"))
  res <- readRDS(results.path)
  thresh <- 1-res$thresh
  cm <- get.acc(res$test, thresh)
  tibble(
    common.name=spec,
    state=st,
    covariate.count=num.of.covariates,
    optimal.threshold=thresh 
  ) %>%
    cbind(as.list(c(cm$overall, cm$byClass)) %>% as_tibble())
}) %>%
  select(common.name:Accuracy,Sensitivity,Specificity,F1)

```

### Overall Metric Summaries

```
# Reshape the data to long format
long.data <- ipp.models %>%
  select(common.name, Accuracy:F1) %>%
  gather(key = "Metric", value = "Value", -common.name)

long.data.old <- ipp.models.old %>%
  select(common.name, Accuracy:F1) %>%
  gather(key = "Metric", value = "Value", -common.name)

# Combine the old and new data frames, adding an identifier column
long.data$new.old <- "Updated Models"
long.data.old$new.old <- "Original Models"
combined.long <- rbind(long.data, long.data.old)

ggplot(combined.long, 
         aes(x = Metric, y = Value, fill = new.old)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(aes(group = new.old), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75),) +
    scale_fill_manual(values = c("#FF9999", "#00659F"),
                      name = "Model Status",
                      labels = c("Original Models", "Updated Models")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "lightgray"),
          plot.background = element_rect(fill = "white")) + 
    coord_flip() +
    facet_wrap(~Metric, ncol=1, scales="free_y")
```

<img src="{{ site.baseurl }}/assets/plots/sdm-7-1.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 500px; height: auto;">

#### Original Models

```
summary(ipp.models.old)
```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
##  common.name           state           covariate.count optimal.threshold
##  Length:32          Length:32          Min.   : 0.00   Min.   :0.2649   
##  Class :character   Class :character   1st Qu.: 9.00   1st Qu.:0.4804   
##  Mode  :character   Mode  :character   Median :42.50   Median :0.6114   
##                                        Mean   :30.72   Mean   :0.6669   
##                                        3rd Qu.:50.00   3rd Qu.:0.9460   
##                                        Max.   :50.00   Max.   :1.0000   
##                                                                         
##     Accuracy       Sensitivity       Specificity            F1          
##  Min.   :0.2766   Min.   :0.00000   Min.   :0.04218   Min.   :0.009479  
##  1st Qu.:0.5487   1st Qu.:0.00000   1st Qu.:0.99603   1st Qu.:0.058968  
##  Median :0.6258   Median :0.03321   Median :0.99932   Median :0.183740  
##  Mean   :0.6462   Mean   :0.15896   Mean   :0.96510   Mean   :0.288013  
##  3rd Qu.:0.7455   3rd Qu.:0.27887   3rd Qu.:1.00000   3rd Qu.:0.549047  
##  Max.   :0.9412   Max.   :0.98264   Max.   :1.00000   Max.   :0.744135  
##                                                       NA's   :9
</code></pre></div><br>

#### Updated Models 

```
summary(ipp.models)
```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
##  common.name           state           covariate.count optimal.threshold
##  Length:32          Length:32          Min.   :-1.00   Min.   :0.05115  
##  Class :character   Class :character   1st Qu.:14.50   1st Qu.:0.08485  
##  Mode  :character   Mode  :character   Median :42.50   Median :0.13564  
##                                        Mean   :33.81   Mean   :0.26414  
##                                        3rd Qu.:50.00   3rd Qu.:0.27369  
##                                        Max.   :50.00   Max.   :0.94994  
##                                                                         
##     Accuracy       Sensitivity      Specificity            F1        
##  Min.   :0.7090   Min.   :0.0000   Min.   :0.05846   Min.   :0.1739  
##  1st Qu.:0.8651   1st Qu.:0.5941   1st Qu.:0.95457   1st Qu.:0.8503  
##  Median :0.8916   Median :0.8018   Median :0.98729   Median :0.8690  
##  Mean   :0.8953   Mean   :0.6992   Mean   :0.91913   Mean   :0.8456  
##  3rd Qu.:0.9581   3rd Qu.:0.9511   3rd Qu.:0.99338   3rd Qu.:0.9596  
##  Max.   :0.9874   Max.   :0.9942   Max.   :1.00000   Max.   :0.9870  
##                                                      NA's   :4
</code></pre></div><br>

### Summaries by Species

```
ggplot(combined.long, 
         aes(x = Metric, y = Value, fill = new.old)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(aes(group = new.old), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75),) +
    scale_fill_manual(values = c("#FF9999", "#00659F"),
                      name = "Model Status",
                      labels = c("Original Models", "Updated Models")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "lightgray"),
          plot.background = element_rect(fill = "white")) + 
    coord_flip() +
    facet_wrap(~common.name + Metric, ncol=2, scales="free_y")

```

<img src="{{ site.baseurl }}/assets/plots/sdm-7-2.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 800px; height: auto;">

### Summaries by State

```
# Reshape the data to long format
long.data <- ipp.models %>%
  select(state, Accuracy:F1) %>%
  gather(key = "Metric", value = "Value", -state)

long.data.old <- ipp.models.old %>%
  select(state, Accuracy:F1) %>%
  gather(key = "Metric", value = "Value", -state)

# Combine the old and new data frames, adding an identifier column
long.data$new.old <- "Updated Models"
long.data.old$new.old <- "Original Models"
combined.long <- rbind(long.data, long.data.old)

# Box plots for each state
ggplot(combined.long, 
         aes(x = Metric, y = Value, fill = new.old)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(aes(group = new.old), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75),) +
    scale_fill_manual(values = c("#FF9999", "#00659F"),
                      name = "Model Status",
                      labels = c("Original Models", "Updated Models")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "lightgray"),
          plot.background = element_rect(fill = "white")) + 
    coord_flip() +
    facet_wrap(~state + Metric, ncol=2, scales="free_y")

```
<img src="{{ site.baseurl }}/assets/plots/sdm-7-3.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 800px; height: auto;">

### Discussion

The revised models, following adjustments from [part 5](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-5-fitting-and-testing-ipp-models), exhibit improved performance with average accuracy elevating to ~90% and F1 score to ~0.85. Unlike before, Sensitivity and Specificity metrics are no longer polarized, indicating better identification of the positive class. This progress suggests that modifications to the pseudo-absence selection and model-fitting processes have positively impacted the model's predictive accuracy. Further refinements, especially concerning pseudo-absence data, may continue to enhance the model's performance in distinguishing between observation/absence points, fostering a more accurate representation for future analysis.