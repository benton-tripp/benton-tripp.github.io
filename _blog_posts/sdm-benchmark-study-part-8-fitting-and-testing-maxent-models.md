---
layout: post
title: "SDM Benchmark Study Part 8: Fitting and Testing MaxEnt Models"
permalink: /blog_posts/sdm-benchmark-study-part-8-fitting-and-testing-maxent-models
---


## Overview

This part of the project fits additional "baseline" models using the MaxEnt
algorithm, as described in [Part 4](https://benton-tripp.github.io/posts/2023-10-05-sdm-benchmark-study-part-4-baseline-species-distribution-models.html) 
of the study. Because an effective process for fitting the models was set up in parts 5-7, a detailed outline is mostly unnecessary 
at this point. Similarly, there will not be a detailed description for each section 
and block of code. The primary write-up will be in the summarization and discussion of
the results at the end of this post.

## Setup

Moving forward, many of the functions/code that have been used in prior parts of the 
study will not be re-written out. Rather, they will be sourced. The scripts being
sourced from external sources outside of this notebook can be found in the 
[Github repository](https://github.com/benton-tripp/presence-only-sdm/tree/86674646f9ff4411984c576228d3ddc2f1b5d3cd/R) used for this project.

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
library(predicts)

# Set seed 
set.seed(19)

# Load pre-processed data
source("R/load_preprocessed_data.R")
# Load other utility function
source("R/get_objects.R")
# Load Variable Importance
source("R/load_variable_importance.R")
# Load model accuracy function
source("R/model_accuracy.R")

```

## MaxEnt Models

```
purrr::walk(1:nrow(spec.state), function(i) {
  covariates.keep <- 50
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  # Set paths for model/results
  me.results.path <- file.path("artifacts/test_results/maxent",
                               paste0(spec, "_", st, "_maxent.rds"))
  me.path <- file.path("artifacts/models/maxent", 
                       paste0(spec, "_", st, "_maxent.rds"))
  if (!file.exists(me.results.path)) {
    cat("Starting MaxEnt for", spec, "in", st, "\n")
    # Load rasters and point data
    r <- rasters[[st]]
      
    Q <- readRDS(file.path("artifacts", "train_spatstat_Q_2",
                           paste0(st, "_", spec, "_Q.rds")))
    Q.test <- readRDS(file.path("artifacts", "test_spatstat_Q_2",
                                paste0(st, "_", spec, "_Q.rds")))
    p.train <- tibble(x=Q$data$x, y=Q$data$y) 
    # %>% sf::st_as_sf(coords = c("x", "y"), crs=4326)
    a.train <- tibble(x=Q$dummy$x, y=Q$dummy$y)
    p.test <- tibble(x=Q.test$data$x, y=Q.test$data$y) 
    a.test <- tibble(x=Q.test$dummy$x, y=Q.test$dummy$y)
      
    no.model.err <- F
    spec.sens.check <- F
    while (!no.model.err | !spec.sens.check) {
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
        covariates <- r[[covariates]]
        
        if (length(names(covariates)) < covariates.keep) {
          covariates.keep <- length(names(covariates)) 
        } else if (length(names(covariates)) > covariates.keep) {
          covariates <- covariates[[1:covariates.keep]]
        }
        
      } else {
        stop("\tThere are no specified covariates for", spec, st, "\n")
      }
      
      # Fit the MaxEnt Model
      cat("\tFitting MaxEnt Model...\n")
      fit.me <- tryCatch({
        predicts::MaxEnt(x=covariates, p=p.train) %>%
          get.object(
            obj=.,
            file.name=paste0(spec, "_", st, "_maxent.rds"), 
            obj.path="artifacts/models/maxent")}, 
        error=function(e) NULL)
      no.model.err <- !is.null(fit.me)
      if (!no.model.err) {
        file.remove(me.path)
        covariates.keep <- covariates.keep - 1
        if (covariates.keep < 1) {
          stop("\tUnable to successfully fit a model given the data.\n")
        } 
        next
      } 
      locations.train <- data.table::rbindlist(
        list(
          data.table(x=Q$data$x, y=Q$data$y, obs=T), 
          data.table(x=Q$dummy$x, y=Q$dummy$y, obs=F)
        )
      ) 
      purrr::walk(names(covariates), function(n) {
        locations.train[, (n) := 
                          terra::extract(covariates[[n]], 
                                         cbind(locations.train$x, 
                                               locations.train$y))]
      })
      me.pred <- predict(fit.me, locations.train)
      me.train <- cbind(locations.train, data.table(p.obs = me.pred))
      optimal.threshold <- optimize.f1(me.train)
      cm <- get.acc(me.train, optimal.threshold)
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
        covariates.keep <- covariates.keep - 1
        if (covariates.keep < 1) {
          stop("\tUnable to successfully fit a model given the data.\n")
        } 
        next
      } else {
        spec.sens.check <- T
      }
    }
    
    me.results <- get.object(
      obj={
        locations.test <- data.table::rbindlist(
          list(
            data.table(x=Q.test$data$x, y=Q.test$data$y, obs=T), 
            data.table(x=Q.test$dummy$x, y=Q.test$dummy$y, obs=F)
          )
        ) 
        purrr::walk(names(covariates), function(n) {
          locations.test[, (n) := 
                           terra::extract(covariates[[n]], 
                                          cbind(locations.test$x,
                                                locations.test$y))]
        })
        me.pred <- predict(fit.me, locations.test)
        me.test <- cbind(locations.test, data.table(p.obs=me.pred))
        cm <- get.acc(me.test, optimal.threshold)
        test.acc <- tibble(
          common.name=spec,
          state=st,
          covariate.count=covariates.keep,
          optimal.threshold=optimal.threshold 
        ) %>%
          cbind(as.list(c(cm$overall, cm$byClass)) %>% 
                  as_tibble()) %>%
          select(common.name:Accuracy, Sensitivity, Specificity, F1)
        all.predictions <- predict(fit.me, covariates)
        list(
          test=me.test,
          train=me.train,
          all.preds=all.predictions,
          thresh=optimal.threshold,
          train.accuracy=acc,
          test.accuracy=test.acc
        )
      },
      file.name=paste0(spec, "_", st, "_maxent.rds"),
      obj.path="artifacts/test_results/maxent"
    )
    cat("\tFinished MaxEnt model for", spec, "in", st, "\n")
  }
  gc()
})
```



## Results

```
me.models <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  results.path <- file.path("artifacts/test_results/maxent",
                        paste0(spec, "_", st, "_maxent.rds"))
  readRDS(results.path)$test.accuracy
})

DT::datatable(
  me.models,
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
  DT::formatStyle(columns=names(me.models), 
                  `font-size`="13px") %>%
  DT::formatSignif(4:ncol(me.models), digits=2)

```

### Load IPP Models for Comparison

```
ipp.models <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  results.path <- file.path("artifacts/test_results/ipp_glm_mpl_2",
                        paste0(spec, "_", st, "_ipp_glm_mpl_2.rds"))
  readRDS(results.path)$test.accuracy
})

```

### Overall Metric Summaries

```
# Reshape the data to long format
long.data <- me.models %>%
  select(common.name, Accuracy:F1) %>%
  tidyr::gather(key = "Metric", value = "Value", -common.name)

long.data.ipp <- ipp.models %>%
  select(common.name, Accuracy:F1) %>%
  tidyr::gather(key = "Metric", value = "Value", -common.name)

# Combine the old and new data frames, adding an identifier column
long.data$new.old <- "MaxEnt Models"
long.data.ipp$new.old <- "IPP Models"
combined.long <- rbind(long.data, long.data.ipp)

ggplot(combined.long, 
         aes(x = Metric, y = Value, fill = new.old)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(aes(group = new.old), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75),) +
    scale_fill_manual(values = c("#AAAAFF", "#FF9999"),
                      name = "Model Status",
                      labels = c("IPP Models", "MaxEnt Models")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "lightgray"),
          plot.background = element_rect(fill = "white")) + 
    coord_flip() +
    facet_wrap(~Metric, ncol=1, scales="free_y")


```

#### MaxEnt Models
```
summary(me.models)
```


#### IPP Models (For Reference)

```
summary(ipp.models)
```

### Summaries by Species

```
ggplot(combined.long, 
         aes(x = Metric, y = Value, fill = new.old)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(aes(group = new.old), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75),) +
    scale_fill_manual(values = c("#AAAAFF", "#FF9999"),
                      name = "Model Status",
                      labels = c("IPP Models", "MaxEnt Models")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "lightgray"),
          plot.background = element_rect(fill = "white")) + 
    coord_flip() +
    facet_wrap(~common.name + Metric, ncol=2, scales="free_y")

```


### Summaries by State

```
# Reshape the data to long format
long.data <- me.models %>%
  select(state, Accuracy:F1) %>%
  tidyr::gather(key = "Metric", value = "Value", -state)

long.data.ipp <- ipp.models %>%
  select(state, Accuracy:F1) %>%
  tidyr::gather(key = "Metric", value = "Value", -state)

# Combine the IPP and MaxEnt data frames, adding an identifier column
long.data$new.old <- "MaxEnt Models"
long.data.ipp$new.old <- "IPP Models"
combined.long <- rbind(long.data, long.data.ipp)

# Box plots for each state
ggplot(combined.long, 
         aes(x = Metric, y = Value, fill = new.old)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(aes(group = new.old), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75),) +
    scale_fill_manual(values = c("#AAAAFF", "#FF9999"),
                      name = "Model Status",
                      labels = c("IPP Models", "MaxEnt Models"))  +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "lightgray"),
          plot.background = element_rect(fill = "white")) + 
    coord_flip() +
    facet_wrap(~state + Metric, ncol=2, scales="free_y")

```

### Discussion

Comparing the MaxEnt and IPP model metrics of the fitted models on the
test data provides the following insights:

- Accuracy: The MaxEnt models generally exhibit high accuracy across most species and states, with several instances of perfect scores (1.00) particularly for Sanderling in Vermont. IPP models have varying accuracy, with a few instances comparable to MaxEnt but also showing lower values, especially for the Belted Kingfisher and Wild Turkey in North Carolina and Vermont.
- F1 Score: MaxEnt models demonstrate consistently high F1 scores, with the lowest being 0.735 for the Sandhill Crane in North Carolina. In contrast, IPP models show a wider range of F1 scores, with notable gaps (NA) where the model failed to predict certain species presence, such as for Sanderling in Colorado and Vermont.
- Specificity: Both model types generally exhibit high specificity, with many scores above 0.95. MaxEnt models achieve perfect specificity (1.00) for Sanderling in multiple states. IPP models have a few instances of lower specificity, particularly for Cedar Waxwing and Downy Woodpecker in Vermont.
- Sensitivity: MaxEnt models show strong sensitivity, with many values above 0.90. However, there are notable exceptions, such as the lower sensitivity for Belted Kingfisher in Vermont. IPP models have a broader range of sensitivity, with some species like Sanderling in North Carolina and Vermont, and Sandhill Crane in Vermont showing very low sensitivity, indicating a significant number of false negatives.

Overall, MaxEnt models seem to provide more consistent and reliable performance across all metrics, with particular strengths in specificity and sensitivity. IPP models, while effective in some cases, display a more variable performance and are notably less consistent in sensitivity.
