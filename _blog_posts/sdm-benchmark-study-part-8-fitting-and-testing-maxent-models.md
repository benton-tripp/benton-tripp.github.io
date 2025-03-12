---
layout: post
title: "SDM Benchmark Study Part 8: Fitting and Testing MaxEnt Models"
permalink: /blog_posts/sdm-benchmark-study-part-8-fitting-and-testing-maxent-models
---


## Overview

This part of the project fits additional "baseline" models using the MaxEnt
algorithm, as described in [Part 4](https://benton-tripp.github.io/blog/blog_posts/sdm-benchmark-study-part-4-baseline-species-distribution-models.html) 
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

<div style="min-width: 320px; overflow-x: auto; border: 1px solid #fff;">
  <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
  <script src="https://code.jquery.com/jquery-3.7.0.js"></script>
  <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>

  <div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item html-widget-static-bound" 
       id="htmlwidget-MaxEnt" style="min-width: 865px; height:auto; margin: 2px;">
    <table id="MaxEntResultsDT" class="display dataTable no-footer table table-condensed" style="width: 865px;">
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
    var MaxEntResults = [
      {
        "common.name": "Belted Kingfisher",
        "state": "CO",
        "covariate.count": 16,
        "optimal.threshold": 0.184552058800106,
        "Accuracy": 0.962564862861379,
        "Sensitivity": 0.944403261675315,
        "Specificity": 0.980726464047443,
        "F1": 0.961872404681012
      },
      {
        "common.name": "Cedar Waxwing",
        "state": "CO",
        "covariate.count": 16,
        "optimal.threshold": 0.185780935235688,
        "Accuracy": 0.962106299212598,
        "Sensitivity": 0.93996062992126,
        "Specificity": 0.984251968503937,
        "F1": 0.961248112732763
      },
      {
        "common.name": "Downy Woodpecker",
        "state": "CO",
        "covariate.count": 16,
        "optimal.threshold": 0.245153169127432,
        "Accuracy": 0.957513416815742,
        "Sensitivity": 0.936969155118462,
        "Specificity": 0.978076062639821,
        "F1": 0.956640803286171
      },
      {
        "common.name": "Ruddy Duck",
        "state": "CO",
        "covariate.count": 15,
        "optimal.threshold": 0.103966175774577,
        "Accuracy": 0.961923847695391,
        "Sensitivity": 0.959919839679359,
        "Specificity": 0.963927855711423,
        "F1": 0.961847389558233
      },
      {
        "common.name": "Sanderling",
        "state": "CO",
        "covariate.count": 12,
        "optimal.threshold": 0.0603805222313601,
        "Accuracy": 0.969565217391304,
        "Sensitivity": 0.766666666666667,
        "Specificity": 1,
        "F1": 0.867924528301887
      },
      {
        "common.name": "Sandhill Crane",
        "state": "CO",
        "covariate.count": 15,
        "optimal.threshold": 0.116401213500849,
        "Accuracy": 0.959550561797753,
        "Sensitivity": 0.961797752808989,
        "Specificity": 0.957303370786517,
        "F1": 0.959641255605381
      },
      {
        "common.name": "Sharp-shinned Hawk",
        "state": "CO",
        "covariate.count": 14,
        "optimal.threshold": 0.213252072043473,
        "Accuracy": 0.907984790874525,
        "Sensitivity": 0.882978723404255,
        "Specificity": 0.933028919330289,
        "F1": 0.905689789555729
      },
      {
        "common.name": "Wild Turkey",
        "state": "CO",
        "covariate.count": 15,
        "optimal.threshold": 0.211716730331628,
        "Accuracy": 0.963446475195822,
        "Sensitivity": 0.945169712793734,
        "Specificity": 0.981723237597911,
        "F1": 0.962765957446809
      },
      {
        "common.name": "Belted Kingfisher",
        "state": "NC",
        "covariate.count": 17,
        "optimal.threshold": 0.407736362711412,
        "Accuracy": 0.896188158961882,
        "Sensitivity": 0.815573770491803,
        "Specificity": 0.975120385232745,
        "F1": 0.886019590382903
      },
      {
        "common.name": "Cedar Waxwing",
        "state": "NC",
        "covariate.count": 16,
        "optimal.threshold": 0.329683799122103,
        "Accuracy": 0.90784155214228,
        "Sensitivity": 0.9,
        "Specificity": 0.915594855305466,
        "F1": 0.906633906633907
      },
      {
        "common.name": "Downy Woodpecker",
        "state": "NC",
        "covariate.count": 17,
        "optimal.threshold": 0.497636362632678,
        "Accuracy": 0.875960651706117,
        "Sensitivity": 0.802404438964242,
        "Specificity": 0.949110974862048,
        "F1": 0.865790786628971
      },
      {
        "common.name": "Ruddy Duck",
        "state": "NC",
        "covariate.count": 15,
        "optimal.threshold": 0.215651000066796,
        "Accuracy": 0.932707355242567,
        "Sensitivity": 0.907936507936508,
        "Specificity": 0.95679012345679,
        "F1": 0.930081300813008
      },
      {
        "common.name": "Sanderling",
        "state": "NC",
        "covariate.count": 15,
        "optimal.threshold": 0.175051035224783,
        "Accuracy": 0.977186311787072,
        "Sensitivity": 0.90625,
        "Specificity": 1,
        "F1": 0.950819672131147
      },
      {
        "common.name": "Sandhill Crane",
        "state": "NC",
        "covariate.count": 13,
        "optimal.threshold": 0.502912587845881,
        "Accuracy": 0.943478260869565,
        "Sensitivity": 0.6,
        "Specificity": 0.995,
        "F1": 0.73469387755102
      },
      {
        "common.name": "Sharp-shinned Hawk",
        "state": "NC",
        "covariate.count": 15,
        "optimal.threshold": 0.348160113501607,
        "Accuracy": 0.883977900552486,
        "Sensitivity": 0.821727019498607,
        "Specificity": 0.945205479452055,
        "F1": 0.875370919881306
      },
      {
        "common.name": "Wild Turkey",
        "state": "NC",
        "covariate.count": 15,
        "optimal.threshold": 0.466129369737495,
        "Accuracy": 0.870343839541547,
        "Sensitivity": 0.80863309352518,
        "Specificity": 0.931526390870185,
        "F1": 0.861302681992337
      },
      {
        "common.name": "Belted Kingfisher",
        "state": "OR",
        "covariate.count": 17,
        "optimal.threshold": 0.128586279826558,
        "Accuracy": 0.981353383458647,
        "Sensitivity": 0.981749528005035,
        "Specificity": 0.980990783410138,
        "F1": 0.98051539912005
      },
      {
        "common.name": "Cedar Waxwing",
        "state": "OR",
        "covariate.count": 17,
        "optimal.threshold": 0.182821607204132,
        "Accuracy": 0.977147520914099,
        "Sensitivity": 0.966890192791283,
        "Specificity": 0.986878727634195,
        "F1": 0.976301311891663
      },
      {
        "common.name": "Downy Woodpecker",
        "state": "OR",
        "covariate.count": 17,
        "optimal.threshold": 0.129584284881713,
        "Accuracy": 0.986161251504212,
        "Sensitivity": 0.983970406905055,
        "Specificity": 0.988249118683901,
        "F1": 0.9857936998147
      },
      {
        "common.name": "Ruddy Duck",
        "state": "OR",
        "covariate.count": 17,
        "optimal.threshold": 0.0810755499960402,
        "Accuracy": 0.981432360742706,
        "Sensitivity": 0.971929824561404,
        "Specificity": 0.991087344028521,
        "F1": 0.981399468556245
      },
      {
        "common.name": "Sanderling",
        "state": "OR",
        "covariate.count": 15,
        "optimal.threshold": 0.0884704644999837,
        "Accuracy": 0.995884773662551,
        "Sensitivity": 0.978260869565217,
        "Specificity": 1,
        "F1": 0.989010989010989
      },
      {
        "common.name": "Sandhill Crane",
        "state": "OR",
        "covariate.count": 16,
        "optimal.threshold": 0.118033459916917,
        "Accuracy": 0.964071856287425,
        "Sensitivity": 0.963276836158192,
        "Specificity": 0.964968152866242,
        "F1": 0.96600566572238
      },
      {
        "common.name": "Sharp-shinned Hawk",
        "state": "OR",
        "covariate.count": 17,
        "optimal.threshold": 0.132232673017086,
        "Accuracy": 0.964539007092199,
        "Sensitivity": 0.972332015810277,
        "Specificity": 0.957070707070707,
        "F1": 0.96407576747224
      },
      {
        "common.name": "Wild Turkey",
        "state": "OR",
        "covariate.count": 17,
        "optimal.threshold": 0.0755181175948422,
        "Accuracy": 0.991428571428571,
        "Sensitivity": 0.990196078431373,
        "Specificity": 0.992711370262391,
        "F1": 0.991584852734923
      },
      {
        "common.name": "Belted Kingfisher",
        "state": "VT",
        "covariate.count": 16,
        "optimal.threshold": 0.587385043837882,
        "Accuracy": 0.796511627906977,
        "Sensitivity": 0.629815745393635,
        "Specificity": 0.960461285008237,
        "F1": 0.754262788365095
      },
      {
        "common.name": "Cedar Waxwing",
        "state": "VT",
        "covariate.count": 15,
        "optimal.threshold": 0.150388443623785,
        "Accuracy": 0.782666666666667,
        "Sensitivity": 1,
        "Specificity": 0.00911854103343465,
        "F1": 0.877811094452774
      },
      {
        "common.name": "Downy Woodpecker",
        "state": "VT",
        "covariate.count": 16,
        "optimal.threshold": 0.141576534512486,
        "Accuracy": 0.754045307443366,
        "Sensitivity": 1,
        "Specificity": 0.0480167014613779,
        "F1": 0.857766687461011
      },
      {
        "common.name": "Ruddy Duck",
        "state": "VT",
        "covariate.count": 12,
        "optimal.threshold": 0.0935634581684863,
        "Accuracy": 0.981132075471698,
        "Sensitivity": 0.666666666666667,
        "Specificity": 1,
        "F1": 0.8
      },
      {
        "common.name": "Sanderling",
        "state": "VT",
        "covariate.count": 11,
        "optimal.threshold": 0.131230830996797,
        "Accuracy": 1,
        "Sensitivity": 1,
        "Specificity": 1,
        "F1": 1
      },
      {
        "common.name": "Sandhill Crane",
        "state": "VT",
        "covariate.count": 14,
        "optimal.threshold": 0.262539061372081,
        "Accuracy": 0.977168949771689,
        "Sensitivity": 0.789473684210526,
        "Specificity": 0.995,
        "F1": 0.857142857142857
      },
      {
        "common.name": "Sharp-shinned Hawk",
        "state": "VT",
        "covariate.count": 14,
        "optimal.threshold": 0.336412842244478,
        "Accuracy": 0.861678004535147,
        "Sensitivity": 0.821917808219178,
        "Specificity": 0.900900900900901,
        "F1": 0.855106888361045
      },
      {
        "common.name": "Wild Turkey",
        "state": "VT",
        "covariate.count": 16,
        "optimal.threshold": 0.511744152912645,
        "Accuracy": 0.821401077752117,
        "Sensitivity": 0.76661514683153,
        "Specificity": 0.875766871165644,
        "F1": 0.810457516339869
      }
    ];

    // Initialize the DataTable using the JSON dataset.
    jQuery(document).ready(function($){
      $('#MaxEntResultsDT').DataTable({
        data: MaxEntResults,
        columns: [
          { data: function(row) { return row["common.name"]; } },
          { data: "state" },
          { data: function(row) { return row["covariate.count"]; } },
          { 
            data: function(row) { return row["optimal.threshold"]; },
            render: function(data, type, row, meta) {
              return type === 'display' && !isNaN(data) ? parseFloat(data).toFixed(2) : data;
            }
          },
          { 
            data: "Accuracy",
            render: function(data, type, row, meta) {
              return type === 'display' && !isNaN(data) ? parseFloat(data).toFixed(2) : data;
            }
          },
          { 
            data: "Sensitivity",
            render: function(data, type, row, meta) {
              return type === 'display' && !isNaN(data) ? parseFloat(data).toFixed(2) : data;
            }
          },
          { 
            data: "Specificity",
            render: function(data, type, row, meta) {
              return type === 'display' && !isNaN(data) ? parseFloat(data).toFixed(2) : data;
            }
          },
          { 
            data: function(row) { return row["F1"] === null ? "" : row["F1"]; },
            render: function(data, type, row, meta) {
              return type === 'display' && data !== "" && !isNaN(data) ? parseFloat(data).toFixed(2) : data;
            }
          }
        ],
        scrollX: true,
        scrollY: true,
        paging: false,
        searching: false
      });
    });
  </script>
</div>

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

<img src="{{ site.baseurl }}/assets/plots/sdm-8-1.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 500px; height: auto;">


#### MaxEnt Models
```
summary(me.models)
```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
##  common.name           state           covariate.count optimal.threshold
##  Length:32          Length:32          Min.   :11.00   Min.   :0.06038  
##  Class :character   Class :character   1st Qu.:15.00   1st Qu.:0.12595  
##  Mode  :character   Mode  :character   Median :15.50   Median :0.18369  
##                                        Mean   :15.28   Mean   :0.23173  
##                                        3rd Qu.:16.25   3rd Qu.:0.33137  
##                                        Max.   :17.00   Max.   :0.58739  
##     Accuracy       Sensitivity      Specificity             F1        
##  Min.   :0.7540   Min.   :0.6000   Min.   :0.009118   Min.   :0.7347  
##  1st Qu.:0.8931   1st Qu.:0.8138   1st Qu.:0.948135   1st Qu.:0.8647  
##  Median :0.9620   Median :0.9385   Median :0.976598   Median :0.9405  
##  Mean   :0.9297   Mean   :0.8870   Mean   :0.909331   Mean   :0.9124  
##  3rd Qu.:0.9772   3rd Qu.:0.9720   3rd Qu.:0.991493   3rd Qu.:0.9646  
##  Max.   :1.0000   Max.   :1.0000   Max.   :1.000000   Max.   :1.0000
</code></pre></div><br>

#### IPP Models (For Reference)

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

<img src="{{ site.baseurl }}/assets/plots/sdm-8-2.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 800px; height: auto;">

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

<img src="{{ site.baseurl }}/assets/plots/sdm-8-3.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 800px; height: auto;">

### Discussion

Comparing the MaxEnt and IPP model metrics of the fitted models on the
test data provides the following insights:

- Accuracy: The MaxEnt models generally exhibit high accuracy across most species and states, with several instances of perfect scores (1.00) particularly for Sanderling in Vermont. IPP models have varying accuracy, with a few instances comparable to MaxEnt but also showing lower values, especially for the Belted Kingfisher and Wild Turkey in North Carolina and Vermont.
- F1 Score: MaxEnt models demonstrate consistently high F1 scores, with the lowest being 0.735 for the Sandhill Crane in North Carolina. In contrast, IPP models show a wider range of F1 scores, with notable gaps (NA) where the model failed to predict certain species presence, such as for Sanderling in Colorado and Vermont.
- Specificity: Both model types generally exhibit high specificity, with many scores above 0.95. MaxEnt models achieve perfect specificity (1.00) for Sanderling in multiple states. IPP models have a few instances of lower specificity, particularly for Cedar Waxwing and Downy Woodpecker in Vermont.
- Sensitivity: MaxEnt models show strong sensitivity, with many values above 0.90. However, there are notable exceptions, such as the lower sensitivity for Belted Kingfisher in Vermont. IPP models have a broader range of sensitivity, with some species like Sanderling in North Carolina and Vermont, and Sandhill Crane in Vermont showing very low sensitivity, indicating a significant number of false negatives.

Overall, MaxEnt models seem to provide more consistent and reliable performance across all metrics, with particular strengths in specificity and sensitivity. IPP models, while effective in some cases, display a more variable performance and are notably less consistent in sensitivity.
