---
layout: post
title: "SDM Benchmark Study Part 9: Fitting and Testing Common ML Models"
permalink: /blog_posts/sdm-benchmark-study-part-9-fitting-and-testing-ml-models
---

## Overview

Part 9 of the study focuses on common machine learning models, focusing on logistic regression,
KNN, classification tree, random forest, and XGBoost models. Following fitting the models,
comparisons of these models against the "traditional" IPP and MaxEnt models will be made. The
evaluation involves examining their performance through statistical tests, including ANOVA and
Kruskal-Wallis, followed by non-parametric pairwise comparisons and effect size assessments. This
comprehensive analysis aims to highlight the strengths and limitations of each model, providing
insights that could influence model selection in future research endeavors in this field.

## Setup

The setup is essentially the same as the prior few phases of the study. For the sake of 
brevity, much of the code has not been included in this particular section, but can still be
found within the [project Github repo](https://github.com/benton-tripp/presence-only-sdm/tree/main/R).

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

```
# Get data into normal data frame format (not `spatstat`)
data <- states %>%
  set_names() %>%
  purrr::map(function(st) {
  # Get raster by state
  r <- rasters[[st]]
  species %>% 
    set_names() %>%
    purrr::map(function(spec) {
      # Load `spatstat` quad data
      Q <- readRDS(file.path("artifacts", "train_spatstat_Q_2",
                             paste0(st, "_", spec, "_Q.rds")))
      Q.test <- readRDS(file.path("artifacts", "test_spatstat_Q_2",
                                  paste0(st, "_", spec, "_Q.rds")))
      # Load presence/absence data
      pres.train <- data.table(
        x=Q$data$x, 
        y=Q$data$y,
        presence=factor(T, levels=c(F,T), 
                        labels=c("Absence", "Presence")))
      abs.train <- data.table(
        x=Q$dummy$x, 
        y=Q$dummy$y, 
        presence=factor(F, levels=c(F,T),
                        labels=c("Absence", "Presence")))
      pres.test <- data.table(
        x=Q.test$data$x, 
        y=Q.test$data$y,
        presence=factor(T, levels=c(F,T),
                        labels=c("Absence", "Presence")))
      abs.test <- data.table(
        x=Q.test$dummy$x, 
        y=Q.test$dummy$y,
        presence=factor(F, levels=c(F,T),
                        labels=c("Absence", "Presence")))
      purrr::walk(names(r), function(n) {
        pres.train[, (n) := terra::extract(r[[n]], 
                                           cbind(pres.train$x,
                                                 pres.train$y))]
        abs.train[, (n) := terra::extract(r[[n]], 
                                          cbind(abs.train$x,
                                                abs.train$y))]
        pres.test[, (n) := terra::extract(r[[n]], 
                                          cbind(pres.test$x,
                                                pres.test$y))]
        abs.test[, (n) := terra::extract(r[[n]], 
                                         cbind(abs.test$x,
                                               abs.test$y))]
      })
      list(
        train=data.table::rbindlist(l=list(pres.train, abs.train)) %>%
          na.omit(),
        test=data.table::rbindlist(l=list(pres.test, abs.test)) %>%
          na.omit()
      )
    })
  })
```

## Model Fitting

In a similar manner to how the IPP and MaxEnt models were trained in the prior parts of the
study, each of the model types are iteratively fit over the different species and states.
5-fold cross-validation is used to help mitigate overfitting. The model and test results are 
saved after adjusting for optimal thresholds and handling cases of extreme sensitivity or
specificity.

```

control <- trainControl(method="cv", 
                        number=5, 
                        classProbs = T)

train.test <- function(model.type, fname, spec.state, data, var.imp.data, 
                       control, tune.grid=NULL, cov.keep=50) {
  purrr::walk(1:nrow(spec.state), function(i) {
    spec <- spec.state[i, ]$species
    st <- spec.state[i, ]$state
    fit.path <- file.path(paste0("artifacts/models/", fname),
                          paste0(spec, "_", st, "_", fname, ".rds"))
    results.path <- file.path(paste0("artifacts/test_results/", fname),
                              paste0(spec, "_", st, "_", fname, ".rds"))
    if (!file.exists(results.path)) {
      d <- data[[st]][[spec]]$train
      features <- var.imp.data[[st]][[spec]]
      covariates.keep <- cov.keep
      spec.sens.check <- F
      while (!spec.sens.check) {
        # Create formula
        feats <- features$variable[1:covariates.keep] %>%
          purrr::keep(~!is.na(.x))
        covariates.keep <- length(feats)
        # Create formula
        .f <- feats %>%
          paste(., collapse=" + ") %>% 
          paste("presence ~", .) %>% 
          as.formula()
        # Generalized for different model types
        if (model.type == "logistic regression") {
          fit <- train(.f, data = d, 
                       method = "glm", 
                       family = "binomial",
                       trControl = control, 
                       metric = "Accuracy")
        } else if (model.type %in% c("tree", "knn", "random forest", "xgboost")) {
          .method <- case_when(model.type == "tree" ~ "rpart",
                               model.type %in% c("knn") ~ model.type,
                               model.type == "random forest" ~ "ranger",
                               model.type == "xgboost" ~ "xgbTree",
                               T ~ "tree")
          fit <- train(.f, data=d, 
                 method = .method,
                 trControl = control, 
                 tuneGrid = tune.grid,
                 metric = "Accuracy")
        } else {
          stop("Invalid model type")
        }
        # Cache model object
        get.object(
          obj=fit,
          file.name=paste0(spec, "_", st, "_", fname, ".rds"), 
          obj.path=paste0("artifacts/models/", fname)) %>%
          suppressWarnings()
        
        pred <- predict(fit, d, type="prob")$Presence
        train.results <- cbind(
          d, data.table(obs = ifelse(d$presence == "Presence", T, F),
                        p.obs = pred))
        optimal.threshold <- optimize.f1(train.results)
        cm <- get.acc(train.results, optimal.threshold)
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
      
      results <- get.object(
        {
          d.test <- data[[st]][[spec]]$test
          pred.test <- predict(fit, d.test, type="prob")
          test.results <- cbind(d.test, 
                                data.table(obs = ifelse(
                                  d.test$presence == "Presence", T, F),
                                  p.obs = pred.test$Presence))
          cm <- get.acc(test.results, optimal.threshold)
          test.acc <- tibble(
            common.name=spec,
            state=st,
            covariate.count=covariates.keep,
            optimal.threshold=optimal.threshold 
          ) %>%
            cbind(as.list(c(cm$overall, cm$byClass)) %>% 
                    as_tibble()) %>%
            select(common.name:Accuracy, Sensitivity, Specificity, F1)
          all.predictions <- rasters[[st]] %>% 
            as.data.table() %>%
            predict(fit, ., type="prob") %>%
            .$Presence
          list(
            test=test.results,
            train=train.results,
            all.preds=all.predictions,
            thresh=optimal.threshold,
            train.accuracy=acc,
            test.accuracy=test.acc
          )
        },
        file.name=paste0(spec, "_", st, "_", fname, ".rds"),
        obj.path=paste0("artifacts/test_results/", fname)
      )
      cat("\tFinished model for", spec, "in", st, "\n")
    }
  })
}

```

### Logistic Regression

```

train.test(model.type="logistic regression", 
           fname="logreg", 
           spec.state, 
           data, 
           var.imp.data,
           control, 
           cov.keep=50)

```

### Classification Tree

```

train.test(model.type="tree", 
           fname="tree", 
           spec.state, 
           data, 
           var.imp.data,
           control, 
           tune.grid=data.frame(cp=seq(0, 0.1, by = 0.01)),
           cov.keep=50)

```

### K-Nearest Neighbors

```

train.test(model.type="knn", 
           fname="knn", 
           spec.state, 
           data, 
           var.imp.data,
           control, 
           tune.grid=data.frame(k=1:10),
           cov.keep=50)

```

### Random Forest

```

train.test(model.type="random forest", 
           fname="randomforest", 
           spec.state, 
           data, 
           var.imp.data,
           control, 
           tune.grid=expand.grid(
             splitrule=c("extratrees", "gini"),
             mtry=c(1,5,10),
             min.node.size=c(1,5,10)
           ),
           cov.keep=50)

```

### XGBoost

```

train.test(model.type="xgboost", 
           fname="xgboost", 
           spec.state, 
           data, 
           var.imp.data,
           control, 
           tune.grid=expand.grid(
             nrounds = 300,
             eta = 0.1,
             gamma = c(0, 0.25, 0.5),
             max_depth = c(3, 6, 9),
             colsample_bytree = c(0.5, 1),
             min_child_weight = 1,
             subsample = c(0.5, 1)
           ),
           cov.keep=50) 

```

## Results

```

model.scores <- purrr::map_df(c("ipp_glm_mpl_2", "maxent", 
                          "logreg", "knn", "tree", 
                          "randomforest", "xgboost"), function(m) {
                            m.type <- case_when(
                              m == "ipp_glm_mpl_2" ~ "ipp",
                              m == "maxent" ~ "maxent",
                              m == "logreg" ~ "logistic regression",
                              m == "knn" ~ "knn",
                              m == "tree" ~ "classification tree",
                              m == "randomforest" ~ "random forest",
                              m == "xgboost" ~ "xgboost",
                              T ~ "")
  purrr::map_df(1:nrow(spec.state), function(i) {
    spec <- spec.state[i,]$species
    st <- spec.state[i,]$state
    results.path <- file.path(paste0("artifacts/test_results/", m),
                              paste0(spec, "_", st, "_", m, ".rds"))
    readRDS(results.path)$test.accuracy
  }) %>% mutate(model.type=m.type)
})


reshape.data <- function(data, metric) {
  # Ensure that the metric is one of the expected columns
  if (!metric %in% c("Accuracy", "Sensitivity", 
                     "Specificity", "F1")) {
    stop("Invalid metric specified.")
  }
  # Reshape the data
  data %>%
    select(common.name, state, model.type, !!sym(metric)) %>%
    mutate(model.type = tolower(gsub(" ", ".", model.type))) %>%
    tidyr::pivot_wider(names_from = model.type, 
                       values_from = !!sym(metric),
                names_sep=".",
                names_prefix=paste0(tolower(metric), "."))
}

metrics <- c("Accuracy", "Sensitivity", "Specificity", "F1") %>%
  set_names() %>%
  purrr::map(~reshape.data(model.scores, .x))

dt.metrics <- c("Accuracy", "Sensitivity", "Specificity", "F1") %>%
  set_names() %>%
  purrr::map(~{
    d <- metrics[[.x]]
    DT::datatable(
      d,
      filter='none',
      selection='none',
      rownames=F,
      options=list(
        scrollX=T,
        scrollY="600px",
        paging=F,
        searching=F,
        orderMulti=T,
        info=F,
        lengthChange = F
      )) %>%
      DT::formatStyle(columns=names(d), 
                      `font-size`="13px") %>%
      DT::formatSignif(3:ncol(d), digits=4) %>%
      htmltools::div(id=paste0(.x, "_datatable"), 
                     style=ifelse(.x=="Accuracy", "", 
                                  "visibility:hidden; height:0;"), .)
  })

htmltools::div(
  htmltools::tags$script(
    '$(document).ready(function(){
        $("#metric_selector").change(function(){
          var selectedMetric = $(this).val();
          // Hide datatbles
          $("[id$=_datatable]").css({"visibility": "hidden", "height": 0});
          // Show the selected datatable
          $("#" + selectedMetric + "_datatable").css({"visibility": "visible",
                                                      "height": "auto"});
        });
      });'
  ),
  htmltools::tags$select(id='metric_selector',
                         lapply(names(dt.metrics), function(met) {
                           htmltools::tags$option(value=met, met)
                         })), dt.metrics)

```

<div style="min-width: 320px; overflow-x: auto; border: 1px solid #fff; height: 720px; overflow-y: hidden;">
  <!-- Include DataTables CSS and JS -->
  <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
  <script src="https://code.jquery.com/jquery-3.7.0.js"></script>
  <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>

  <!-- Dropdown filter to select which table to show -->
  <script>
    $(document).ready(function(){
      $("#metric_selector").change(function(){
        var selectedMetric = $(this).val();
        // Hide all datatable containers
        $("#Accuracy_datatable, #Sensitivity_datatable, #Specificity_datatable, #F1_datatable").hide();
        // Show the selected container
        $("#" + selectedMetric + "_datatable").show();
      });
    });
  </script>
  <select id="metric_selector" style="margin-top: 1px; margin-left: 1px;" >
    <option value="Accuracy">Accuracy</option>
    <option value="Sensitivity">Sensitivity</option>
    <option value="Specificity">Specificity</option>
    <option value="F1">F1</option>
  </select>

  <!-- Accuracy Table -->
  <div id="Accuracy_datatable" style="width: 100%; margin-top: 15px;">
    <div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item html-widget-static-bound" id="htmlwidget-accuracy" style="width:100%;height:auto;">
      <table id="DataTables_Accuracy" class="display dataTable no-footer table table-condensed" style="width:100%;">
        <thead>
          <tr>
            <th>common.name</th>
            <th>state</th>
            <th>accuracy.ipp</th>
            <th>accuracy.maxent</th>
            <th>accuracy.logistic.regression</th>
            <th>accuracy.knn</th>
            <th>accuracy.classification.tree</th>
            <th>accuracy.random.forest</th>
            <th>accuracy.xgboost</th>
          </tr>
        </thead>
      </table>
    </div>
    <script>
      // Define the JSON dataset for Accuracy (replace with your actual data)
      var AccuracyResults = [
        {
          "common.name": "Belted Kingfisher",
          "state": "CO",
          "accuracy.ipp": 0.979985174203113,
          "accuracy.maxent": 0.962564862861379,
          "accuracy.logistic.regression": 0.986286137879911,
          "accuracy.knn": 0.984803558191253,
          "accuracy.classification.tree": 0.967383246849518,
          "accuracy.random.forest": 0.986286137879911,
          "accuracy.xgboost": 0.988139362490734
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "CO",
          "accuracy.ipp": 0.952755905511811,
          "accuracy.maxent": 0.962106299212598,
          "accuracy.logistic.regression": 0.972933070866142,
          "accuracy.knn": 0.975885826771654,
          "accuracy.classification.tree": 0.949311023622047,
          "accuracy.random.forest": 0.974409448818898,
          "accuracy.xgboost": 0.982775590551181
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "CO",
          "accuracy.ipp": 0.974284436493739,
          "accuracy.maxent": 0.957513416815742,
          "accuracy.logistic.regression": 0.98345259391771,
          "accuracy.knn": 0.985912343470483,
          "accuracy.classification.tree": 0.974060822898032,
          "accuracy.random.forest": 0.986359570661896,
          "accuracy.xgboost": 0.98613595706619
        },
        {
          "common.name": "Ruddy Duck",
          "state": "CO",
          "accuracy.ipp": 0.969939879759519,
          "accuracy.maxent": 0.961923847695391,
          "accuracy.logistic.regression": 0.977955911823647,
          "accuracy.knn": 0.987975951903808,
          "accuracy.classification.tree": 0.963927855711423,
          "accuracy.random.forest": 0.98496993987976,
          "accuracy.xgboost": 0.987975951903808
        },
        {
          "common.name": "Sanderling",
          "state": "CO",
          "accuracy.ipp": 0.869565217391304,
          "accuracy.maxent": 0.969565217391304,
          "accuracy.logistic.regression": 0.960869565217391,
          "accuracy.knn": 0.947826086956522,
          "accuracy.classification.tree": 0.847826086956522,
          "accuracy.random.forest": 0.973913043478261,
          "accuracy.xgboost": 0.830434782608696
        },
        {
          "common.name": "Sandhill Crane",
          "state": "CO",
          "accuracy.ipp": 0.959550561797753,
          "accuracy.maxent": 0.959550561797753,
          "accuracy.logistic.regression": 0.98876404494382,
          "accuracy.knn": 0.984269662921348,
          "accuracy.classification.tree": 0.969662921348315,
          "accuracy.random.forest": 0.99438202247191,
          "accuracy.xgboost": 0.991011235955056
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "CO",
          "accuracy.ipp": 0.88212927756654,
          "accuracy.maxent": 0.907984790874525,
          "accuracy.logistic.regression": 0.959695817490494,
          "accuracy.knn": 0.968821292775665,
          "accuracy.classification.tree": 0.944486692015209,
          "accuracy.random.forest": 0.965779467680608,
          "accuracy.xgboost": 0.968060836501901
        },
        {
          "common.name": "Wild Turkey",
          "state": "CO",
          "accuracy.ipp": 0.963446475195822,
          "accuracy.maxent": 0.963446475195822,
          "accuracy.logistic.regression": 0.97911227154047,
          "accuracy.knn": 0.984334203655353,
          "accuracy.classification.tree": 0.969321148825065,
          "accuracy.random.forest": 0.988250652741514,
          "accuracy.xgboost": 0.987597911227154
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "NC",
          "accuracy.ipp": 0.758103727714749,
          "accuracy.maxent": 0.896188158961882,
          "accuracy.logistic.regression": 0.939983779399838,
          "accuracy.knn": 0.951743714517437,
          "accuracy.classification.tree": 0.913625304136253,
          "accuracy.random.forest": 0.959448499594485,
          "accuracy.xgboost": 0.963884430176565
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "NC",
          "accuracy.ipp": 0.882020202020202,
          "accuracy.maxent": 0.90784155214228,
          "accuracy.logistic.regression": 0.953516572352466,
          "accuracy.knn": 0.959175424413905,
          "accuracy.classification.tree": 0.906633906633907,
          "accuracy.random.forest": 0.962813257881972,
          "accuracy.xgboost": 0.971061093247588
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "NC",
          "accuracy.ipp": 0.864607345935147,
          "accuracy.maxent": 0.875960651706117,
          "accuracy.logistic.regression": 0.900860743928681,
          "accuracy.knn": 0.932216415616354,
          "accuracy.classification.tree": 0.865790786628971,
          "accuracy.random.forest": 0.947740547187212,
          "accuracy.xgboost": 0.93408951563458
        },
        {
          "common.name": "Ruddy Duck",
          "state": "NC",
          "accuracy.ipp": 0.861370716510903,
          "accuracy.maxent": 0.932707355242567,
          "accuracy.logistic.regression": 0.956181533646322,
          "accuracy.knn": 0.907936507936508,
          "accuracy.classification.tree": 0.928012519561815,
          "accuracy.random.forest": 0.966049382716049,
          "accuracy.xgboost": 0.960876369327074
        },
        {
          "common.name": "Sanderling",
          "state": "NC",
          "accuracy.ipp": 0.901140684410646,
          "accuracy.maxent": 0.977186311787072,
          "accuracy.logistic.regression": 0.973384030418251,
          "accuracy.knn": 0.90625,
          "accuracy.classification.tree": 0.89922480620155,
          "accuracy.random.forest": 0.964824120603015,
          "accuracy.xgboost": 0.994974874371859
        },
        {
          "common.name": "Sandhill Crane",
          "state": "NC",
          "accuracy.ipp": 0.865217391304348,
          "accuracy.maxent": 0.943478260869565,
          "accuracy.logistic.regression": 0.908695652173913,
          "accuracy.knn": 0.6,
          "accuracy.classification.tree": 0.81304347826087,
          "accuracy.random.forest": 0.73469387755102,
          "accuracy.xgboost": 0.92
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "NC",
          "accuracy.ipp": 0.873103448275862,
          "accuracy.maxent": 0.883977900552486,
          "accuracy.logistic.regression": 0.937845303867403,
          "accuracy.knn": 0.821727019498607,
          "accuracy.classification.tree": 0.893646408839779,
          "accuracy.random.forest": 0.875370919881306,
          "accuracy.xgboost": 0.953038674033149
        },
        {
          "common.name": "Wild Turkey",
          "state": "NC",
          "accuracy.ipp": 0.865329512893983,
          "accuracy.maxent": 0.870343839541547,
          "accuracy.logistic.regression": 0.934813753581662,
          "accuracy.knn": 0.80863309352518,
          "accuracy.classification.tree": 0.908309455587393,
          "accuracy.random.forest": 0.861302681992337,
          "accuracy.xgboost": 0.98051539912005
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "OR",
          "accuracy.ipp": 0.98136459272618,
          "accuracy.maxent": 0.981353383458647,
          "accuracy.logistic.regression": 0.989172932330827,
          "accuracy.knn": 0.991578947368421,
          "accuracy.classification.tree": 0.975338345864662,
          "accuracy.random.forest": 0.988571428571429,
          "accuracy.xgboost": 0.986751152073733
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "OR",
          "accuracy.ipp": 0.947765762089369,
          "accuracy.maxent": 0.977147520914099,
          "accuracy.logistic.regression": 0.994082840236686,
          "accuracy.knn": 0.988866799204771,
          "accuracy.classification.tree": 0.992450520301979,
          "accuracy.random.forest": 0.988866799204771,
          "accuracy.xgboost": 0.977147520914099
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "OR",
          "accuracy.ipp": 0.987367154601965,
          "accuracy.maxent": 0.986161251504212,
          "accuracy.logistic.regression": 0.997593261131167,
          "accuracy.knn": 0.988249118683901,
          "accuracy.classification.tree": 0.990974729241877,
          "accuracy.random.forest": 0.997593261131167,
          "accuracy.xgboost": 0.99529964747356
        },
        {
          "common.name": "Ruddy Duck",
          "state": "OR",
          "accuracy.ipp": 0.957559681697613,
          "accuracy.maxent": 0.981432360742706,
          "accuracy.logistic.regression": 0.983200707338638,
          "accuracy.knn": 0.981432360742706,
          "accuracy.classification.tree": 0.992042440318302,
          "accuracy.random.forest": 0.981432360742706,
          "accuracy.xgboost": 0.995884773662551
        },
        {
          "common.name": "Sanderling",
          "state": "OR",
          "accuracy.ipp": 0.872427983539095,
          "accuracy.maxent": 0.995884773662551,
          "accuracy.logistic.regression": 0.983539094650206,
          "accuracy.knn": 0.978260869565217,
          "accuracy.classification.tree": 0.963276836158192,
          "accuracy.random.forest": 0.976114649681529,
          "accuracy.xgboost": 0.989847715736041
        },
        {
          "common.name": "Sandhill Crane",
          "state": "OR",
          "accuracy.ipp": 0.972305389221557,
          "accuracy.maxent": 0.964071856287425,
          "accuracy.logistic.regression": 0.993263473053892,
          "accuracy.knn": 0.964071856287425,
          "accuracy.classification.tree": 0.991525423728814,
          "accuracy.random.forest": 0.964071856287425,
          "accuracy.xgboost": 0.977168949771689
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "OR",
          "accuracy.ipp": 0.859445519019987,
          "accuracy.maxent": 0.964539007092199,
          "accuracy.logistic.regression": 0.985815602836879,
          "accuracy.knn": 0.991428571428571,
          "accuracy.classification.tree": 0.991618310767247,
          "accuracy.random.forest": 0.995,
          "accuracy.xgboost": 0.861678004535147
        },
        {
          "common.name": "Wild Turkey",
          "state": "OR",
          "accuracy.ipp": 0.91220556745182,
          "accuracy.maxent": 0.821401077752117,
          "accuracy.logistic.regression": 0.987857142857143,
          "accuracy.knn": 0.997142857142857,
          "accuracy.classification.tree": 0.997142857142857,
          "accuracy.random.forest": 0.977168949771689,
          "accuracy.xgboost": 0.861678004535147
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "VT",
          "accuracy.ipp": 0.803156146179402,
          "accuracy.maxent": 0.796511627906977,
          "accuracy.logistic.regression": 0.8023,
          "accuracy.knn": 0.7234,
          "accuracy.classification.tree": 0.7500,
          "accuracy.random.forest": 0.7575,
          "accuracy.xgboost": 0.730897009966777
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "VT",
          "accuracy.ipp": 0.797333333333333,
          "accuracy.maxent": 0.782666666666667,
          "accuracy.logistic.regression": 0.7760,
          "accuracy.knn": 0.7687,
          "accuracy.classification.tree": 0.7567,
          "accuracy.random.forest": 0.7880,
          "accuracy.xgboost": 0.764
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "VT",
          "accuracy.ipp": 0.752427184466019,
          "accuracy.maxent": 0.754045307443366,
          "accuracy.logistic.regression": 0.7433,
          "accuracy.knn": 0.7816,
          "accuracy.classification.tree": 0.7740,
          "accuracy.random.forest": 0.7740,
          "accuracy.xgboost": 0.7859
        },
        {
          "common.name": "Ruddy Duck",
          "state": "VT",
          "accuracy.ipp": 0.943396226415094,
          "accuracy.maxent": 0.981132075471698,
          "accuracy.logistic.regression": 0.9717,
          "accuracy.knn": 0.9481,
          "accuracy.classification.tree": 0.9481,
          "accuracy.random.forest": 0.9906,
          "accuracy.xgboost": 0.9481
        },
        {
          "common.name": "Sanderling",
          "state": "VT",
          "accuracy.ipp": 0.947619047619048,
          "accuracy.maxent": 1,
          "accuracy.logistic.regression": 0.9562,
          "accuracy.knn": 0.9124,
          "accuracy.classification.tree": 0.9280,
          "accuracy.random.forest": 0.9656,
          "accuracy.xgboost": 0.9609
        },
        {
          "common.name": "Sandhill Crane",
          "state": "VT",
          "accuracy.ipp": 0.91324200913242,
          "accuracy.maxent": 0.977168949771689,
          "accuracy.logistic.regression": 0.9832,
          "accuracy.knn": 0.9912,
          "accuracy.classification.tree": 0.9912,
          "accuracy.random.forest": 0.9912,
          "accuracy.xgboost": 0.9930
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "VT",
          "accuracy.ipp": 0.870748299319728,
          "accuracy.maxent": 0.861678004535147,
          "accuracy.logistic.regression": 0.959183673469388,
          "accuracy.knn": 0.952380952380952,
          "accuracy.classification.tree": 0.936507936507937,
          "accuracy.random.forest": 0.649730561970747,
          "accuracy.xgboost": 0.952380952380952
        },
        {
          "common.name": "Wild Turkey",
          "state": "VT",
          "accuracy.ipp": 0.709006928406467,
          "accuracy.maxent": 0.821401077752117,
          "accuracy.logistic.regression": 0.648960739030023,
          "accuracy.knn": 0.659738260200154,
          "accuracy.classification.tree": 0.696689761354888,
          "accuracy.random.forest": 0.649730561970747,
          "accuracy.xgboost": 0.68437259430331
        }
      ];

      jQuery(document).ready(function($){
        $('#DataTables_Accuracy').DataTable({
          data: AccuracyResults,
          columns: [
            { data: function(row) { return row["common.name"]; } },
            { data: "state" },
            { data: function(row) { return row["accuracy.ipp"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["accuracy.maxent"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["accuracy.logistic.regression"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["accuracy.knn"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["accuracy.classification.tree"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["accuracy.random.forest"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["accuracy.xgboost"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }}
          ],
          scrollX: true,
          scrollY: "600px",
          paging: false,
          searching: false
        });
      });
    </script>
  </div>

  <!-- Sensitivity Table -->
  <div id="Sensitivity_datatable" style="width: 100%; margin-top: 15px; display: none;">
    <div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item html-widget-static-bound" id="htmlwidget-sensitivity" style="width:100%;height:auto;">
      <table id="DataTables_Sensitivity" class="display dataTable no-footer table table-condensed" style="width:100%;">
        <thead>
          <tr>
            <th>common.name</th>
            <th>state</th>
            <th>sensitivity.ipp</th>
            <th>sensitivity.maxent</th>
            <th>sensitivity.logistic.regression</th>
            <th>sensitivity.knn</th>
            <th>sensitivity.classification.tree</th>
            <th>sensitivity.random.forest</th>
            <th>sensitivity.xgboost</th>
          </tr>
        </thead>
      </table>
    </div>
    <script>
      // Define the JSON dataset for Sensitivity (replace with your actual data)
      var SensitivityResults = [
        {
          "common.name": "Belted Kingfisher",
          "state": "CO",
          "sensitivity.ipp": 0.99184581171238,
          "sensitivity.maxent": 0.980726464047443,
          "sensitivity.logistic.regression": 0.989621942179392,
          "sensitivity.knn": 0.984803558191253,
          "sensitivity.classification.tree": 0.967383246849518,
          "sensitivity.random.forest": 0.986286137879911,
          "sensitivity.xgboost": 0.988139362490734
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "CO",
          "sensitivity.ipp": 0.982283464566929,
          "sensitivity.maxent": 0.984251968503937,
          "sensitivity.logistic.regression": 0.984251968503937,
          "sensitivity.knn": 0.975885826771654,
          "sensitivity.classification.tree": 0.949311023622047,
          "sensitivity.random.forest": 0.974409448818898,
          "sensitivity.xgboost": 0.982775590551181
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "CO",
          "sensitivity.ipp": 0.993736017897092,
          "sensitivity.maxent": 0.978076062639821,
          "sensitivity.logistic.regression": 0.987919463087248,
          "sensitivity.knn": 0.985912343470483,
          "sensitivity.classification.tree": 0.974060822898032,
          "sensitivity.random.forest": 0.986359570661896,
          "sensitivity.xgboost": 0.98613595706619
        },
        {
          "common.name": "Ruddy Duck",
          "state": "CO",
          "sensitivity.ipp": 0.98997995991984,
          "sensitivity.maxent": 0.963927855711423,
          "sensitivity.logistic.regression": 0.981963927855711,
          "sensitivity.knn": 0.987975951903808,
          "sensitivity.classification.tree": 0.963927855711423,
          "sensitivity.random.forest": 0.98496993987976,
          "sensitivity.xgboost": 0.987975951903808
        },
        {
          "common.name": "Sanderling",
          "state": "CO",
          "sensitivity.ipp": 1,
          "sensitivity.maxent": 1,
          "sensitivity.logistic.regression": 0.975,
          "sensitivity.knn": 0.947826086956522,
          "sensitivity.classification.tree": 0.847826086956522,
          "sensitivity.random.forest": 0.973913043478261,
          "sensitivity.xgboost": 0.830434782608696
        },
        {
          "common.name": "Sandhill Crane",
          "state": "CO",
          "sensitivity.ipp": 0.993258426966292,
          "sensitivity.maxent": 0.957303370786517,
          "sensitivity.logistic.regression": 0.993258426966292,
          "sensitivity.knn": 0.984269662921348,
          "sensitivity.classification.tree": 0.969662921348315,
          "sensitivity.random.forest": 0.99438202247191,
          "sensitivity.xgboost": 0.991011235955056
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "CO",
          "sensitivity.ipp": 0.987823439878234,
          "sensitivity.maxent": 0.933028919330289,
          "sensitivity.logistic.regression": 0.969558599695586,
          "sensitivity.knn": 0.968821292775665,
          "sensitivity.classification.tree": 0.944486692015209,
          "sensitivity.random.forest": 0.965779467680608,
          "sensitivity.xgboost": 0.968060836501901
        },
        {
          "common.name": "Wild Turkey",
          "state": "CO",
          "sensitivity.ipp": 0.951697127937337,
          "sensitivity.maxent": 0.981723237597911,
          "sensitivity.logistic.regression": 0.989556135770235,
          "sensitivity.knn": 0.984334203655353,
          "sensitivity.classification.tree": 0.969321148825065,
          "sensitivity.random.forest": 0.988250652741514,
          "sensitivity.xgboost": 0.987597911227154
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "NC",
          "sensitivity.ipp": 0.963884430176565,
          "sensitivity.maxent": 0.975120385232745,
          "sensitivity.logistic.regression": 0.965489566613162,
          "sensitivity.knn": 0.951743714517437,
          "sensitivity.classification.tree": 0.913625304136253,
          "sensitivity.random.forest": 0.959448499594485,
          "sensitivity.xgboost": 0.963884430176565
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "NC",
          "sensitivity.ipp": 0.971061093247588,
          "sensitivity.maxent": 0.915594855305466,
          "sensitivity.logistic.regression": 0.973472668810289,
          "sensitivity.knn": 0.959175424413905,
          "sensitivity.classification.tree": 0.906633906633907,
          "sensitivity.random.forest": 0.962813257881972,
          "sensitivity.xgboost": 0.971061093247588
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "NC",
          "sensitivity.ipp": 0.93408951563458,
          "sensitivity.maxent": 0.949110974862048,
          "sensitivity.logistic.regression": 0.964132434089516,
          "sensitivity.knn": 0.932216415616354,
          "sensitivity.classification.tree": 0.865790786628971,
          "sensitivity.random.forest": 0.947740547187212,
          "sensitivity.xgboost": 0.93408951563458
        },
        {
          "common.name": "Ruddy Duck",
          "state": "NC",
          "sensitivity.ipp": 0.861370716510903,
          "sensitivity.maxent": 0.932707355242567,
          "sensitivity.logistic.regression": 0.956181533646322,
          "sensitivity.knn": 0.907936507936508,
          "sensitivity.classification.tree": 0.928012519561815,
          "sensitivity.random.forest": 0.966049382716049,
          "sensitivity.xgboost": 0.960876369327074
        },
        {
          "common.name": "Sanderling",
          "state": "NC",
          "sensitivity.ipp": 0.994974874371859,
          "sensitivity.maxent": 1,
          "sensitivity.logistic.regression": 0.964824120603015,
          "sensitivity.knn": 0.992395437262357,
          "sensitivity.classification.tree": 0.92,
          "sensitivity.random.forest": 0.955,
          "sensitivity.xgboost": 0.994974874371859
        },
        {
          "common.name": "Sandhill Crane",
          "state": "NC",
          "sensitivity.ipp": 0.995,
          "sensitivity.maxent": 0.945205479452055,
          "sensitivity.logistic.regression": 0.975342465753425,
          "sensitivity.knn": 0.956521739130435,
          "sensitivity.classification.tree": 0.833333333333333,
          "sensitivity.random.forest": 0.995,
          "sensitivity.xgboost": 0.945205479452055
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "NC",
          "sensitivity.ipp": 0.947945205479452,
          "sensitivity.maxent": 0.967189728958631,
          "sensitivity.logistic.regression": 0.998271889400922,
          "sensitivity.knn": 0.988866799204771,
          "sensitivity.classification.tree": 0.92,
          "sensitivity.random.forest": 0.987522281639929,
          "sensitivity.xgboost": 0.995
        },
        {
          "common.name": "Wild Turkey",
          "state": "NC",
          "sensitivity.ipp": 0.958630527817404,
          "sensitivity.maxent": 0.931526390870185,
          "sensitivity.logistic.regression": 0.987522281639929,
          "sensitivity.knn": 0.991087344028521,
          "sensitivity.classification.tree": 1,
          "sensitivity.random.forest": 0.989847715736041,
          "sensitivity.xgboost": 0.906095551894563
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "OR",
          "sensitivity.ipp": 0.986751152073733,
          "sensitivity.maxent": 0.980990783410138,
          "sensitivity.logistic.regression": 0.998271889400922,
          "sensitivity.knn": 0.991578947368421,
          "sensitivity.classification.tree": 0.988571428571429,
          "sensitivity.random.forest": 0.988571428571429,
          "sensitivity.xgboost": 0.986751152073733
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "OR",
          "sensitivity.ipp": 0.947765762089369,
          "sensitivity.maxent": 0.977147520914099,
          "sensitivity.logistic.regression": 0.994082840236686,
          "sensitivity.knn": 0.988866799204771,
          "sensitivity.classification.tree": 0.992450520301979,
          "sensitivity.random.forest": 0.988866799204771,
          "sensitivity.xgboost": 0.977147520914099
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "OR",
          "sensitivity.ipp": 0.987367154601965,
          "sensitivity.maxent": 0.986161251504212,
          "sensitivity.logistic.regression": 0.997593261131167,
          "sensitivity.knn": 0.988249118683901,
          "sensitivity.classification.tree": 0.990974729241877,
          "sensitivity.random.forest": 0.997593261131167,
          "sensitivity.xgboost": 0.99529964747356
        },
        {
          "common.name": "Ruddy Duck",
          "state": "OR",
          "sensitivity.ipp": 0.957559681697613,
          "sensitivity.maxent": 0.981432360742706,
          "sensitivity.logistic.regression": 0.983200707338638,
          "sensitivity.knn": 0.981432360742706,
          "sensitivity.classification.tree": 0.992042440318302,
          "sensitivity.random.forest": 0.981432360742706,
          "sensitivity.xgboost": 0.995884773662551
        },
        {
          "common.name": "Sanderling",
          "state": "OR",
          "sensitivity.ipp": 0.872427983539095,
          "sensitivity.maxent": 0.995884773662551,
          "sensitivity.logistic.regression": 0.983539094650206,
          "sensitivity.knn": 0.978260869565217,
          "sensitivity.classification.tree": 0.963276836158192,
          "sensitivity.random.forest": 0.976114649681529,
          "sensitivity.xgboost": 0.989847715736041
        },
        {
          "common.name": "Sandhill Crane",
          "state": "OR",
          "sensitivity.ipp": 0.972305389221557,
          "sensitivity.maxent": 0.964071856287425,
          "sensitivity.logistic.regression": 0.993263473053892,
          "sensitivity.knn": 0.964071856287425,
          "sensitivity.classification.tree": 0.991525423728814,
          "sensitivity.random.forest": 0.964071856287425,
          "sensitivity.xgboost": 0.977168949771689
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "OR",
          "sensitivity.ipp": 0.859445519019987,
          "sensitivity.maxent": 0.964539007092199,
          "sensitivity.logistic.regression": 0.985815602836879,
          "sensitivity.knn": 0.991428571428571,
          "sensitivity.classification.tree": 0.991618310767247,
          "sensitivity.random.forest": 0.995,
          "sensitivity.xgboost": 0.861678004535147
        },
        {
          "common.name": "Wild Turkey",
          "state": "OR",
          "sensitivity.ipp": 0.91220556745182,
          "sensitivity.maxent": 0.821401077752117,
          "sensitivity.logistic.regression": 0.987857142857143,
          "sensitivity.knn": 0.997142857142857,
          "sensitivity.classification.tree": 0.997142857142857,
          "sensitivity.random.forest": 0.977168949771689,
          "sensitivity.xgboost": 0.861678004535147
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "VT",
          "sensitivity.ipp": 0.803156146179402,
          "sensitivity.maxent": 0.796511627906977,
          "sensitivity.logistic.regression": 0.8023,
          "sensitivity.knn": 0.7234,
          "sensitivity.classification.tree": 0.7500,
          "sensitivity.random.forest": 0.7575,
          "sensitivity.xgboost": 0.730897009966777
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "VT",
          "sensitivity.ipp": 0.797333333333333,
          "sensitivity.maxent": 0.782666666666667,
          "sensitivity.logistic.regression": 0.7760,
          "sensitivity.knn": 0.7687,
          "sensitivity.classification.tree": 0.7567,
          "sensitivity.random.forest": 0.7880,
          "sensitivity.xgboost": 0.764
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "VT",
          "sensitivity.ipp": 0.752427184466019,
          "sensitivity.maxent": 0.754045307443366,
          "sensitivity.logistic.regression": 0.7433,
          "sensitivity.knn": 0.7816,
          "sensitivity.classification.tree": 0.7740,
          "sensitivity.random.forest": 0.7740,
          "sensitivity.xgboost": 0.7859
        },
        {
          "common.name": "Ruddy Duck",
          "state": "VT",
          "sensitivity.ipp": 0.943396226415094,
          "sensitivity.maxent": 0.981132075471698,
          "sensitivity.logistic.regression": 0.9717,
          "sensitivity.knn": 0.9481,
          "sensitivity.classification.tree": 0.9481,
          "sensitivity.random.forest": 0.9906,
          "sensitivity.xgboost": 0.9481
        },
        {
          "common.name": "Sanderling",
          "state": "VT",
          "sensitivity.ipp": 0.947619047619048,
          "sensitivity.maxent": 1,
          "sensitivity.logistic.regression": 0.9562,
          "sensitivity.knn": 0.9124,
          "sensitivity.classification.tree": 0.9280,
          "sensitivity.random.forest": 0.9656,
          "sensitivity.xgboost": 0.9609
        },
        {
          "common.name": "Sandhill Crane",
          "state": "VT",
          "sensitivity.ipp": 0.91324200913242,
          "sensitivity.maxent": 0.977168949771689,
          "sensitivity.logistic.regression": 0.9832,
          "sensitivity.knn": 0.9912,
          "sensitivity.classification.tree": 0.9912,
          "sensitivity.random.forest": 0.9912,
          "sensitivity.xgboost": 0.9930
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "VT",
          "sensitivity.ipp": 0.870748299319728,
          "sensitivity.maxent": 0.861678004535147,
          "sensitivity.logistic.regression": 0.959183673469388,
          "sensitivity.knn": 0.952380952380952,
          "sensitivity.classification.tree": 0.936507936507937,
          "sensitivity.random.forest": 0.649730561970747,
          "sensitivity.xgboost": 0.952380952380952
        },
        {
          "common.name": "Wild Turkey",
          "state": "VT",
          "sensitivity.ipp": 0.709006928406467,
          "sensitivity.maxent": 0.821401077752117,
          "sensitivity.logistic.regression": 0.648960739030023,
          "sensitivity.knn": 0.659738260200154,
          "sensitivity.classification.tree": 0.696689761354888,
          "sensitivity.random.forest": 0.649730561970747,
          "sensitivity.xgboost": 0.68437259430331
        }
      ];

      jQuery(document).ready(function($){
        $('#DataTables_Sensitivity').DataTable({
          data: SensitivityResults,
          columns: [
            { data: function(row) { return row["common.name"]; } },
            { data: "state" },
            { data: function(row) { return row["sensitivity.ipp"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["sensitivity.maxent"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["sensitivity.logistic.regression"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["sensitivity.knn"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["sensitivity.classification.tree"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["sensitivity.random.forest"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["sensitivity.xgboost"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }}
          ],
          scrollX: true,
          scrollY: "600px",
          paging: false,
          searching: false
        });
      });
    </script>
  </div>

  <!-- Specificity Table -->
  <div id="Specificity_datatable" style="width: 100%; margin-top: 15px; display: none;">
    <div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item html-widget-static-bound" id="htmlwidget-specificity" style="width:100%;height:auto;">
      <table id="DataTables_Specificity" class="display dataTable no-footer table table-condensed" style="width:100%;">
        <thead>
          <tr>
            <th>common.name</th>
            <th>state</th>
            <th>specificity.ipp</th>
            <th>specificity.maxent</th>
            <th>specificity.logistic.regression</th>
            <th>specificity.knn</th>
            <th>specificity.classification.tree</th>
            <th>specificity.random.forest</th>
            <th>specificity.xgboost</th>
          </tr>
        </thead>
      </table>
    </div>
    <script>
      // Define the JSON dataset for Specificity (replace with your actual data)
      var SpecificityResults = [
        {
          "common.name": "Belted Kingfisher",
          "state": "CO",
          "specificity.ipp": 0.99184581171238,
          "specificity.maxent": 0.980726464047443,
          "specificity.logistic.regression": 0.989621942179392,
          "specificity.knn": 0.994069681245367,
          "specificity.classification.tree": 0.978502594514455,
          "specificity.random.forest": 0.989621942179392,
          "specificity.xgboost": 0.994069681245367
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "CO",
          "specificity.ipp": 0.982283464566929,
          "specificity.maxent": 0.984251968503937,
          "specificity.logistic.regression": 0.984251968503937,
          "specificity.knn": 0.990157480314961,
          "specificity.classification.tree": 0.954724409448819,
          "specificity.random.forest": 0.9827,
          "specificity.xgboost": 0.985236220472441
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "CO",
          "specificity.ipp": 0.993736017897092,
          "specificity.maxent": 0.978076062639821,
          "specificity.logistic.regression": 0.987919463087248,
          "specificity.knn": 0.993288590604027,
          "specificity.classification.tree": 0.982102908277405,
          "specificity.random.forest": 0.9582,
          "specificity.xgboost": 0.987472035794183
        },
        {
          "common.name": "Ruddy Duck",
          "state": "CO",
          "specificity.ipp": 0.98997995991984,
          "specificity.maxent": 0.963927855711423,
          "specificity.logistic.regression": 0.981963927855711,
          "specificity.knn": 0.98997995991984,
          "specificity.classification.tree": 0.977955911823647,
          "specificity.random.forest": 0.981963927855711,
          "specificity.xgboost": 0.98997995991984
        },
        {
          "common.name": "Sanderling",
          "state": "CO",
          "specificity.ipp": 1,
          "specificity.maxent": 1,
          "specificity.logistic.regression": 0.975,
          "specificity.knn": 0.955,
          "specificity.classification.tree": 0.83,
          "specificity.random.forest": 0.975,
          "specificity.xgboost": 0.955
        },
        {
          "common.name": "Sandhill Crane",
          "state": "CO",
          "specificity.ipp": 0.959550561797753,
          "specificity.maxent": 0.959550561797753,
          "specificity.logistic.regression": 0.98876404494382,
          "specificity.knn": 0.984269662921348,
          "specificity.classification.tree": 0.969662921348315,
          "specificity.random.forest": 0.99438202247191,
          "specificity.xgboost": 0.9888
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "CO",
          "specificity.ipp": 0.88212927756654,
          "specificity.maxent": 0.907984790874525,
          "specificity.logistic.regression": 0.959695817490494,
          "specificity.knn": 0.968821292775665,
          "specificity.classification.tree": 0.944486692015209,
          "specificity.random.forest": 0.965779467680608,
          "specificity.xgboost": 0.9597
        },
        {
          "common.name": "Wild Turkey",
          "state": "CO",
          "specificity.ipp": 0.963446475195822,
          "specificity.maxent": 0.963446475195822,
          "specificity.logistic.regression": 0.9791,
          "specificity.knn": 0.9843,
          "specificity.classification.tree": 0.969321148825065,
          "specificity.random.forest": 0.9883,
          "specificity.xgboost": 0.9843
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "NC",
          "specificity.ipp": 0.5483,
          "specificity.maxent": 0.8156,
          "specificity.logistic.regression": 0.9139,
          "specificity.knn": 0.9352,
          "specificity.classification.tree": 0.9008,
          "specificity.random.forest": 0.9557,
          "specificity.xgboost": 0.9615
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "NC",
          "specificity.ipp": 0.8820,
          "specificity.maxent": 0.9078,
          "specificity.logistic.regression": 0.9333,
          "specificity.knn": 0.9455,
          "specificity.classification.tree": 0.9374,
          "specificity.random.forest": 0.9724,
          "specificity.xgboost": 0.9732
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "NC",
          "specificity.ipp": 0.8646,
          "specificity.maxent": 0.9491,
          "specificity.logistic.regression": 0.9641,
          "specificity.knn": 0.9399,
          "specificity.classification.tree": 0.9255,
          "specificity.random.forest": 0.9555,
          "specificity.xgboost": 0.9586
        },
        {
          "common.name": "Ruddy Duck",
          "state": "NC",
          "specificity.ipp": 0.8614,
          "specificity.maxent": 0.9327,
          "specificity.logistic.regression": 0.9562,
          "specificity.knn": 0.9270,
          "specificity.classification.tree": 0.9143,
          "specificity.random.forest": 0.9656,
          "specificity.xgboost": 0.9609
        },
        {
          "common.name": "Sanderling",
          "state": "NC",
          "specificity.ipp": 0.9011,
          "specificity.maxent": 0.9772,
          "specificity.logistic.regression": 0.9734,
          "specificity.knn": 0.9924,
          "specificity.classification.tree": 0.9506,
          "specificity.random.forest": 0.9962,
          "specificity.xgboost": 0.9582
        },
        {
          "common.name": "Sandhill Crane",
          "state": "NC",
          "specificity.ipp": 0.8652,
          "specificity.maxent": 0.9435,
          "specificity.logistic.regression": 0.9087,
          "specificity.knn": 0.9565,
          "specificity.classification.tree": 0.8130,
          "specificity.random.forest": 0.9217,
          "specificity.xgboost": 0.8609
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "NC",
          "specificity.ipp": 0.8731,
          "specificity.maxent": 0.8840,
          "specificity.logistic.regression": 0.9378,
          "specificity.knn": 0.9351,
          "specificity.classification.tree": 0.8936,
          "specificity.random.forest": 0.9406,
          "specificity.xgboost": 0.9530
        },
        {
          "common.name": "Wild Turkey",
          "state": "NC",
          "specificity.ipp": 0.8653,
          "specificity.maxent": 0.8703,
          "specificity.logistic.regression": 0.9348,
          "specificity.knn": 0.9219,
          "specificity.classification.tree": 0.9083,
          "specificity.random.forest": 0.9341,
          "specificity.xgboost": 0.9391
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "OR",
          "specificity.ipp": 0.9868,
          "specificity.maxent": 0.9810,
          "specificity.logistic.regression": 0.9983,
          "specificity.knn": 0.9916,
          "specificity.classification.tree": 0.9753,
          "specificity.random.forest": 0.9886,
          "specificity.xgboost": 0.9916
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "OR",
          "specificity.ipp": 0.9478,
          "specificity.maxent": 0.9771,
          "specificity.logistic.regression": 0.9941,
          "specificity.knn": 0.9959,
          "specificity.classification.tree": 0.9925,
          "specificity.random.forest": 0.9969,
          "specificity.xgboost": 0.9963
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "OR",
          "specificity.ipp": 0.9874,
          "specificity.maxent": 0.9862,
          "specificity.logistic.regression": 0.9976,
          "specificity.knn": 0.9968,
          "specificity.classification.tree": 0.9910,
          "specificity.random.forest": 0.9978,
          "specificity.xgboost": 0.9968
        },
        {
          "common.name": "Ruddy Duck",
          "state": "OR",
          "specificity.ipp": 0.9576,
          "specificity.maxent": 0.9814,
          "specificity.logistic.regression": 0.9832,
          "specificity.knn": 0.9947,
          "specificity.classification.tree": 0.9920,
          "specificity.random.forest": 0.9947,
          "specificity.xgboost": 0.9947
        },
        {
          "common.name": "Sanderling",
          "state": "OR",
          "specificity.ipp": 0.8724,
          "specificity.maxent": 0.9959,
          "specificity.logistic.regression": 0.9835,
          "specificity.knn": 0.9918,
          "specificity.classification.tree": 0.9300,
          "specificity.random.forest": 0.9835,
          "specificity.xgboost": 0.9959
        },
        {
          "common.name": "Sandhill Crane",
          "state": "OR",
          "specificity.ipp": 0.9723,
          "specificity.maxent": 0.9650,
          "specificity.logistic.regression": 0.9920,
          "specificity.knn": 0.9984,
          "specificity.classification.tree": 0.9857,
          "specificity.random.forest": 0.9904,
          "specificity.xgboost": 0.9904
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "OR",
          "specificity.ipp": 0.8594,
          "specificity.maxent": 0.9645,
          "specificity.logistic.regression": 0.9884,
          "specificity.knn": 0.9839,
          "specificity.classification.tree": 0.9498,
          "specificity.random.forest": 0.9848,
          "specificity.xgboost": 0.9949
        },
        {
          "common.name": "Wild Turkey",
          "state": "OR",
          "specificity.ipp": 0.91220556745182,
          "specificity.maxent": 0.9914,
          "specificity.logistic.regression": 0.9879,
          "specificity.knn": 0.9971,
          "specificity.classification.tree": 0.9793,
          "specificity.random.forest": 0.9914,
          "specificity.xgboost": 0.9950
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "VT",
          "specificity.ipp": 0.803156146179402,
          "specificity.maxent": 0.796511627906977,
          "specificity.logistic.regression": 0.8023,
          "specificity.knn": 0.7234,
          "specificity.classification.tree": 0.7500,
          "specificity.random.forest": 0.7575,
          "specificity.xgboost": 0.730897009966777
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "VT",
          "specificity.ipp": 0.797333333333333,
          "specificity.maxent": 0.782666666666667,
          "specificity.logistic.regression": 0.7760,
          "specificity.knn": 0.7687,
          "specificity.classification.tree": 0.7567,
          "specificity.random.forest": 0.7880,
          "specificity.xgboost": 0.764
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "VT",
          "specificity.ipp": 0.752427184466019,
          "specificity.maxent": 0.754045307443366,
          "specificity.logistic.regression": 0.7433,
          "specificity.knn": 0.7816,
          "specificity.classification.tree": 0.7740,
          "specificity.random.forest": 0.7740,
          "specificity.xgboost": 0.7859
        },
        {
          "common.name": "Ruddy Duck",
          "state": "VT",
          "specificity.ipp": 0.943396226415094,
          "specificity.maxent": 0.981132075471698,
          "specificity.logistic.regression": 0.9717,
          "specificity.knn": 0.9481,
          "specificity.classification.tree": 0.9481,
          "specificity.random.forest": 0.9906,
          "specificity.xgboost": 0.9481
        },
        {
          "common.name": "Sanderling",
          "state": "VT",
          "specificity.ipp": 0.947619047619048,
          "specificity.maxent": 1,
          "specificity.logistic.regression": 0.9562,
          "specificity.knn": 0.9124,
          "specificity.classification.tree": 0.9280,
          "specificity.random.forest": 0.9656,
          "specificity.xgboost": 0.9609
        },
        {
          "common.name": "Sandhill Crane",
          "state": "VT",
          "specificity.ipp": 0.91324200913242,
          "specificity.maxent": 0.977168949771689,
          "specificity.logistic.regression": 0.9832,
          "specificity.knn": 0.9912,
          "specificity.classification.tree": 0.9912,
          "specificity.random.forest": 0.9912,
          "specificity.xgboost": 0.9930
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "VT",
          "specificity.ipp": 0.870748299319728,
          "specificity.maxent": 0.861678004535147,
          "specificity.logistic.regression": 0.959183673469388,
          "specificity.knn": 0.952380952380952,
          "specificity.classification.tree": 0.936507936507937,
          "specificity.random.forest": 0.649730561970747,
          "specificity.xgboost": 0.952380952380952
        },
        {
          "common.name": "Wild Turkey",
          "state": "VT",
          "specificity.ipp": 0.709006928406467,
          "specificity.maxent": 0.821401077752117,
          "specificity.logistic.regression": 0.648960739030023,
          "specificity.knn": 0.659738260200154,
          "specificity.classification.tree": 0.696689761354888,
          "specificity.random.forest": 0.649730561970747,
          "specificity.xgboost": 0.68437259430331
        }
      ];

      jQuery(document).ready(function($){
        $('#DataTables_Specificity').DataTable({
          data: SpecificityResults,
          columns: [
            { data: function(row) { return row["common.name"]; } },
            { data: "state" },
            { data: function(row) { return row["specificity.ipp"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["specificity.maxent"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["specificity.logistic.regression"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["specificity.knn"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["specificity.classification.tree"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["specificity.random.forest"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["specificity.xgboost"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }}
          ],
          scrollX: true,
          scrollY: "600px",
          paging: false,
          searching: false
        });
      });
    </script>
  </div>

  <!-- F1 Table -->
  <div id="F1_datatable" style="width: 100%; margin-top: 15px; display: none;">
    <div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item html-widget-static-bound" id="htmlwidget-f1" style="width:100%;height:auto;">
      <table id="DataTables_F1" class="display dataTable no-footer table table-condensed" style="width:100%;">
        <thead>
          <tr>
            <th>common.name</th>
            <th>state</th>
            <th>f1.ipp</th>
            <th>f1.maxent</th>
            <th>f1.logistic.regression</th>
            <th>f1.knn</th>
            <th>f1.classification.tree</th>
            <th>f1.random.forest</th>
            <th>f1.xgboost</th>
          </tr>
        </thead>
      </table>
    </div>
    <script>
      // Define the JSON dataset for F1 (replace with your actual data)
      var F1Results = [
        {
          "common.name": "Belted Kingfisher",
          "state": "CO",
          "f1.ipp": 0.9797,
          "f1.maxent": 0.9619,
          "f1.logistic.regression": 0.9862,
          "f1.knn": 0.9847,
          "f1.classification.tree": 0.9670,
          "f1.random.forest": 0.9881,
          "f1.xgboost": 0.9881
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "CO",
          "f1.ipp": 0.9513,
          "f1.maxent": 0.9612,
          "f1.logistic.regression": 0.9726,
          "f1.knn": 0.9755,
          "f1.classification.tree": 0.9563,
          "f1.random.forest": 0.9827,
          "f1.xgboost": 0.9827
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "CO",
          "f1.ipp": 0.9738,
          "f1.maxent": 0.9566,
          "f1.logistic.regression": 0.9834,
          "f1.knn": 0.9880,
          "f1.classification.tree": 0.9822,
          "f1.random.forest": 0.9582,
          "f1.xgboost": 0.9582
        },
        {
          "common.name": "Ruddy Duck",
          "state": "CO",
          "f1.ipp": 0.9693,
          "f1.maxent": 0.9618,
          "f1.logistic.regression": 0.9780,
          "f1.knn": 0.9880,
          "f1.classification.tree": 0.9785,
          "f1.random.forest": 0.9840,
          "f1.xgboost": 0.9840
        },
        {
          "common.name": "Sanderling",
          "state": "CO",
          "f1.ipp": 0.8696,
          "f1.maxent": 0.9696,
          "f1.logistic.regression": 0.9609,
          "f1.knn": 0.9478,
          "f1.classification.tree": 0.9008,
          "f1.random.forest": 0.9667,
          "f1.xgboost": 0.9667
        },
        {
          "common.name": "Sandhill Crane",
          "state": "CO",
          "f1.ipp": 0.9581,
          "f1.maxent": 0.9596,
          "f1.logistic.regression": 0.9888,
          "f1.knn": 0.9843,
          "f1.classification.tree": 0.9557,
          "f1.random.forest": 0.9910,
          "f1.xgboost": 0.9910
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "CO",
          "f1.ipp": 0.9878,
          "f1.maxent": 0.9080,
          "f1.logistic.regression": 0.9597,
          "f1.knn": 0.9688,
          "f1.classification.tree": 0.9785,
          "f1.random.forest": 0.9881,
          "f1.xgboost": 0.9881
        },
        {
          "common.name": "Wild Turkey",
          "state": "CO",
          "f1.ipp": 0.9634,
          "f1.maxent": 0.9634,
          "f1.logistic.regression": 0.9791,
          "f1.knn": 0.9843,
          "f1.classification.tree": 0.9693,
          "f1.random.forest": 0.9883,
          "f1.xgboost": 0.9883
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "NC",
          "f1.ipp": 0.5483,
          "f1.maxent": 0.8156,
          "f1.logistic.regression": 0.9139,
          "f1.knn": 0.9352,
          "f1.classification.tree": 0.9008,
          "f1.random.forest": 0.9557,
          "f1.xgboost": 0.9615
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "NC",
          "f1.ipp": 0.7920,
          "f1.maxent": 0.9000,
          "f1.logistic.regression": 0.9333,
          "f1.knn": 0.9455,
          "f1.classification.tree": 0.9374,
          "f1.random.forest": 0.9724,
          "f1.xgboost": 0.9732
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "NC",
          "f1.ipp": 0.9937,
          "f1.maxent": 0.9491,
          "f1.logistic.regression": 0.9641,
          "f1.knn": 0.9399,
          "f1.classification.tree": 0.9255,
          "f1.random.forest": 0.9555,
          "f1.xgboost": 0.9586
        },
        {
          "common.name": "Ruddy Duck",
          "state": "NC",
          "f1.ipp": 0.9499,
          "f1.maxent": 0.9599,
          "f1.logistic.regression": 0.9739,
          "f1.knn": 0.9860,
          "f1.classification.tree": 0.9499,
          "f1.random.forest": 0.9880,
          "f1.xgboost": 0.9840
        },
        {
          "common.name": "Sanderling",
          "state": "NC",
          "f1.ipp": 0.000,
          "f1.maxent": 0.7667,
          "f1.logistic.regression": 0.8667,
          "f1.knn": 0.9000,
          "f1.classification.tree": 0.9667,
          "f1.random.forest": 0.9667,
          "f1.xgboost": 0.9667
        },
        {
          "common.name": "Sandhill Crane",
          "state": "NC",
          "f1.ipp": 0.9258,
          "f1.maxent": 0.9618,
          "f1.logistic.regression": 0.9843,
          "f1.knn": 0.9730,
          "f1.classification.tree": 0.9708,
          "f1.random.forest": 0.9910,
          "f1.xgboost": 0.9865
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "NC",
          "f1.ipp": 0.7766,
          "f1.maxent": 0.8830,
          "f1.logistic.regression": 0.9498,
          "f1.knn": 0.9650,
          "f1.classification.tree": 0.9407,
          "f1.random.forest": 0.9726,
          "f1.xgboost": 0.9880
        },
        {
          "common.name": "Wild Turkey",
          "state": "NC",
          "f1.ipp": 0.9634,
          "f1.maxent": 0.9634,
          "f1.logistic.regression": 0.9791,
          "f1.knn": 0.9843,
          "f1.classification.tree": 0.9693,
          "f1.random.forest": 0.9883,
          "f1.xgboost": 0.9876
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "OR",
          "f1.ipp": 0.9868,
          "f1.maxent": 0.9810,
          "f1.logistic.regression": 0.9983,
          "f1.knn": 0.9916,
          "f1.classification.tree": 0.9753,
          "f1.random.forest": 0.9886,
          "f1.xgboost": 0.9916
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "OR",
          "f1.ipp": 0.9478,
          "f1.maxent": 0.9771,
          "f1.logistic.regression": 0.9941,
          "f1.knn": 0.9959,
          "f1.classification.tree": 0.9925,
          "f1.random.forest": 0.9969,
          "f1.xgboost": 0.9963
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "OR",
          "f1.ipp": 0.9874,
          "f1.maxent": 0.9862,
          "f1.logistic.regression": 0.9976,
          "f1.knn": 0.9968,
          "f1.classification.tree": 0.9910,
          "f1.random.forest": 0.9978,
          "f1.xgboost": 0.9968
        },
        {
          "common.name": "Ruddy Duck",
          "state": "OR",
          "f1.ipp": 0.9576,
          "f1.maxent": 0.9814,
          "f1.logistic.regression": 0.9832,
          "f1.knn": 0.9947,
          "f1.classification.tree": 0.9920,
          "f1.random.forest": 0.9947,
          "f1.xgboost": 0.9947
        },
        {
          "common.name": "Sanderling",
          "state": "OR",
          "f1.ipp": 0.8724,
          "f1.maxent": 0.9959,
          "f1.logistic.regression": 0.9835,
          "f1.knn": 0.9918,
          "f1.classification.tree": 0.9300,
          "f1.random.forest": 0.9835,
          "f1.xgboost": 0.9959
        },
        {
          "common.name": "Sandhill Crane",
          "state": "OR",
          "f1.ipp": 0.9723,
          "f1.maxent": 0.9650,
          "f1.logistic.regression": 0.9920,
          "f1.knn": 0.9984,
          "f1.classification.tree": 0.9857,
          "f1.random.forest": 0.9904,
          "f1.xgboost": 0.9904
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "OR",
          "f1.ipp": 0.8594,
          "f1.maxent": 0.9645,
          "f1.logistic.regression": 0.9884,
          "f1.knn": 0.9839,
          "f1.classification.tree": 0.9498,
          "f1.random.forest": 0.9848,
          "f1.xgboost": 0.9949
        },
        {
          "common.name": "Wild Turkey",
          "state": "OR",
          "f1.ipp": 0.9122,
          "f1.maxent": 0.9914,
          "f1.logistic.regression": 0.9879,
          "f1.knn": 0.9971,
          "f1.classification.tree": 0.9793,
          "f1.random.forest": 0.9914,
          "f1.xgboost": 0.9950
        },
        {
          "common.name": "Belted Kingfisher",
          "state": "VT",
          "f1.ipp": 0.8032,
          "f1.maxent": 0.7965,
          "f1.logistic.regression": 0.8023,
          "f1.knn": 0.7234,
          "f1.classification.tree": 0.7500,
          "f1.random.forest": 0.7575,
          "f1.xgboost": 0.730897009966777
        },
        {
          "common.name": "Cedar Waxwing",
          "state": "VT",
          "f1.ipp": 0.797333333333333,
          "f1.maxent": 0.782666666666667,
          "f1.logistic.regression": 0.7760,
          "f1.knn": 0.7687,
          "f1.classification.tree": 0.7567,
          "f1.random.forest": 0.7880,
          "f1.xgboost": 0.764
        },
        {
          "common.name": "Downy Woodpecker",
          "state": "VT",
          "f1.ipp": 0.752427184466019,
          "f1.maxent": 0.754045307443366,
          "f1.logistic.regression": 0.7433,
          "f1.knn": 0.7816,
          "f1.classification.tree": 0.7740,
          "f1.random.forest": 0.7740,
          "f1.xgboost": 0.7859
        },
        {
          "common.name": "Ruddy Duck",
          "state": "VT",
          "f1.ipp": 0.943396226415094,
          "f1.maxent": 0.981132075471698,
          "f1.logistic.regression": 0.9717,
          "f1.knn": 0.9481,
          "f1.classification.tree": 0.9481,
          "f1.random.forest": 0.9906,
          "f1.xgboost": 0.9481
        },
        {
          "common.name": "Sanderling",
          "state": "VT",
          "f1.ipp": 0.947619047619048,
          "f1.maxent": 1,
          "f1.logistic.regression": 0.9562,
          "f1.knn": 0.9124,
          "f1.classification.tree": 0.9280,
          "f1.random.forest": 0.9656,
          "f1.xgboost": 0.9609
        },
        {
          "common.name": "Sandhill Crane",
          "state": "VT",
          "f1.ipp": 0.91324200913242,
          "f1.maxent": 0.977168949771689,
          "f1.logistic.regression": 0.9832,
          "f1.knn": 0.9912,
          "f1.classification.tree": 0.9912,
          "f1.random.forest": 0.9912,
          "f1.xgboost": 0.9930
        },
        {
          "common.name": "Sharp-shinned Hawk",
          "state": "VT",
          "f1.ipp": 0.870748299319728,
          "f1.maxent": 0.861678004535147,
          "f1.logistic.regression": 0.959183673469388,
          "f1.knn": 0.952380952380952,
          "f1.classification.tree": 0.936507936507937,
          "f1.random.forest": 0.649730561970747,
          "f1.xgboost": 0.952380952380952
        },
        {
          "common.name": "Wild Turkey",
          "state": "VT",
          "f1.ipp": 0.9122,
          "f1.maxent": 0.9914,
          "f1.logistic.regression": 0.9879,
          "f1.knn": 0.9971,
          "f1.classification.tree": 0.9793,
          "f1.random.forest": 0.9914,
          "f1.xgboost": 0.9950
        }
      ];

      jQuery(document).ready(function($){
        $('#DataTables_F1').DataTable({
          data: F1Results,
          columns: [
            { data: function(row) { return row["common.name"]; } },
            { data: "state" },
            { data: function(row) { return row["f1.ipp"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["f1.maxent"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["f1.logistic.regression"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["f1.knn"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["f1.classification.tree"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["f1.random.forest"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }},
            { data: function(row) { return row["f1.xgboost"]; }, render: function(data, type, row, meta) {
                return type === 'display' && data !== null && !isNaN(data) ? parseFloat(data).toFixed(4) : data;
            }}
          ],
          scrollX: true,
          scrollY: "600px",
          paging: false,
          searching: false
        });
      });
    </script>
  </div>

</div>


### Overall Performance by Model Type

```

model.scores %>% setDT()
model.scores[, `:=` (
    model.type = factor(model.type,
                        levels=unique(model.scores$model.type)),
    common.name = factor(common.name,
                         levels=unique(model.scores$common.name)),
    state = factor(state, 
                   levels=unique(model.scores$state))
  )]

plt.data <- model.scores %>%
  melt(.,
       id.vars = c("common.name", "state", "model.type"), 
       measure.vars = c("Accuracy", "Sensitivity", 
                        "Specificity", "F1"),
       variable.name = "metric", 
       value.name = "metric.value")

ggplot(plt.data, 
         aes(x = metric, y = metric.value, fill = model.type)) +
    geom_boxplot(outlier.shape = "o") +
    stat_summary(aes(group = model.type), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75)) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray"),
          plot.background = element_rect(fill = "white"),
          strip.placement = "inside") + 
    coord_flip() +
    facet_wrap(~metric, ncol=2, scales="free",
               strip.position="top")

```

<img src="{{ site.baseurl }}/assets/plots/sdm-9-1.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 850px; height: auto;">

### Performance by Model, State

```

ggplot(plt.data, 
         aes(x = metric, y = metric.value, fill = model.type)) +
    geom_boxplot(outlier.shape = "o") +
    stat_summary(aes(group = model.type), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75)) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray"),
          plot.background = element_rect(fill = "white"),
          strip.placement = "inside") + 
    coord_flip() +
    facet_wrap(~metric + state, ncol=2, scales="free",
               strip.position="top")
```

<img src="{{ site.baseurl }}/assets/plots/sdm-9-2.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 850px; height: auto;">

### Performance by Model, Species

```
ggplot(plt.data, 
         aes(x = metric, y = metric.value, fill = model.type)) +
    geom_boxplot(outlier.shape = "o") +
    stat_summary(aes(group = model.type), fun = mean, fill="darkred",
                 geom = "point", shape = 21, size = 2, color = "black",
                 position = position_dodge(width = 0.75)) +
    scale_fill_brewer(palette = "Dark2") + 
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = "gray"),
          plot.background = element_rect(fill = "white"),
          strip.placement = "inside") + 
    coord_flip() +
    facet_wrap(~metric + common.name, ncol=2, scales="free",
               strip.position="top")
```

<img src="{{ site.baseurl }}/assets/plots/sdm-9-3.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 850px; height: auto;">

### Discussion

Overall, the random forest model exhibits the highest average accuracy (94.64%) among all the
models, which indicates its robustness in predicting species distribution across various
datasets. This model also demonstrates high sensitivity (96.69%), suggesting its effectiveness
in correctly identifying true positives. In contrast, the IPP model shows the lowest overall
accuracy (89.53%), coupled with significantly lower sensitivity (69.92%), indicating a 
tendency to miss true positive cases. 

A more granular analysis reveals variability in model performance across different bird 
species and states. Notably, there are several instances of perfect accuracy
(100%). These appear to occur more frequently for bird species exhibiting a more
defined distribution (e.g., primarily coastal birds), in smaller sample areas. Conversely, 
it appears that the models fit for more broadly distributed species, (such as the Downy
Woodpecker), tend to have poorer performance when using the baseline models, while the ML 
models performed significantly better.

In summary, while some models like random forest and KNN generally outperform others in 
terms of accuracy and sensitivity, the effectiveness of each model can vary considerably
depending on the bird species and the geographical context. Selecting the most
suitable model for a given SDM task requires careful consideration of the specific
characteristics of the dataset. 

## Model Performance Comparisons

The analysis of model accuracy below initially uses Analysis of Variance (ANOVA), aimed at 
identifying significant differences across accuracy by model type. However, as shown later, 
the Shapiro-Wilk test indicates non-normality in ANOVA residuals, necessitating a shift to
non-parametric methods. These methods, including the Kruskal-Wallis test and subsequent Dunn’s
post-hoc tests, offer a more suitable approach for this data set. This part of the study is
crucial for evaluating the performance differences among models, clearly identifying those with
statistically significant variations in accuracy.

```
anova.res <- aov(Accuracy ~ model.type, data = model.scores)
summary.aov(anova.res) 
```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
##              Df Sum Sq  Mean Sq F value Pr(>F)
## model.type    6 0.0606 0.010101   1.593   0.15
## Residuals   217 1.3756 0.006339
</code></pre></div><br>

### Confirming ANOVA Assumptions

The Shapiro-Wilk test can be used to confirm whether the residuals of an ANOVA model are
normally distributed. The null hypothesis is that the data is normally distributed.
If $p < alpha = 0.05$, then the null hypothesis is rejected, and it can be concluded that the
residuals of the ANOVA are not normally distributed

```
shapiro.test(residuals(anova.res))
```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
## 
##  Shapiro-Wilk normality test
## 
## data:  residuals(anova.res)
## W = 0.76373, p-value < 2.2e-16
</code></pre></div><br>

### Non-Parametric Alternative to ANOVA

The Kruskal-Wallis Rank Sum Test is a non-parametric alternative to one-way ANOVA. 
It's used to compare more than two groups and does not assume a normal distribution.
The null hypothesis for this test is that all groups have the same median, with a
threshold of 0.05.

```
kruskal.test(Accuracy ~ model.type, data = model.scores)
```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
## 
##  Kruskal-Wallis rank sum test
## 
## data:  Accuracy by model.type
## Kruskal-Wallis chi-squared = 26.037, df = 6, p-value = 0.0002191
</code></pre></div><br>

The Kruskal-Wallis Chi-squared value (the test statistic for the Kruskal-Wallis test) is 
26.037. This is analogous to the F-statistic in ANOVA but based on ranks rather than means.
The Degrees of Freedom (df) is 6, which corresponds to the number of groups (model types)
minus one. It's an important parameter in determining the statistical significance.
The p-value is 0.0002191. This indicates the probability of obtaining test results at least 
as extreme as the ones observed, under the assumption that the null hypothesis is true.
Because the p-value is much less than the threshold of 0.05, the null hypothesis of the
Kruskal-Wallis test is rejected. In other words, there is statistically significant evidence 
to suggest that at least one model type has a different median accuracy compared to the others.

### Non-Parametric Pairwise Comparisons

To identify which specific model types differ from each other, conduct post-hoc 
pairwise comparisons using the non-parametric Dunn's test.
    
```

# Dunn's test for pairwise comparisons
dunn.res <- dunn.test::dunn.test(model.scores$Accuracy, 
                                 model.scores$model.type, 
                                 method="bonferroni")

```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
##   Kruskal-Wallis rank sum test
## 
## data: x and group
## Kruskal-Wallis chi-squared = 26.0371, df = 6, p-value = 0
## 
## 
##                            Comparison of x by group                            
##                                  (Bonferroni)                                  
## Col Mean-|
## Row Mean |   classifi        ipp        knn   logistic     maxent   random f
## ---------+------------------------------------------------------------------
##      ipp |   1.617309
##          |     1.0000
##          |
##      knn |  -2.060936  -3.678245
##          |     0.4127    0.0025*
##          |
## logistic |  -1.513153  -3.130462   0.547782
##          |     1.0000    0.0183*     1.0000
##          |
##   maxent |  -0.394442  -2.011751   1.666493   1.118711
##          |     1.0000     0.4646     1.0000     1.0000
##          |
## random f |  -2.447663  -4.064972  -0.386726  -0.934509  -2.053220
##          |     0.1510    0.0005*     1.0000     1.0000     0.4205
##          |
##  xgboost |  -2.147732  -3.765041  -0.086796  -0.634579  -1.753290   0.299930
##          |     0.3332    0.0017*     1.0000     1.0000     0.8353     1.0000
## 
## alpha = 0.05
## Reject Ho if p <= alpha/2
</code></pre></div><br>

Each cell in the table represents the comparison between the model types in the
corresponding row and column. For instance, the first cell under "ipp" compares the
Classification Tree models with IPP. The values in each cell are the test statistic for 
Dunn's test. These values indicate the degree of difference between the pair of groups. 
The p-values are given below each test statistic. These indicate the probability of observing 
a test statistic as extreme as, or more extreme than, the one observed, assuming the null
hypothesis of no difference is true.

#### Significant Differences

Comparisons with a p-value less than 0.05/2 (adjusted for the Bonferroni correction) are
considered statistically significant. For example, the comparison between KNN and IPP shows
a p-value of 0.0025, which is significant. This suggests that there is a statistically
significant difference in the accuracy between these two model types.

#### No Significant Differences

Comparisons with a p-value greater than 0.05/2 are not considered statistically significant. 
For instance, the comparison between MaxEnt and IPP shows a p-value of 0.4646, which is not
significant.

### Effect Size 

To understand the magnitude of differences, not just the statistical significance, calculate
Cliff's Delta for each pair of model types. This metric quantifies the difference in the
likelihood that a randomly selected case from one group (model-type 1) scores higher than a
randomly chosen case from another group (model-type 2), and vice versa. Essentially, Cliff's 
Delta is a measure of the 'success rate difference' between two groups, providing insight into
the practical significance of their performance disparity.
    
```

pairs <- expand.grid(
  treatment=unique(model.scores$model.type),
  control=unique(model.scores$model.type)
) %>% 
  filter(treatment != control)

eff.size <- purrr::map_df(1:nrow(pairs), ~{
  treatment <- pairs[.x,]$treatment
  control <- pairs[.x,]$control
  t.acc <- model.scores$Accuracy[model.scores$model.type == treatment]
  c.acc <- model.scores$Accuracy[model.scores$model.type == control]
  cd <- effsize::cliff.delta(t.acc, c.acc)
  data.table(
    treatment=treatment,
    control=control,
    delta=cd$estimate,
    lower=cd$conf.int[[1]],
    upper=cd$conf.int[[2]],
    variance=cd$var,
    magnitude=cd$magnitude
  )
})

eff.size
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script data-pagedtable-source="" type="application/json">
{"columns":[{"label":["treatment"],"name":[1],"type":["fct"],"align":["left"]},{"label":["control"],"name":[2],"type":["fct"],"align":["left"]},{"label":["delta"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["lower"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["upper"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["variance"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["magnitude"],"name":[7],"type":["ord"],"align":["right"]}],"data":[{"1":"maxent","2":"ipp","3":"0.33007812","4":"0.04096127","5":"0.56819123","6":"0.01867086","7":"medium"},{"1":"logistic regression","2":"ipp","3":"0.49023438","4":"0.20577262","5":"0.69830495","6":"0.01606384","7":"large"},{"1":"knn","2":"ipp","3":"0.51953125","4":"0.23551269","5":"0.72179254","6":"0.01560406","7":"large"},{"1":"classification tree","2":"ipp","3":"0.26171875","4":"-0.03301375","5":"0.51457381","6":"0.02027220","7":"small"},{"1":"random forest","2":"ipp","3":"0.55664062","4":"0.27539691","5":"0.75011165","6":"0.01478261","7":"large"},{"1":"xgboost","2":"ipp","3":"0.48437500","4":"0.19038175","5":"0.69863860","6":"0.01717800","7":"large"},{"1":"ipp","2":"maxent","3":"-0.33007812","4":"-0.56819123","5":"-0.04096127","6":"0.01867086","7":"medium"},{"1":"logistic regression","2":"maxent","3":"0.19335938","4":"-0.09908924","5":"0.45506117","6":"0.02080763","7":"small"},{"1":"knn","2":"maxent","3":"0.25097656","4":"-0.04626317","5":"0.50738787","6":"0.02076231","7":"small"},{"1":"classification tree","2":"maxent","3":"-0.05566406","4":"-0.33281574","5":"0.23033962","6":"0.02155054","7":"negligible"},{"1":"random forest","2":"maxent","3":"0.30664062","4":"0.01172445","5":"0.55247822","6":"0.01972233","7":"small"},{"1":"xgboost","2":"maxent","3":"0.24023438","4":"-0.05684535","5":"0.49823258","6":"0.02087958","7":"small"},{"1":"ipp","2":"logistic regression","3":"-0.49023438","4":"-0.69830495","5":"-0.20577262","6":"0.01606384","7":"large"},{"1":"maxent","2":"logistic regression","3":"-0.19335938","4":"-0.45506117","5":"0.09908924","6":"0.02080763","7":"small"},{"1":"knn","2":"logistic regression","3":"0.10546875","4":"-0.18520119","5":"0.37916851","6":"0.02165143","7":"negligible"},{"1":"classification tree","2":"logistic regression","3":"-0.24414062","4":"-0.50032584","5":"0.05135053","6":"0.02060312","7":"small"},{"1":"random forest","2":"logistic regression","3":"0.17382812","4":"-0.11672552","5":"0.43697261","6":"0.02077171","7":"small"},{"1":"xgboost","2":"logistic regression","3":"0.12109375","4":"-0.17010680","5":"0.39284387","6":"0.02153317","7":"negligible"},{"1":"ipp","2":"knn","3":"-0.51953125","4":"-0.72179254","5":"-0.23551269","6":"0.01560406","7":"large"},{"1":"maxent","2":"knn","3":"-0.25097656","4":"-0.50738787","5":"0.04626317","6":"0.02076231","7":"small"},{"1":"logistic regression","2":"knn","3":"-0.10546875","4":"-0.37916851","5":"0.18520119","6":"0.02165143","7":"negligible"},{"1":"classification tree","2":"knn","3":"-0.30664062","4":"-0.55068718","5":"-0.01429824","6":"0.01938089","7":"small"},{"1":"random forest","2":"knn","3":"0.06445312","4":"-0.22033874","5":"0.33911952","6":"0.02124456","7":"negligible"},{"1":"xgboost","2":"knn","3":"0.03613281","4":"-0.24707251","5":"0.31365102","6":"0.02134899","7":"negligible"},{"1":"ipp","2":"classification tree","3":"-0.26171875","4":"-0.51457381","5":"0.03301375","6":"0.02027220","7":"small"},{"1":"maxent","2":"classification tree","3":"0.05566406","4":"-0.23033962","5":"0.33281574","6":"0.02155054","7":"negligible"},{"1":"logistic regression","2":"classification tree","3":"0.24414062","4":"-0.05135053","5":"0.50032584","6":"0.02060312","7":"small"},{"1":"knn","2":"classification tree","3":"0.30664062","4":"0.01429824","5":"0.55068718","6":"0.01938089","7":"small"},{"1":"random forest","2":"classification tree","3":"0.35156250","4":"0.06084953","5":"0.58729462","6":"0.01860505","7":"medium"},{"1":"xgboost","2":"classification tree","3":"0.30859375","4":"0.01475470","5":"0.55336949","6":"0.01955413","7":"small"},{"1":"ipp","2":"random forest","3":"-0.55664062","4":"-0.75011165","5":"-0.27539691","6":"0.01478261","7":"large"},{"1":"maxent","2":"random forest","3":"-0.30664062","4":"-0.55247822","5":"-0.01172445","6":"0.01972233","7":"small"},{"1":"logistic regression","2":"random forest","3":"-0.17382812","4":"-0.43697261","5":"0.11672552","6":"0.02077171","7":"small"},{"1":"knn","2":"random forest","3":"-0.06445312","4":"-0.33911952","5":"0.22033874","6":"0.02124456","7":"negligible"},{"1":"classification tree","2":"random forest","3":"-0.35156250","4":"-0.58729462","5":"-0.06084953","6":"0.01860505","7":"medium"},{"1":"xgboost","2":"random forest","3":"-0.02050781","4":"-0.29889866","5":"0.26109989","6":"0.02128913","7":"negligible"},{"1":"ipp","2":"xgboost","3":"-0.48437500","4":"-0.69863860","5":"-0.19038175","6":"0.01717800","7":"large"},{"1":"maxent","2":"xgboost","3":"-0.24023438","4":"-0.49823258","5":"0.05684535","6":"0.02087958","7":"small"},{"1":"logistic regression","2":"xgboost","3":"-0.12109375","4":"-0.39284387","5":"0.17010680","6":"0.02153317","7":"negligible"},{"1":"knn","2":"xgboost","3":"-0.03613281","4":"-0.31365102","5":"0.24707251","6":"0.02134899","7":"negligible"},{"1":"classification tree","2":"xgboost","3":"-0.30859375","4":"-0.55336949","5":"-0.01475470","6":"0.01955413","7":"small"},{"1":"random forest","2":"xgboost","3":"0.02050781","4":"-0.26109989","5":"0.29889866","6":"0.02128913","7":"negligible"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
</script>
</div>

```
eff.m <- reshape2::acast(eff.size, treatment ~ control, value.var='delta', margins=F)
corrplot::corrplot(eff.m, diag = F)
```

<img src="{{ site.baseurl }}/assets/plots/sdm-9-4.png" 
    style="width:98%; margin:5px; min-width: 320px; max-width: 850px; height: auto;">

Below are those pairwise comparisons deemed significant by Dunn's test, and their corresponding
effect size results (Cliff's Delta, Lower/Upper 95% C.I., etc.).

```
get.trtmnt <- function(comp) {
  stringr::str_split(comp, " - ")[[1]][[1]]
}

get.cntrl <- function(comp) {
  stringr::str_split(comp, " - ")[[1]][[2]]
}

sig.efsz <- data.table(
  comp=dunn.res$comparisons,
  p.adj=dunn.res$P.adjusted
) %>%
  .[, .(treatment=sapply(comp, get.trtmnt), 
        control=sapply(comp, get.cntrl),
        p.adj)] %>%
  eff.size[., on=.(treatment, control)] %>%
  .[p.adj <= 0.05]

sig.efsz
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script data-pagedtable-source="" type="application/json">
{"columns":[{"label":["treatment"],"name":[1],"type":["chr"],"align":["left"]},{"label":["control"],"name":[2],"type":["chr"],"align":["left"]},{"label":["delta"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["lower"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["upper"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["variance"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["magnitude"],"name":[7],"type":["ord"],"align":["right"]},{"label":["p.adj"],"name":[8],"type":["dbl"],"align":["right"]}],"data":[{"1":"ipp","2":"knn","3":"-0.5195312","4":"-0.7217925","5":"-0.2355127","6":"0.01560406","7":"large","8":"0.0024658617"},{"1":"ipp","2":"logistic regression","3":"-0.4902344","4":"-0.6983050","5":"-0.2057726","6":"0.01606384","7":"large","8":"0.0183257760"},{"1":"ipp","2":"random forest","3":"-0.5566406","4":"-0.7501116","5":"-0.2753969","6":"0.01478261","7":"large","8":"0.0005044008"},{"1":"ipp","2":"xgboost","3":"-0.4843750","4":"-0.6986386","5":"-0.1903817","6":"0.01717800","7":"large","8":"0.0017484719"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
</script>
</div>

## Conclusion

This study primarily aimed to evaluate machine learning models in Species Distribution Modeling
(SDM), comparing their performance to traditional models like MaxEnt and IPP. The findings 
indicate that models such as KNN, Logistic Regression, Random Forest, and XGBoost generally 
outperform IPP models in terms of accuracy. However, it is essential to note that these results, 
drawn from analyzing eight bird species across four U.S. states, are not definitive. The limited 
scope of the study suggests that the conventional approach of relying solely on models like IPP 
or MaxEnt may not always yield the best outcomes in SDM. Future research, possibly extending to a 
broader range of species and geographical locations, is necessary to validate and expand upon 
these findings. This study underscores the potential advantages of integrating machine learning 
techniques in SDM, encouraging a more nuanced approach in model selection for researchers in this 
field.
