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

### Confirming ANOVA Assumptions

The Shapiro-Wilk test can be used to confirm whether the residuals of an ANOVA model are
normally distributed. The null hypothesis is that the data is normally distributed.
If $p < alpha = 0.05$, then the null hypothesis is rejected, and it can be concluded that the
residuals of the ANOVA are not normally distributed

```
shapiro.test(residuals(anova.res))
```

### Non-Parametric Alternative to ANOVA

The Kruskal-Wallis Rank Sum Test is a non-parametric alternative to one-way ANOVA. 
It's used to compare more than two groups and does not assume a normal distribution.
The null hypothesis for this test is that all groups have the same median, with a
threshold of 0.05.

```
kruskal.test(Accuracy ~ model.type, data = model.scores)
```

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

```
eff.m <- reshape2::acast(eff.size, treatment ~ control, value.var='delta', margins=F)
corrplot::corrplot(eff.m, diag = F)
```

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
