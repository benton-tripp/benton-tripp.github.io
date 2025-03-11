---
layout: page
order: 2
title: Species Distribution Model Benchmark Study
permalink: /portfolio/species-distribution-model-benchmark-study
---

*December, 2023*

Github Page: <a href="https://github.com/benton-tripp/presence-only-sdm" target="_blank">https://github.com/benton-tripp/presence-only-sdm</a>

## Problem

The objective of this project was to benchmark various species distribution modeling (SDM) techniques using bird observation data from eBird alongside comprehensive environmental predictors. The goal was to accurately predict the distribution of eight bird species across four U.S. states (Colorado, North Carolina, Oregon, and Vermont). Existing methods have limitations due to biased observational data and variability in species-environment relationships, necessitating robust modeling frameworks that effectively handle these complexities.

The project specifically aimed to compare traditional SDM approaches, such as logistic regression and Bioclim models, with more advanced methods including Inhomogeneous Poisson Process (IPP) models, Maximum Entropy (MaxEnt), and machine learning algorithms (Random Forests, Gradient Boosting Machines). By systematically evaluating these approaches, the project sought to identify best practices for species distribution modeling that could guide conservation planning and ecological research.

## Analysis Procedures

Data Acquisition and Preprocessing

- Collected spatial boundaries and bird observation data (eBird) from 2016-2019.
- Retrieved diverse environmental rasters, including DEM, NDVI, climate variables, land cover, and hydrography.
- Applied R packages (sf, terra, stars) for spatial data manipulation and transformation (EPSG:5070).
- Performed rigorous data cleaning and spatial joining to integrate environmental variables at bird observation locations.

Exploratory Data Analysis (EDA)

- Conducted visualization of spatial and temporal observation patterns with ggplot2.
- Summarized data to identify potential observation biases and seasonal trends.
- Evaluated relationships between species occurrences and environmental variables to inform model selection.

Modeling Methods

- Implemented baseline logistic regression and Bioclim models using climatic predictors.
- Developed IPP models leveraging maximum likelihood estimation to handle presence-only data effectively.
- Explored pseudo-absence resampling methods to assess impacts on model accuracy.
- Utilized MaxEnt models via the dismo package, adjusting parameters to prevent overfitting.
- Applied machine learning models, specifically Random Forests and Gradient Boosting Machines (caret, xgboost), performing hyperparameter tuning and validation through cross-validation.

## Results

The benchmarking study provided clear comparisons across diverse SDM techniques. Initial baseline models demonstrated moderate predictive ability, with logistic regression and Bioclim approaches setting a performance standard. Subsequent implementation of Inhomogeneous Poisson Process (IPP) and MaxEnt models improved predictive accuracy substantially, particularly through strategic sampling of pseudo-absence data points. IPP models highlighted enhanced performance through refined pseudo-absence sampling methods, with noticeable gains in predictive accuracy.

Machine learning approaches, particularly Gradient Boosting Machines, delivered superior performance with the highest overall accuracy and predictive robustness (Figure 1). Fine-tuning of these models significantly improved results, demonstrating their strength in handling complex ecological relationships and observational biases. Metrics such as Accuracy, Sensitivity, Specificity, and F1 scores clearly illustrated the relative strengths and weaknesses of each method, underscoring Gradient Boosting Machines and Random Forests as consistently high-performing methods. The methodological insights emphasized the importance of strategic pseudo-absence sampling, feature selection, and rigorous model validation procedures.

<img src="{{ site.baseurl }}/assets/images/programming-img-1.png" style="width:80%; height:auto; max-width: 700px; min-width: 300px;">

*Figure 1: Comparative performance metrics for multiple species distribution modeling techniques (IPP, MaxEnt, Logistic Regression, k-Nearest Neighbors, Classification Tree, Random Forest, XGBoost). Metrics presented include Accuracy, Sensitivity, Specificity, and F1 score. Red points indicate mean values.*

## Reflection

This project reinforced critical lessons regarding data quality, methodological rigor, and model interpretability in species distribution modeling. Data preprocessing and exploratory analysis emerged as fundamental steps for ensuring meaningful results, underscoring the importance of carefully prepared inputs in spatial ecological modeling. Challenges encountered included effectively addressing biases inherent in citizen science data and selecting appropriate pseudo-absence sampling techniques. Future studies might explore incorporating additional dynamic covariates or adopting ensemble modeling approaches to further enhance predictive performance and ecological insight.