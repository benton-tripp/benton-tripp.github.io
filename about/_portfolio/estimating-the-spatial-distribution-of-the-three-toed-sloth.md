---
layout: page
order: 3
title: Estimating the Spatial Distribution of the Three-Toed Sloth
permalink: /portfolio/estimating-the-spatial-distribution-of-the-three-toed-sloth
---

*July, 2023*

Github Page: <a href="https://github.com/benton-tripp/slothful-seer" target="_blank">https://github.com/benton-tripp/slothful-seer</a>


## Problem

This” project addresses the challenge of predicting the spatial distribution of the three-toed sloth using presence-only data. In many ecological applications, only positive observations are available; thus, the project implements presence-only modeling techniques to estimate the probability of occurrence across a study area. The primary objective is to use spatial statistical methods to derive meaningful predictions from sparse and geographically biased data.

The project was developed as part of a Masters of Geospatial Information Science and Technology and leverages multiple modeling approaches. Baseline methods such as the Inhomogeneous Poisson Process (IPP) and Maximum Entropy (MaxEnt) are implemented alongside more advanced models including logistic regression (GLM) with regularization, decision trees, and random forests. This suite of models facilitates both methodological comparisons and an understanding of how different analytical approaches perform under similar data conditions

## Analysis Procedure

A core component of the project is the rigorous data partitioning and stratification process. The dataset is split into training and testing subsets using stratified sampling methods that ensure both spatial representativeness and categorical consistency (for example, based on biome classifications). This approach mitigates the risk of spatial bias and allows for reliable cross-validation. The use of stratification techniques was critical in managing the inherent imbalance in presence-only data.

The modeling workflow integrates multiple analytical procedures. Baseline models are first established using IPP and MaxEnt, which serve as reference points. Advanced models, including a GLM with elastic net regularization, a classification tree, and a random forest, are then trained using a user-defined grid of hyperparameters. Cross-validation is employed to optimize model performance and prevent overfitting. The evaluation process involves comparing metrics such as recall, specificity, and overall accuracy, all of which are visualized through probability rasters, confusion matrices, and other diagnostic plots.

## Results

The implementation yielded a robust modeling framework capable of generating spatial predictions in the form of probability rasters and performance tables. Each model was evaluated based on key performance metrics, which allowed for direct comparison between baseline and advanced methods. The visualization outputs include interactive maps and static plots that clearly delineate predicted versus observed occurrences, offering valuable insights into the spatial dynamics of the species.

The integration of multiple modeling approaches demonstrated that while baseline models provide a solid reference, advanced models can capture non-linear relationships and improve prediction accuracy when appropriately tuned. The project’s evaluation module effectively highlights areas of strong predictive performance and points out potential weaknesses in the model configurations. This systematic evaluation reinforces the importance of model selection and parameter tuning in spatial predictive analytics.

<img src="{{ site.baseurl }}/assets/images/modeling-img1.png" style="width:80%; max-width: 800px; min-width: 300px;">

*Figure 1. Observed occurrence points for the three-toed sloth across its known distribution in Central and South America. Data points represent confirmed presence locations.*

<img src="{{ site.baseurl }}/assets/images/modeling-img2.png" style="width:80%; max-width: 800px; min-width: 300px;">

*Figure 2. Comparison of predicted probability rasters for the distribution of the three-toed sloth derived from five models: Inhomogeneous Poisson Process (IPP), MaxEnt, Generalized Additive Model (GAM), Decision Tree (Tree), and Random Forest. Probability values range from low (white) to high (dark green).*

<img src="{{ site.baseurl }}/assets/images/modeling-img3.png" style="width:80%; max-width: 800px; min-width: 300px;">

*Figure 3. Interactive prediction results at a selected geographic location within the study area. The displayed table provides model-derived probabilities and associated environmental covariates at the selected location (marked in blue).*

## Reflection

The project provided a comprehensive learning experience in the application of spatial statistics to real-world ecological data. The rigorous approach to data partitioning and cross-validation underscored the importance of stratified sampling in addressing spatial biases. Working with a combination of baseline and advanced models highlighted both the strengths and limitations of different predictive techniques when applied to presence-only data. Reflecting on the process, the project reinforced key analytical competencies such as model tuning, validation, and evaluation. The iterative workflow—integrating data management, model training, and detailed evaluation—served as an effective framework for assessing predictive performance. This experience has contributed significantly to a deeper understanding of geospatial predictive modeling and has laid a strong foundation for further research in spatial analytics.