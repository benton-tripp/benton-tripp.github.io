Sharing is Caring - A News Article Popularity Analysis
================

## Using R to Analyze News Article Data and Predict Shares

I recently completed an analysis of an online news popularity data set
from UCI, using R (see <https://doi.org/10.24432/C5NS3V>). The objective
of this analysis was to gain an understanding of the factors influencing
online news popularity, and develop a predictive model capable of
accurately estimating the total shares of a given article given these
factors.

The analysis can be accessed via the following locations:

- [Online News Popularity Analysis - Main Page](https://benton-tripp.github.io/news-popularity-analysis/)
  (published via Github Pages)
  - Published from the project README.md file.
  - Links to each of the six published reports of the analysis (one for
    each channel type) can be found at this location, along with a brief
    description of the code used to automate the rendering process.
- [Github
  Repository](https://github.com/benton-tripp/news-popularity-analysis)
  - This is the repository containing all of the code/files used in the
    project (and through which the Github Pages posts are published)

## Analysis Reflections

### What Would I do Differently?

Initially, I didn’t really do much to handle multi-collinearity or
linear dependence among the predictors. Although I eventually went back
and did more to account for these issues (described in further detail in
the published reports), I would probably do this much sooner during
EDA/model training in the future.

### What was the Most Difficult Part?

One of the largest challenges was the nature of the data itself. There
are a lot of outliers in the data, which can make building an accurate
predictive model difficult. For example, consider the distribution of
the business articles:

| Measure        | Value  |
|:---------------|:-------|
| min            | 1      |
| lower.quartile | 952    |
| median         | 1400   |
| mean           | 3063   |
| upper.quartile | 2500   |
| max            | 690400 |

The maximum number of shares was so far above the normal range of the
data, I was unable to accurately predict that value (which was in my
test data) despite all of my efforts.

### What are my Key Takeaways?

I spent a bit of time exploring variations of linear regression that are
more robust to outliers. I first tested using different regularization
methods (e.g., LASSO, Ridge Regression, Elastic Net), but didn’t have
much luck. In the end, I settled on using quantile regression. At a high
level, quantile regression estimates the conditional quantile of the
response variable, as opposed to traditional OLS which models the mean
of the response variable (given the predictors). Consider how a median
of a distribution tends to be more robust to outliers than a mean of the
same distribution. Quantile regression uses this to its advantage to
help reduce some of the bias introduced by outliers in the training
data.
