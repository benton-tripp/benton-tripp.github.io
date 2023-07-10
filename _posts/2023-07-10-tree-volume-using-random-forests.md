Tree Volume using Random Forests
================

## Introduction

The Random Forest algorithm is a versatile machine learning algorithm
capable of performing both regression and classification tasks. It is a
type of ensemble learning method, where a group of weak models (in this
case decision trees) combine to form a strong model. In Random Forest,
we grow multiple decision trees as opposed to a single tree. To classify
a new object based on attributes, each tree gives a classification. The
forest chooses the classification having the most votes (for all the
trees in the forest) and in case of regression, it takes the average of
outputs by different trees.

I often prefer Random Forests over other methods like linear regression,
gradient boosting, or k-nearest neighbors for several reasons:

- Robustness: Random Forests are less likely to overfit than other
  methods due to the randomness injected into the model building
  process.
- Versatility: They can handle both numerical and categorical data, and
  they can handle missing data. It is also unnecessary to scale data
  prior to model training.
- Parallelizable: The process of building trees is easily parallelizable
  as each tree is built independently of the others.
- Performance: They often provide very good predictive performance with
  little tuning required.
- Feature Importance: They provide a built-in method for estimating
  feature importance.

## Random Forest Algorithm

The Random Forest algorithm works by creating multiple decision trees
and merging them together to get a more accurate and stable prediction.
The main concept behind Random Forest is the idea of bagging, which is
used to reduce the variation in the predictions by combining the result
of multiple decision trees on different samples of the data set.

Here is a simplified version of the Random Forest algorithm:

1.  Select random samples from a given dataset.
2.  Construct a decision tree for each sample and get a prediction
    result from each decision tree.
3.  Perform a vote for each predicted result.
4.  Select the prediction result with the most votes as the final
    prediction.

## Tree Dataset

The `trees` dataset used in this example come from the built in R
`datasets` package:

> *This data set provides measurements of the diameter, height and
> volume of timber in 31 felled black cherry trees. Note that the
> diameter (in inches) is erroneously labelled Girth in the data. It is
> measured at 4 ft 6 in above the ground* (R Documentation \| trees
> {datasets}).

See also Ryan, T. A., Joiner, B. L. and Ryan, B. F. (1976) *The Minitab
Student Handbook*. Duxbury Press.

## Random Forest Example

In the following example, we will use the `trees` dataset to predict the
volume of a tree based on its girth and height using a Random Forest
model. We will use the `caret` package in R to perform the analysis.

The below code block first splits the data into training and testing
sets. Then, it defines a control function and a tuning grid for the
Random Forest model. The model is then fitted on the training data and
used to make predictions on the test data. The performance of the model
on the test data is then evaluated.

``` r
# Import libraries
library(caret)
library(tidyverse)

# Set the seed
set.seed(123)

# Create the index for the train split using the `caret` package
train.idx <- createDataPartition(seq(1, nrow(trees)), p=0.7, list=F, times=1)
  
# Create the training set
df.train <- trees[train.idx,]
  
# Create the testing set
df.test <- trees[-train.idx,]

# Define the control for the train function
ctrl <- trainControl(method = "cv", number = 10)

# Define the tuning grid
rf.tune.grid <- expand.grid(
  mtry=c(1, 2),
  splitrule=c("variance", "extratrees", "maxstat"),
  min.node.size=c(1, 3, 5)
)

# Fit the random forest model
fit <- train(Volume ~ ., 
                data = df.train, 
                method = "ranger", 
                trControl = ctrl, 
                tuneGrid = rf.tune.grid,
                importance = "impurity",
                metric = "RMSE")

# Generate predictions on the test data
yhat <- predict(fit, newdata = df.test)
metrics <- postResample(yhat, obs=df.test$Volume)

# View the metrics of the model on the test data
metrics
```

    ##      RMSE  Rsquared       MAE 
    ## 2.8753340 0.9723468 2.6712143

The model performed well with a high R-squared value of 0.97, indicating
that ~97% of the variance in the tree volume can be explained by the
girth and height of the tree. The Root Mean Square Error (RMSE) was
~2.9, which gives us an idea of how much our predictions deviate, on
average, from the actual values in the dataset. The lower the RMSE, the
better our model’s predictions are.

Next, a scatter plot of the predicted vs. actual volume of the trees in
the test set is generated. The red line represents the line of perfect
prediction (i.e., where the predicted volume equals the actual volume).

``` r
# Plot predicted vs. actual volume (using the test set)
ggplot(data.frame(actual=df.test$Volume, predicted=yhat), 
       aes(x = actual, y = predicted)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(x = "Actual", y = "Predicted", 
       title = paste0("Predicted vs Actual Volume")) +
  theme_minimal()
```

![](C:/Users/bento/st558/benton-tripp.github.io/_posts/2023-07-10-tree-volume-using-random-forests_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

Looking at the actual vs predicted values, we can see that the model’s
predictions are quite close to the actual values. This indicates that
our model is doing a good job of predicting tree volume based on girth
and height.

Finally, the importance of each feature in the model is calculated and
used to generate a bar plot. The importance of a feature is measured by
the increase in the model’s prediction error after permuting the
feature. A higher increase in error indicates a more important feature.

``` r
# Get feature importance of model
feature.imp <- varImp(fit, scale = F)$importance %>% 
  arrange(desc(Overall))

ggplot(data.frame(Variable = rownames(feature.imp), Importance = feature.imp$Overall),
       aes(x=reorder(Variable, Importance), y=Importance)) +
  geom_col(show.legend = F, fill="darkblue") +
  labs(x = "Variable", y = "Importance", 
       title = "Variable Importance") +
  theme_minimal() +
  coord_flip() + 
  theme(text = element_text(size = 11), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 13))
```

![](C:/Users/bento/st558/benton-tripp.github.io/_posts/2023-07-10-tree-volume-using-random-forests_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

The feature importance plot shows that the girth of the tree is the most
important feature for predicting the volume, with an importance score of
~5171. This is significantly higher than the importance score of the
height, which is ~922. This suggests that the girth of a tree is a much
stronger predictor of its volume than the height.
