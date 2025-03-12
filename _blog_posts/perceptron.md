---
layout: post
title: The Perceptron Algorithm (in R and Python)
permalink: /blog_posts/perceptron
---

## Introduction
The Perceptron algorithm is possibly the simplest of binary classification
algorithms. Originally invented in 1958 by Frank Rosenblatt, it was quickly 
discovered that the constraints of a single-layer perceptron make the 
classification method impractical. However, multi-layer perceptrons are still
used today as the basic building-block of "feed-forward" neural networks. As a 
result, a keen understanding of this algorithm is critical to any aspiring data
scientist.

This article is meant to serve as a high-level overview of the mathematics
behind the algorithm, as well as an example of how one might go
about programming their own Perceptron algorithm. Although, the code examples
are written in the R language, I have similar functions in Python that you
can find [here](https://github.com/benton-tripp/perceptron-decepticon/blob/3ada285f56ed09a26e77ef74e2800bac16fd94f1/notebooks/perceptron.ipynb).

### R Packages
Before we begin, we can go ahead and import the required R packages. If they are
not currently installed, that can be done first. The packages used in this 
tutorial are: <br>
- [data.table](https://cran.r-project.org/web/packages/data.table/index.html) <br>
- [ggplot](https://ggplot2.tidyverse.org/) <br>
- [plotly](https://plotly.com/r/) <br>
- [matlib](https://cran.r-project.org/web/packages/matlib/index.html) <br>

```
# Install packages if needed

if (!require(data.table)) install.packages('data.table')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(plotly)) install.packages('plotly')
if (!require(matlib)) install.packages('matlib')

library(data.table)
library(ggplot2)
library(plotly)
library(matlib)
```

## Math Fundamentals of the Perceptron Algorithm
In order to understand the Perceptron algorithm, a basic knowledge of the
following mathematical principles must first be understood:
- Vectors (specifically, how to calculate the direction and norm)
- The Dot Product

I am going to assume that majority of people reading this article have prior 
knowledge of each of these fundamentals. If you do not, a great resource can be
found [here](https://www.khanacademy.org/math/precalculus/x9e81a4f98389efdf:vectors). I am not going to attempt to
teach the fundamentals, because there are already so many (better) resources 
available.

## Linearly Separable Data
As mentioned preciously, one of the challenges with the single-layer perceptrons
is that the constraints limit the algorithm's ability to accurately classify
data. Specifically, the primary constraint is that the data *must* be linearly
separable. 

The definition of linearly separable data is pretty self-explanatory.
At a high level, for data to be linearly separable it simply means that there 
exists some plane that can literally "separate" the data. For example, consider 
the following dataset of two classes:

```
DT <- data.table(
  class_ = c(rep(1, 31), rep(-1, 31)),
  x = c(seq(1, 4, .1), seq(7, 10, .1)),
  y = c(runif(31, 1, 4), runif(31, 7, 10))
)

ggplot(DT, aes(x=x, y=y)) + 
  geom_point(colour=c(rep('blue',31), rep('darkred', 31)), size=3) + 
  theme(legend.position='none',
        axis.text.x=element_blank(),
        axis.text.y=element_blank())
```

<img src="{{ site.baseurl }}/assets/plots/perceptron-plt-1.png" style="width: 90%; min-width: 320px; max-width: 600px; height: auto; border: 1px solid #aaa;">

It is pretty obvious that you can draw some arbitrary line that separates the 
two different colored data points:

```
ggplot(DT, aes(x=x, y=y)) + 
  geom_point(colour=c(rep('blue',31), rep('darkred', 31)), size=3) + 
  geom_abline(slope=-1.25, intercept=12) + 
  theme(legend.position='none',
        axis.text.x=element_blank(),
        axis.text.y=element_blank())
```

<img src="{{ site.baseurl }}/assets/plots/perceptron-plt-2.png" style="width: 90%; min-width: 320px; max-width: 600px; height: auto; border: 1px solid #aaa;">

Similarly, consider the following three-dimensional data:

```
DT <- DT[, z := c(runif(31, 1, 4), runif(31, 7, 10))]

scene = list(
  camera = list(
    center = list(x = 0, y = 0, z = 0), 
    eye = list(x = 1.75, y = .75, z = .25)
    ),
  xaxis=list(title='x', showticklabels=F),
  yaxis=list(title='y', showticklabels=F),
  zaxis=list(title='z', showticklabels=F)
)

plot_ly(data=DT,
        x=~x, y=~y, z=~z, 
        color=~as.factor(class_),
        type='scatter3d',
        mode='markers', 
        colors=c('darkred', 'blue'),
        showlegend=F) %>%
  layout(scene=scene) %>%
  config(
    displaylogo = FALSE, 
    modeBarButtonsToRemove = c("hoverclosest", 
                               "resetCameraLastSave", 
                               "orbitRotation", 
                               "turntableRotation", 
                               "zoomIn3d", 
                               "zoomOut3d", 
                               "toImage")
  )
```

<iframe src="{{ site.baseurl }}/assets/plots/perceptron-plt-1.html" style="width: 90%; min-width: 320px; max-width: 600px; height: 500px; border: none;"></iframe>

In this instance, a line would be insufficient when attempting to separate the 
two data classes. Instead, a two-dimensional plane can be used:

```
plot_ly(data=DT,
        colors=c('darkred', 'gray', 'blue')) %>%
  add_trace(x=~seq(min(x), max(x), .1)*1.5, 
            y=~seq(min(y), max(y), .1)*1.5, 
            z=~matrix(6, nrow=62, ncol=62),
            type='surface') %>%
  add_trace(x=~x, y=~y, z=~z, 
            color=~as.factor(class_),
            type='scatter3d',
            mode='markers',
            showlegend=F) %>%
  layout(scene=scene) %>% 
  hide_colorbar() %>%
  config(
    displaylogo = FALSE, 
    modeBarButtonsToRemove = c("hoverclosest", 
                               "resetCameraLastSave", 
                               "orbitRotation", 
                               "turntableRotation", 
                               "zoomIn3d", 
                               "zoomOut3d", 
                               "toImage")
  )
```

<iframe src="{{ site.baseurl }}/assets/plots/perceptron-plt-2.html" style="width: 90%; min-width: 320px; max-width: 600px; height: 500px; border: none;"></iframe>

Of course, this gets a little harder to visualize - or even explain - 
as dimensionality continues to increase. As a result, a more generalized term
to describe this linear separation is used: a *hyperplane*.

### Hyperplanes
According to the [Wikipedia definition](https://en.wikipedia.org/wiki/Hyperplane): 
"In geometry, a hyperplane is a subspace whose dimension is one less than 
that of its ambient space." This makes sense. Looking at our first example, the
dimensionality of the data is two, and the line drawn to separate the data is
one-dimensional. In the second example, the data is three-dimensional and the
separation is two.

## Perceptron Algorithm
Now assuming you have some basic knowledge of vectors and we have reviewed the 
the meaning of linear separability and hyperplanes, we can build on this 
knowledge to understand how the Perceptron algorithm works.

To begin, consider the following:

\[ 
  f(x) = \begin{cases} 
      1 & \text{if } \vec{x} \cdot \vec{w} + b > 0 \\
      -1 & \text{otherwise} 
   \end{cases} \\
  \text{where}\ \vec{w} \ \text{ and } \vec{x} \text{ are real-valued vectors (} \vec{w}\text{ represents weights); } \\
  \vec{x} \cdot \vec{w} = \sum_{i=1}^k \vec{x}\vec{w} \text{, where } k \text{ is the number of inputs; } b \text{ is the bias}
\]

In more simple terms, we are saying that there are two classes in the data (**x**),
hence a binary classification problem. We are calling the classes 1 and -1, but 
they could be anything. Supposing that the two classes are linearly separable, 
then there must exist some **w** such that when the the dot product (plus bias) between 
**x** and **w** is positive, the class is 1. And when the dot product is negative,
the class is negative 1. 

To find the optimal **w**, we can first initialize the vector as random values.
Technically, the values could all be zero. I chose to initialize random values
between zero and one.

In this example, we will work on the two dimensional data given in the 
example earlier in the article demonstrating linear separability. Because the 
data has two dimensions (*x*, *y*), the vector **w** will also have two dimensions,
plus some bias *b*. 

To make things a little easier, we can actually completely remove the bias from
the equation. To do this, first recall that:

$$\vec{x} \cdot \vec{w} = \sum_{i=1}^k \vec{x}\vec{w}$$

This means that:

$$\vec{x} = (x_1, x_2, ... x_k) \text{ and } \vec{w} = (w_1, w_2, ...w_k)$$

To remove the bias *b* we can simple add a zero component to **x** and **w**, 
such that:

$$\vec{x} = (x_0, x_1, ... x_k) \text{ where } x_0 = 1, \text{ and } \vec{w} = (w_0, w_1, ...w_k)$$

Now, the original function can be written as:

$$
  f(x) = \begin{cases} 
      1 & \text{if } \vec{x} \cdot \vec{w} > 0 \\
      -1 & \text{otherwise} 
   \end{cases} \\
   \text{where } \vec{x} = (x_0, x_1, ... x_k) \text{ and } x_0 = 1; \text{ and } \vec{w} = (w_0, w_1, ...w_k)
$$

In our two dimensional example, this means that **w** should be initialized with
3 values rather than the original two. Our matrix **x** should also be augmented
with an additional "column" of data, where each value in the additional column
is equal to 1. In R, these steps look like the following code:

```
# Generate random weights (initialize `w` with 3 values)
w <- runif(3, min=0, max=1)

# Remove `z` since we are using 2D data in this example
# Augment the data such that `x0` = 1
DT[, z := NULL][, x0 := 1]

```

Now that **w** and **x** are both initialized and augmented to account for any
bias *b*, we can make some initial predictions using the formula *f(x)*. This
can be done by taking the sign of the dot product of **x** and **w**. We will
also identify where the initial predictions are incorrect, since those will
play a key role in the next step of the algorithm.

```
# Predict class (either -1 or 1) at each point by taking the dot product
pred <- apply(DT[, 2:4], 1, FUN = function(x) sign(sum(x * w)))

# Get incorrect predictions
incorrect_preds <- DT[class_ != pred]
```

In the next step, we are going to take a random sample from the mis-classified
predictions, multiply each of the three **x** values of the sample by the actual
class (either 1 or -1), and then add each of the three respective values of the
sample to the vector **w**. This gives us an updated vector **w**. Then, the 
original process is repeated. Predictions are made using the new **w**, incorrect
predictions are identified and used to update **w**. The process iterates until
there are no incorrect predictions (or a pre-defined threshold is met).

To understand the reasoning behind this process, it helps to picture what is 
actually happening to the vectors. As an example, consider the following
vectors, where `x_sample1` is a random sample taken from the incorrect 
predictions of some arbitrary vector **x**, and `w1` is the weight vector
used in the calculation:

```
x_sample1 <- c(.3, .8)
w1 <- c(.6, .1)
```

Visualized, we get the following (where *theta* is the angle between `x_sample1`
and `w1`):

```
plot(1, type='n', xlim=c(-.1, 1), ylim=c(-.1, 1), 
     xlab='', ylab='', xaxt='n', yaxt='n')
grid()

vectors(x_sample1, labels=expression(x))
vectors(w1, labels=expression(w))

arc(w1, c(0, 0), x_sample1, d=.1, absolute=T)
text(.15, .15, expression(theta), cex=1.5)
```

<img src="{{ site.baseurl }}/assets/plots/perceptron-plt-3.png" style="width: 90%; min-width: 320px; max-width: 600px; height: auto; border: 1px solid #aaa;">

We can see that the angle between the two vectors `x_sample1` and `w1` is 
less than 90°. We also know that the predicted class of `x_sample1` is incorrect,
since the sample was taken from the mis-classified data. This means that the
actual class of `x_sample1` is the opposite of the dot product between 
`x_sample1` and `w1`. Now, we need to adjust `w1` in some way so 
that the sign of the dot product between `x_sample1` and `w1` is equal to
the class of `x_sample1`.

To do this, let's first see what the predicted value of `x_sample1` is:

```
print(paste0('The predicted class is: ', sign(x_sample1 %*% w1)[1]))
```

Because the predicted class is 1 and we know the prediction in this example is
incorrect, we can conclude that the actual class is -1. To get a predicted value
of 1, we need the dot product between `x_sample1` and our **w** value to be 
negative.

To do this, we can simply subtract `x_sample1` from `w1`:

```
plot(1, type='n', xlim=c(-.05, .7), ylim=c(-.8, .8), 
     xlab='', ylab='', xaxt='n', yaxt='n')
grid()

new_w1 <-  w1 - x_sample1

vectors(x_sample1, labels=expression(x))
vectors(w1, labels=expression(w))
vectors(new_w1, labels=expression(w - x), col='darkred')

arc(w1, c(0, 0), x_sample1, d=.3, absolute=F)
text(.2, .3, expression(theta), cex=1.5)

arc(new_w1, c(0, 0), x_sample1, d=.2, col='darkred')
text(.12, .16, expression(beta), cex=1.5, col='darkred')
```

<img src="{{ site.baseurl }}/assets/plots/perceptron-plt-4.png" style="width: 90%; min-width: 320px; max-width: 600px; height: auto; border: 1px solid #aaa;">

Which produces the following prediction:

```
print(paste0('The predicted class is: ', sign(x_sample1 %*% new_w1)[1]))
```

Now, let's take two more arbitrary vectors `x_sample2` and `w2`. Let the 
angle between `x_sample2` and `w2` is greater than 90°, the actual 
class of `x_sample2` be 1, and the predicted value be -1:


```
plot(1, type='n', xlim=c(-.8, 1), ylim=c(-.1, .6), 
     xlab='', ylab='', xaxt='n', yaxt='n')
grid()

x_sample2 <- c(-.7, .45)
w2 <- c(.9, .05)

vectors(x_sample2, labels=expression(x))
vectors(w2, labels=expression(w))

arc(w1, c(0, 0), x_sample2, d=.1)
text(.05, .07, expression(theta), cex=1.5)
```

<img src="{{ site.baseurl }}/assets/plots/perceptron-plt-5.png" style="width: 90%; min-width: 320px; max-width: 600px; height: auto; border: 1px solid #aaa;">

```
print(paste0('The predicted class is: ', sign(x_sample2 %*% w2)[1]))
```

Because the predicted class is -1 and we need it to be 1, we actually need to
decrease the angle between `x_sample` and `w2`. To do this, we add them together
instead of subtract:

```
plot(1, type='n', xlim=c(-.9, 1.1), ylim=c(-.05, .6), 
     xlab='', ylab='', xaxt='n', yaxt='n')
grid()

new_w2 <-  w2 + x_sample2

vectors(x_sample2, labels=expression(x))
vectors(w2, labels=expression(w))
vectors(new_w2, labels=expression(w + x), col='darkblue')

arc(w2, c(0, 0), x_sample2, d=.1)
text(-.02, .08, expression(theta), cex=1.2)

arc(new_w2, c(0, 0), x_sample2, d=.15, col='darkblue')
text(-.05, .18, expression(beta), cex=1.2, col='darkblue')

```

<img src="{{ site.baseurl }}/assets/plots/perceptron-plt-6.png" style="width: 90%; min-width: 320px; max-width: 600px;height: auto; border: 1px solid #aaa;">

```

print(paste0('The predicted class is: ', sign(x_sample2 %*% new_w2)[1]))

```

As mentioned before, this process iterates with each new **w** until there are
either no incorrect predictions, or some pre-defined threshold is met. If the 
data is not linearly separable, then **w** will never converge.

```

# Iteratively update weights
repeat{
  w <- w + unname(unlist(
    incorrect_preds[sample(1:nrow(incorrect_preds), 1)][, .SD * class_, .SDcols=2:4][1,]
    ))
  pred <- apply(DT[, 2:4], 1, FUN = function(x) sign(sum(x * w)))
  incorrect_preds <- DT[class_ != pred]
    if (nrow(incorrect_preds) == 0) break
}

# Plot
ggplot(data=DT, aes(x=x, y=y)) +
  geom_point(color=c(rep('blue',31), rep('darkred', 31)), size=3) +
  geom_abline(slope=-w[[2]] / w[[1]], intercept=-w[[3]] / w[[1]]) + 
  theme(legend.position='none',
        axis.text.x=element_blank(),
        axis.text.y=element_blank())

```

<img src="{{ site.baseurl }}/assets/plots/perceptron-plt-7.png" style="width: 90%; min-width: 320px; max-width: 600px; height: auto; border: 1px solid #aaa;">

We've done it! Feel free to try all of the code yourself. Below is a full R 
script with all of the steps in the algorithm defined as functions, and if
haven't already feel free to check out the same code written in Python [here](https://github.com/benton-tripp/perceptron-decepticon/blob/3ada285f56ed09a26e77ef74e2800bac16fd94f1/notebooks/perceptron.ipynb).

```
#' Predict class of data by taking the sign of the dot product of x & w
#' 
#' @param data linearly separable data.table (class_, x0, x, y)
#' @param w vector of real-valued weights
#' 
#' @return vector of predicted classes
#' 
predict.perc <- function(data, w) {
  pred <- apply(data[, 2:4], 1, FUN = function(x) sign(sum(x * w)))
  return(pred)
}


#' Iterate over data to generate optimal w vector, slope, intercept, and predicted values
#' 
#' @param data linearly separable data.table (class_, x0, x, y)
#' @param incorrect_preds subset of data containing initial incorrect predictions
#' @param w starting w vector
#'
#' @return list containing w, slope/intercept, and predicted values
#' 
update_preds <- function(data, incorrect_preds, w) {
  repeat{
    w <- w + incorrect_preds[sample(1:nrow(incorrect_preds), 1)][, .SD * cls, .SDcols=2:4][1,]
    pred <- predict.perc(data, w)
    incorrect_preds <- data[cls != pred]
    if (nrow(incorrect_preds) == 0) break
  } 
  return(
    list(
      w = c(unname(unlist(w))),
      slope = -w[[2]] / w[[1]],
      intercept = -w[[3]] / w[[1]],
      pred = pred
    )
  )
}


#' Generate w vector, slope, intercept, and predicted values
#' 
#' @param data linearly separable data.table (class_, x, y, x0)
#' @param verbose whether to print output line properties
#'
#' @return list containing w, slope/intercept, and predicted values
#' 
perceptron <- function(data, verbose=F) {
  w <- runif(3, min=0, max=1)
  preds <- predict.perc(data, w)
  incorrect_preds <- data[class_ != pred]
  p <- update_preds(data, incorrect_preds, w)
  if (verbose == T) {
    message(paste0('The calculation has a slope of ', round(p$slope, 3), 
                   ' and intercept at ', round(p$intercept, 3)))
  }
  return(p)
}

```

## References
https://en.wikipedia.org/wiki/Perceptron <br>
https://bitbucket.org/syncfusiontech/svm-succinctly/src/master/


