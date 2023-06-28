R you ready for this? My favorite features of R
================

## What’s the best part about the R language? Well… It’s hard to choose

I have never been good at “picking favorites.” I absolutely despise
get-to-know-you questions like “What is your favorite food?” Or “Where
is your dream vacation?” I am the type of person who appreciates
variety, practicality, and novelty. These are the very qualities that
make the R language so appealing to me. For this post, (in the interest
of brevity), I will focus in on novelty and highlight a selection of
unique features exclusive to the R language.

### Formula Syntax

While other languages might have similar ways to specify models, the
specific syntax and flexibility of R’s formula interface is quite
unique. It’s one of the features that makes R particularly well-suited
for statistical modeling and data analysis.

For example, a simple linear regression:

``` r
# Load the built-in 'mtcars' dataset
data(mtcars)

# Fit a linear regression model using the formula syntax
model <- lm(mpg ~ cyl + hp, data = mtcars)

# Print the model summary
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = mpg ~ cyl + hp, data = mtcars)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.4948 -2.4901 -0.1828  1.9777  7.2934 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 36.90833    2.19080  16.847  < 2e-16 ***
    ## cyl         -2.26469    0.57589  -3.933  0.00048 ***
    ## hp          -0.01912    0.01500  -1.275  0.21253    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.173 on 29 degrees of freedom
    ## Multiple R-squared:  0.7407, Adjusted R-squared:  0.7228 
    ## F-statistic: 41.42 on 2 and 29 DF,  p-value: 3.162e-09

In this example, `mpg ~ cyl + hp` is a formula that specifies a model
where the miles per gallon (mpg) is explained by the number of cylinders
(cyl) and horsepower (hp).

### Piping

Although the concept itself is not unique (other languages such as
Python or JavaScript have “chaining”), piping is an elegant way of
keeping R code organized, interpretable, and efficient. The pipe
operator `%>%` (which more recently has been updated to alternatively be
`|>`) is simply placed between functions, linking them into a continuous
process. For example:

``` r
# Load the necessary library
library(magrittr) 
# Can also be loaded through dplyr, which requires magrittr

# Use piping to filter and summarize the mtcars dataset
mtcars %>%
  subset(cyl == 6) %>%
  summarise(mean_mpg = mean(mpg), mean_hp = mean(hp))
```

    ##   mean_mpg  mean_hp
    ## 1 19.74286 122.2857

### Built-In Data Frame Data Structure

In R, data frames are built into the base language, and many of R’s
built-in functions and packages are designed to work seamlessly with
them. This makes data frames in R particularly easy (and efficient) to
work with. Although the concept itself is not unique to R (e.g., the
Python language is known for its pandas library), the implementation of
a data frame type within the base language sets R apart as a powerful
language for data analytics and data science.

The built-in `mtcars` R dataset, as used in the prior two examples, is
an example of an R data frame (see the data structure below):

``` r
str(mtcars)
```

    ## 'data.frame':    32 obs. of  11 variables:
    ##  $ mpg : num  21 21 22.8 21.4 18.7 18.1 14.3 24.4 22.8 19.2 ...
    ##  $ cyl : num  6 6 4 6 8 6 8 4 4 6 ...
    ##  $ disp: num  160 160 108 258 360 ...
    ##  $ hp  : num  110 110 93 110 175 105 245 62 95 123 ...
    ##  $ drat: num  3.9 3.9 3.85 3.08 3.15 2.76 3.21 3.69 3.92 3.92 ...
    ##  $ wt  : num  2.62 2.88 2.32 3.21 3.44 ...
    ##  $ qsec: num  16.5 17 18.6 19.4 17 ...
    ##  $ vs  : num  0 0 1 1 0 1 0 1 1 1 ...
    ##  $ am  : num  1 1 1 0 0 0 0 0 0 0 ...
    ##  $ gear: num  4 4 4 3 3 3 3 4 4 4 ...
    ##  $ carb: num  4 4 1 1 2 1 4 2 2 4 ...
