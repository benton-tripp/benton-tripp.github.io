---
layout: post
title: "SDM Benchmark Study Part 4: Baseline Species Distribution Models"
permalink: /blog_posts/sdm-benchmark-study-part-4-baseline-species-distribution-models
---

## Introduction

For this study, two modeling methods will be used as “baseline” models (i.e., they will be fit to the same data, and used as comparisons for other modeling methods later in the study). The two baseline models chosen - Inhomogeneous Poisson Process models and Maximum Entropy models - were selected because they are generally considered the standard for species distribution modeling problems.

MaxEnt models the probability per grid cell and analyzes data after aggregating them into presence/absence grid cells. In contrast, a Poisson Point Process models the limiting expected count or intensity per unit area, rather than per grid cell. This per area basis, compared to per grid cell, is a significant distinction between the two approaches (Renner & Warton, 2013). Interestingly, Renner and Warton (2013) highlight the mathematical equivalence of the MaxEnt procedure and Poisson regression, asserting that both approaches fit the same model and estimate parameters to maximize the same function up to a constant. Moreover, the MaxEnt and Point Process solutions for grid cell data are proportional, with identical estimates of slope parameters (see my [literature review on species distribution modeling]({{ site.baseurl }}/docs/lit_review.pdf) for a more in-depth discussion).

Of course, despite their similarities and mathematical equivalencies, a Poisson Point Process and MaxEnt model will typically yield different results. As such, each will be considered separately within this study.


## Inhomogeneous Poisson Process (IPP)

An Inhomogeneous Poisson Process (IPP) is a statistical model used for events that occur randomly over space or time. Unlike a regular Poisson process, the rate of event occurrence in an IPP can vary.

IPP models often used in presence-only problems to model the intensity of events across different locations or times. An IPP model is defined as $\lambda(x) = N * f(x)$, where $\lambda(x)$ is the intensity function, $N$ is the total number of events, and $f(x)$ is the density function.

By fitting an IPP to presence-only data, we estimate the underlying intensity function, which tells us how the rate of event occurrence changes across different locations or times. For example, we might hypothesize that the presence of the Cedar Waxwing is influenced by factors like canopy cover, land cover type, and temperature. We use the presence-only data to fit the IPP model, estimating the parameters of $\lambda(x)$ that best fit the observed data.

The intensity function $\lambda(x)$ in an IPP is typically defined as a function of some parameters $\theta$ and the environmental variables at location $x$. For example, we might have: $\lambda(x; \theta) = \exp(\theta_1 * canopy + \theta_2 * land\_cover + ...)$, where $canopy$, $land\_cover$, etc. are the environmental variables, and $\theta_1$, $\theta_2$, etc. are the parameters to be estimated.

The likelihood of the observed data given the parameters $\theta$ is given by:

$$L(θ) = [∏ λ(x_i; θ)] * exp(-∫ λ(x; θ) dx)$$ 

Where the product is over all observed presence locations $x_i$, and the integral is over all possible locations $x$. The first part of this formula represents the probability of observing the species at the observed locations, and the second part represents the probability of not observing the species at any other location.

The goal of maximum likelihood estimation is to find the parameters $θ$ that maximize this likelihood. This is typically done using numerical optimization methods, such as gradient descent or Newton's method.

Once we have estimated the parameters $θ$, we can use them to calculate the intensity function $λ(x; θ)$ at any location $x$. This gives us an estimate of the rate of species occurrence at that location, based on the environmental variables at that location.

### IPP Assumptions

As with any parametric statistical model, one of the biggest constraints is that there are certain assumptions that must be met in order for the model to be "valid".

-   Independence: IPP assumes that the events occur independently in space or time, meaning the occurrence of an event at one location or time does not affect the occurrence of an event at another location or time.
-   Inhomogeneity: an IPP model assumes that the intensity function can vary across space or time, as opposed to a homogeneous Poisson process, which assumes a constant rate of event occurrence.
-   Known Intensity Function: IPP models assume that the form of the intensity function is known, although the parameters of the function need to be estimated from the data. This assumption can be violated if the true intensity function is not well-captured by the chosen form.
-   Complete Spatial Coverage: IPP models assume that the entire study area has been surveyed and that presence data is available for all locations where the species is present. This assumption can be violated due to incomplete or biased sampling.

## Maximum Entropy (MaxEnt)

Entropy is a measure of uncertainty, randomness, or chaos in a set of data. In other words, it quantifies the amount of unpredictability or surprise in possible outcomes.

Maximum Entropy (MaxEnt) is a method that selects the most spread-out probability distribution fitting our known data. It's useful in presence-only problems as it minimizes bias in estimating event occurrence based on observed data. *"It agrees with everything that is known, but carefully avoids assuming anything that is not known"* (Jaynes, 1990).

### MaxEnt Assumptions

-   Incomplete Information: MaxEnt assumes that we do not have complete information about the system we are modeling. It uses the available information (presence data and environmental variables) to make the least biased estimate of the probability distribution.
-   Feature Independence: MaxEnt assumes that the features (environmental variables) are independent. In reality, this assumption is often violated as environmental variables can be correlated.
-   Linear Response: MaxEnt assumes a linear response in the log odds of presence with respect to the environmental variables. This assumption can be relaxed by including interaction terms and quadratic terms in the model.
    -   The "linear response" assumption in MaxEnt refers to the idea that the log odds of presence (i.e., the natural logarithm of the odds of a species being present versus absent) is a linear function of the environmental variables. In mathematical terms, if $P$ is the probability of presence, the log odds of presence is $log(P/(1-P))$. The linear response assumption means that this quantity is modeled as a linear function of the environmental variables. For example, if we have two environmental variables $x_1$ and $x_2$, the linear response model would be: $log(P/(1-P)) = β_0 + β_1x_1 + β_2x_2$, where $β_0$, $β_1$, and $β_2$ are parameters to be estimated from the data.
    -   This is a simplifying assumption that makes the model easier to estimate and interpret. However, it may not always be realistic. For example, the relationship between the probability of presence and the environmental variables might be nonlinear, or there might be interactions between the environmental variables.
    -   To accommodate these possibilities, we can include interaction terms and quadratic terms in the model. An interaction term (e.g., $x_1*x_2$) allows the effect of one variable to depend on the value of another variable. A quadratic term (e.g., $x_1^2$ or $x_2^2$) allows for a nonlinear relationship between the probability of presence and the environmental variable.
    -   By including these terms in the model, we can relax the linear response assumption and potentially achieve a better fit to the data. However, this also makes the model more complex and potentially harder to interpret.
-   Sample Representativeness: MaxEnt assumes that the presence data is representative of the species' distribution. This means that the locations where the species is observed are a random sample from the species' actual distribution.
