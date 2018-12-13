---
title: "An Overlooked Regression - Total Least Squares"
categories: category1
header-includes: \usepackage{mathtools}
output: html_document
layout: post
subtitle: 
tags:
- R
- regression
body-class: categoryclass
---

The source for this post is available [here](https://github.com/vladpetyuk/vladpetyuk.github.io/blob/master/blog/_R/2018-12-12-total_least_squares_regression.Rmd)


# Intro & Motivation
When someone needs to perform a linear regression, the first thing that comes to mind is `lm` function from the the `stats` R package. However, this linear model relies on assumption that the independent variable `x` does not contain any measurement error. In case of categorical variables it is probably OK. However, it is still possible that certain points belong to a certain class or factor level with some degree of probability. I won't consider this case here. Let's take a look at the situation when `x` is a continuous variable. If it is a real-world scenario then most likely `x` is measured with some error. Out of my experience I can point to an example of correlating mRNA vs protein relative abundance. Or correlating protein relative abundance vs some clinical measure. Both values have decent (or at least non-negligible) degree of measurement error. I emphasize *both*.

# Web Resources
A couple of good discussions on that topic can be found [here](https://stats.stackexchange.com/questions/13152/how-to-perform-orthogonal-regression-total-least-squares-via-pca) and [here](https://stackoverflow.com/questions/6872928/how-to-calculate-total-least-squares-in-r-orthogonal-regression). 

## Simulation
For simulating the example below I picked t-distribution to represent relative protein/mRNA abundances. Typically they don't just follow normal. Normal component describe mostly measurement noise. Some proteins/mRNA trully change. Thus there is a need for higher then normal probability in the tails. Alternatively one can simulate relative protein abundances with a two-component distribution - combination of normal and uniform.

{% highlight r %}
set.seed(0)
x0 <- rt(2000, 10)
x <- rnorm(length(x0), x0, 1)
y <- rnorm(length(x0), x0, 1)
plot(x,y)
{% endhighlight %}

![figure](/blog/figs/2018-12-12-total_least_squares_regression/unnamed-chunk-1-1.png)

## Regular linear model
Regular linear model `lm` relies on oridinary least squares approach. It assumed that `x` does not have any error. The goodness of fit calculated as squared difference between observed and fitted `y` values.

{% highlight r %}
m1 <- lm(y ~ x)
plot(x,y)
abline(m1, col="red")
{% endhighlight %}

![figure](/blog/figs/2018-12-12-total_least_squares_regression/unnamed-chunk-2-1.png)
The slope of that line doesn't look right. Isn't it?

## Total least squares leveraging PCA
It is fairly intuitive that PCA can be one helpful approach here. The most variance is along the `x vs y` slope. Thus PCA will rotate the scatterplot such that first principal component will be along the slope.

{% highlight r %}
v <- prcomp(cbind(x, y))$rotation
beta <- v[2,1]/v[1,1]

plot(x,y)
abline(m1, col="red")
abline(a = 0, b = beta, col="green")
{% endhighlight %}

![figure](/blog/figs/2018-12-12-total_least_squares_regression/unnamed-chunk-3-1.png)

## TLS regression `pracma::odregress`
`pracma` package has a dedicated function for total least squares or alternatively called orthogonal distance regression.

{% highlight r %}
library(pracma)
m2 <- odregress(x, y)

plot(x,y)
abline(m1, col="red")
abline(a = 0, b = beta, col="green")
abline(b=m2$coeff[1], a=m2$coeff[2], col="cyan")
{% endhighlight %}

![figure](/blog/figs/2018-12-12-total_least_squares_regression/unnamed-chunk-4-1.png)


## Deming regression
The most general approach for hangling such scenarios would be [Deming regression](https://en.wikipedia.org/wiki/Deming_regression). In its simplest form with default settings it is equivalent to total least squares.

{% highlight r %}
library(MethComp)
{% endhighlight %}



{% highlight text %}
## Loading required package: nlme
{% endhighlight %}



{% highlight r %}
m3 <- Deming(x, y)

plot(x,y)
abline(m1, col="red")
abline(a = 0, b = beta, col="green")
abline(b=m2$coeff[1], a=m2$coeff[2], col="cyan")
abline(a=m3["Intercept"], b=m3["Slope"], col="yellow", lwd=2)
{% endhighlight %}

![figure](/blog/figs/2018-12-12-total_least_squares_regression/unnamed-chunk-5-1.png)



