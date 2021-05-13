library(truncnorm)
library(tidyverse)
library(data.table)
source("scripts/markov_functions.R")
# Weibull Markov ----------------------------------------------------------

# Make mixture
location <- 40000
scale <- 100
shape <- 1/4

# kth moment of the weibull distribution
weibull_k_moment <- function(k, scale, shape){
  (scale^k)*gamma(1 + k/shape)
}

p_test <- 1-1e-12
x_test <- qweibull(p_test, shape = 1/shape, scale = scale)

contour_quantile <-  (weibull_k_moment(k = 1, shape = 1/shape, scale = scale) / (1-p_test))
x_test
contour_quantile/x_test



# Gumbel ------------------------------------------------------------------


# Make mixture
location <- 40000
scale <- 1

p_test <- 1-1e-12
x_test <- evd::qgumbel(p_test, scale = scale)


contour_quantile <-  log((gamma(1 + 1)) / (1-p_test))
x_test
contour_quantile/log(x_test)




# Frechet -----------------------------------------------------------------

# Make mixture
location <- 40000
scale <- 100
shape <- 1/8


p_test <- 1-1e-12
x_test <- evd::qfrechet(p_test, shape = 1/shape, scale = scale)


contour_quantile <-  (reverse_weibull_k_moment(k = 1, scale = scale, shape = 1/shape)) / (1-p_test)
x_test
contour_quantile/x_test




