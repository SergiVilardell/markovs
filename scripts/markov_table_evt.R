library(truncnorm)
library(tidyverse)
library(data.table)
source("scripts/markov_functions.R")
# Make mixture
mean1 <- 50
mean2 <- 100
mean3 <- 400
sd1 <- 50
sd2 <- 50
sd3 <- 50

# Range of mixture
x_seq <- seq(1, 900, length.out = 1000)

# COmpute mixture ccdf
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_ccdf <- c(mixture_ccdf, 1 - (
    0.75 * pnorm(x, mean = mean1, sd = sd1) +
      0.2499 * pnorm(x, mean = mean2, sd = sd2) +
      0.0001 * pnorm(x, mean = mean3, sd = sd3)
  ))
}
plot(x_seq, mixture_ccdf, log = "y")

mixture_data <- tibble(ccdf = mixture_ccdf, x = x_seq)

# Get a big sample to compute the expected value^k
mostra <- sort(c(rnorm(75000000, mean = mean1, sd = sd2), rnorm(24990000, mean = mean2, sd = sd2), rnorm(10000, mean = mean3, sd = sd3)))

# Threshold from where to fit the exponential
#thresholds <- c(76.76,123.14,142.59)
thresholds <- c(234.97, 554.43, 612.93)
# Colors for the plot
colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(length(thresholds))

# Make plot empty to draw the line later
plot(x_seq, mixture_ccdf, log = "y", xlab = "Time", ylab = "CCDF", type = "n")
fit_data <- list()
# compute Markov's
for (i in 1:length(thresholds)) {
  max_dist <- max(x_seq)
  trunc <- thresholds[i]
  
  # Get lambda to scale the probability to the threshold
  lambda <- 1 - (0.75 * pnorm(trunc, mean = mean1, sd = sd1) +
                   0.2499 * pnorm(trunc, mean = mean2, sd = sd2) +
                   0.0001 * pnorm(trunc, mean = mean3, sd = sd3))
  x_theo <- seq(trunc, max_dist+300, length.out = 1000)
  
  # For high quantiles, get the truncated normal mean
  if (trunc > 500) {
    rate_mostra_trunc <- 1 / (etruncnorm(a = trunc, b = Inf, mean = mean3, sd = sd3) - trunc)
  } else {
    mostra_trunc <- mostra[mostra > trunc] - trunc
    # The rate for the exp is 1/mean
    rate_mostra_trunc <- 1 / mean(mostra_trunc)
  }
  
  # Compute the ccdf from the threshold
  ccdf_mostra_trunc <- lambda * (1 - pexp(x_theo - trunc, rate_mostra_trunc))
  lines(x_theo, ccdf_mostra_trunc, col = colors_palette[i])
  fit_data[[i]] <- tibble(ccdf_mostra_trunc, x_theo)
}

# Draw the distribution line now so it is on top of the others
points(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF", lwd = 4)

fit_data_3 <- fit_data[[3]]





# Weibull Markov ----------------------------------------------------------

library(ismev)
# Make mixture
location <- 40000
scale <- 100
shape <- 1/8

# kth moment of the weibull distribution
weibull_k_moment <- function(k, scale, shape){
  (scale^k)*gamma(1 + k/shape)
}

x_seq <- seq(from = 1, to = 400, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pweibull(x, shape = 1/shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

a <- qweibull(uni_ccdf, shape = 1/shape, scale = scale)
a <- a[a!=Inf]
plot(x = a, y = 1-uni_ccdf, log= "y")

library(ercv)
ercv::
gev.fit(a)

p_test <- 1-1e-12
x_test <- qweibull(p_test, shape = 1/shape, scale = scale)

ks <- seq(from = 5, to = 200, by = 1)
iter <- 1
data_test_k <- c()
for (k in ks){
  e_k <- weibull_k_moment(k = k, shape = 1/shape, scale = scale)
  cota_k <-  e_k / x_test^(k)
  names(cota_k) <- paste0("k", k)
  data_test_k[iter] <- cota_k
  iter <- iter + 1
  #lines(x_seq, cota_k, col = colors_palette[k], lwd = 2)
}


