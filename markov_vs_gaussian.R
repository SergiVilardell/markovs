library(mixtools)
library(tidyverse)

# Here we will compare how tight MI when changing Gaussian parameters

m <- 100
sd <- 1
gaussian_narrow <- sort(rnorm(1000, mean = m, sd = sd))
translation <- min(gaussian_narrow)
gaussian_narrow_transf <- gaussian_narrow - translation

# compute the expected value
ev <- m - translation
p_a <- c()
markov <- c()
# get a random high quantile
for (i in 1:length(gaussian_narrow_transf)) {
  p_a[i] <- 1 - pnorm(gaussian_narrow_transf[i], mean = m - translation, sd = sd)
  markov[i] <- ev / gaussian_narrow_transf[i]
}

plot(y = p_a, x = gaussian_narrow_transf)
lines(y = markov, x = gaussian_narrow_transf)
