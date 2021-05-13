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

x_seq <- seq(from = 0, to = 200, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pweibull(x, shape = 1/shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

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

data_test_k <- na.omit(data_test_k)
data_test_k <- data_test_k[data_test_k!=0]
data_test_k <- data_test_k[data_test_k!=Inf]
index_min_k <- which(data_test_k==min(data_test_k))
best_k <- ks[index_min_k]
contour_quantile <-  (weibull_k_moment(k = best_k, shape = 1/shape, scale = scale)^(1/best_k) / (1-p_test)^(1/best_k))
x_test
contour_quantile/x_test



# Gumbel ------------------------------------------------------------------


# Make mixture
location <- 40000
scale <- 1


x_seq <- seq(from = 0, to = 60, length.out = 500)
#x_seq <- lseq(from = 1, to = qfrechet(1- 1e-1, shape = shape, scale = scale), length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - evd::pgumbel(x, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

p_test <- 1-1e-12
x_test <- evd::qgumbel(p_test, scale = scale)

data_test_k <- c()
ks <- 1:50
for (k in ks){
  e_k <- gamma(1 + k)
  cota_k <-  e_k / exp(x_test*k)
  names(cota_k) <- paste0("k", k)
  data_test_k[iter] <- cota_k
  iter <- iter + 1
}

data_test_k <- na.omit(data_test_k)
data_test_k <- data_test_k[data_test_k!=0]
data_test_k <- data_test_k[data_test_k!=Inf]
index_min_k <- which(data_test_k==min(data_test_k))
best_k <- ks[index_min_k]


contour_quantile <-  log((gamma(1 + best_k)) / (1-p_test))/best_k
x_test
contour_quantile/log(x_test)



# Gaussian ----------------------------------------------------------------


# Make mixture
mean1 <- 100
sd1 <- 50

p_test <- 1-1e-12
x_test <- qnorm(p_test, mean = mean1, sd = sd1)

data_test_k <- c()
ks <- 1:200
for (k in ks){
  e_k <- tail(gaussian_expected_value_k(k + 1, mean1, sd1), n = 1)
  cota_k <-  e_k / x_test^k
  names(cota_k) <- paste0("k", k)
  data_test_k[iter] <- cota_k
  iter <- iter + 1
}

data_test_k <- na.omit(data_test_k)
data_test_k <- data_test_k[data_test_k!=0]
data_test_k <- data_test_k[data_test_k!=Inf]
index_min_k <- which(data_test_k==min(data_test_k))
best_k <- ks[index_min_k]


contour_quantile <-  (tail(gaussian_expected_value_k(best_k + 1, mean1, sd1), n = 1)^(1/best_k)) / ((1-p_test)^(1/best_k))
x_test
contour_quantile/x_test



# Frechet -----------------------------------------------------------------

# Make mixture
location <- 40000
scale <- 100
shape <- 1/8


p_test <- 1-1e-9
x_test <- evd::qfrechet(p_test, shape = 1/shape, scale = scale)


ks <- 1:(1/shape-1)
data_test_k <- c()
iter <- 1
for (k in ks){
  e_k <- reverse_weibull_k_moment(k = k, scale = scale, shape = 1/shape)
  cota_k <-  e_k / (x_test)^(k)
  data_test_k[iter] <- cota_k
  iter <- iter + 1
}


data_test_k <- na.omit(data_test_k)
data_test_k <- data_test_k[data_test_k!=0]
data_test_k <- data_test_k[data_test_k!=Inf]
index_min_k <- which(data_test_k==min(data_test_k))
best_k <- ks[index_min_k]


contour_quantile <-  (reverse_weibull_k_moment(k = best_k, scale = scale, shape = 1/shape))^(1/best_k) / (1-p_test)^(1/best_k)
x_test
contour_quantile/x_test

