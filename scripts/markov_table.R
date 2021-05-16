library(truncnorm)
library(tidyverse)
library(data.table)
library(evd)
source("scripts/markov_functions.R")
# Weibull Markov ----------------------------------------------------------

# Make mixture
location <- 40000
scale <- 80
shape <- 8

# kth moment of the weibull distribution
weibull_k_moment <- function(k, scale, shape){
  (scale^k)*gamma(1 + k/shape)
}

x_seq <- -seq(from = 0, to = 200, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - evd::prweibull(x, shape = shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

p_test <- 1-1e-6
x_test <- qweibull(p_test, shape = shape, scale = scale)

ks <- seq(from = 5, to = 200, by = 1)
iter <- 1
data_test_k <- c()
for (k in ks){
  e_k <- weibull_k_moment(k = k, shape = shape, scale = scale)
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
contour_quantile <-  (weibull_k_moment(k = best_k, shape = shape, scale = scale)^(1/best_k) / (1-p_test)^(1/best_k))
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



# Beta  -------------------------------------------------------------------


# Make mixture
location <- 40000
scale <- 100
shape1 <- 1/8
shape2 <- 8
# kth moment of the weibull distribution
beta_k_moment <- function(k, sh1, sh2){
  ev_k <- c()
  ev_k[1] <- sh1/(sh1 + sh2)
  if(k > 1){
    for(i in 2:k){
      ev_k[i] <- ev_k[i - 1]*(sh1 + i - 1)/(sh1 + sh2 + i - 1)
    }
  }
 return(tail(ev_k, n = 1))
}

x_seq <- seq(from = 0, to = 1, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pbeta(x, shape1 = shape1, shape2 = shape2))
}
plot(x_seq, uni_ccdf, log = "y")

p_test <-  c(1-10^(-6))
x_test <- qbeta(p_test, shape1 = shape1, shape2 = shape2)

ks <- seq(from = 5, to = 200, by = 1)
iter <- 1
data_test_k <- c()
for (k in ks){
  e_k <- beta_k_moment(k = k, sh1 = shape1, sh2 = shape2)
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
contour_quantile <-  (beta_k_moment(k = best_k, sh1 = shape1, sh2 = shape2))^(1/best_k) / ((1-p_test)^(1/best_k))
x_test
contour_quantile/x_test





# Gamma -------------------------------------------------------------------



# Make mixture
location <- 40000
shape <- 1/8
scale  <- 100
# kth moment of the weibull distribution
gamma_k_moment <- function(k, sh, sc){
  (sc^k)*(gamma(k + sh))/(gamma(sh))
}

x_seq <- seq(from = 0, to = 4000, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pgamma(x, shape = shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

p_test <- 1-1e-12
x_test <- qgamma(p_test, shape = shape, scale = scale)

ks <- seq(from = 5, to = 200, by = 1)
iter <- 1
data_test_k <- c()
for (k in ks){
  e_k <- gamma_k_moment(k = k, sh = shape, sc = scale)
  cota_k <-  e_k / x_test^(k)
  names(cota_k) <- paste0("k", k)
  data_test_k[iter] <- cota_k
  iter <- iter + 1
  lines(x_seq,gamma_k_moment(k,shape,scale)/(x_seq^k))
  #lines(x_seq, cota_k, col = colors_palette[k], lwd = 2)
}


data_test_k <- na.omit(data_test_k)
data_test_k <- data_test_k[data_test_k!=0]
data_test_k <- data_test_k[data_test_k!=Inf]
index_min_k <- which(data_test_k==min(data_test_k))
best_k <- ks[index_min_k]
contour_quantile <-  (gamma_k_moment(k = best_k, sh = shape, sc = scale))^(1/best_k) / ((1-p_test)^(1/best_k))
x_test
contour_quantile/x_test



# Gaussian mixture --------------------------------------------------------

source("scripts/markov_functions.R")
library(CharFun)
# Make mixture
mean1 <- 5
mean2 <- 50
mean3 <- 100
sd1 <- 10
sd2 <- 10
sd3 <- 10

x_seq <- seq(1, 170, 1)
mixture_density <- c()
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_density <- c(
    mixture_density,
    0.6 * dnorm(x, mean = mean1, sd = sd1) +
      0.399 * dnorm(x, mean = mean2, sd = sd2) +
      0.001 * dnorm(x, mean = mean3, sd = sd3)
  )
  mixture_ccdf <- c(mixture_ccdf, 1 - (
    0.6 * pnorm(x, mean = mean1, sd = sd1) +
      0.399 * pnorm(x, mean = mean2, sd = sd2) +
      0.001 * pnorm(x, mean = mean3, sd = sd3)
  ))
}

plot_data <- tibble(theo =  mixture_ccdf, x = x_seq)

p_test <- 1-1e-12
x_test <- 160

iter <- 1
data_test_k <- c()
ks <- seq(from = 5, to = 200, by = 1)
iter <- 1
for (k in ks){
  e_k <- 0.6 * (tail(gaussian_expected_value_k(k+1, mean1, sd1), n = 1) )/(1-pnorm(0, mean = mean1, sd = sd1)) + 
    0.399 * (tail(gaussian_expected_value_k(k+1,mean2, sd2), n = 1))/(1-pnorm(0, mean = mean2, sd = sd2)) + 
    0.001 * (tail(gaussian_expected_value_k(k+1, mean3, sd3), n = 1))/(1-pnorm(0, mean = mean3, sd = sd3))
  cota_k <-  e_k / (x_test)^(k)
  data_test_k[iter] <- cota_k
  iter <- iter + 1
}


data_test_k <- na.omit(data_test_k)
data_test_k <- data_test_k[data_test_k!=0]
data_test_k <- data_test_k[data_test_k!=Inf]
index_min_k <- which(data_test_k==min(data_test_k))
best_k <- ks[index_min_k]
contour_quantile <-  (0.6 * (tail(gaussian_expected_value_k(best_k+1, mean1, sd1), n = 1) )/(1-pnorm(0, mean = mean1, sd = sd1)) + 
                        0.399 * (tail(gaussian_expected_value_k(best_k+1,mean2, sd2), n = 1))/(1-pnorm(0, mean = mean2, sd = sd2)) + 
                        0.001 * (tail(gaussian_expected_value_k(best_k+1, mean3, sd3), n = 1))/(1-pnorm(0, mean = mean3, sd = sd3)))^(1/best_k) / ((1-p_test)^(1/best_k))
x_test
contour_quantile/x_test
