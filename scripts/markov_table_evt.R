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
weights <- c(0.6, 0.399, 0.001)
# Range of mixture
x_seq <- seq(1, 900, length.out = 1000)

# COmpute mixture ccdf
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_ccdf <- c(mixture_ccdf, 1 - (
    weights[1] * pnorm(x, mean = mean1, sd = sd1) +
      weights[2] * pnorm(x, mean = mean2, sd = sd2) +
      weights[3] * pnorm(x, mean = mean3, sd = sd3)
  ))
}
plot(x_seq, mixture_ccdf, log = "y")

mixture_data <- tibble(ccdf = mixture_ccdf, x = x_seq)

# Get a big sample to compute the expected value^k
mostra <- sort(c(rnorm(weights[1]*10000000, mean = mean1, sd = sd2), 
                 rnorm(weights[2]*10000000, mean = mean2, sd = sd2),
                 rnorm(weights[3]*10000000, mean = mean3, sd = sd3)))

# Threshold from where to fit the exponential
#thresholds <- c(76.76,123.14,142.59)
thresholds <- c(69.39, 107, 141)

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
  lambda <- 1 - (weights[1] * pnorm(trunc, mean = mean1, sd = sd1) +
                   weights[2] * pnorm(trunc, mean = mean2, sd = sd2) +
                   weights[3] * pnorm(trunc, mean = mean3, sd = sd3))
  x_theo <- seq(trunc, max_dist, length.out = 1000)
  
  # For high quantiles, get the truncated normal mean
  if (trunc > 400) {
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

which(fit_data_3$ccdf_mostra_trunc > 1.2e-12 && fit_data_3$ccdf_mostra_trunc < 9e-13)

868/699

# Weibull Markov ----------------------------------------------------------

library(ismev)
library(ercv)
library(evir)

# Make mixture
location <- 40000
scale <- 10
shape <- 4

x_seq <- seq(from = 0, to = 20, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pweibull(x, shape = shape, scale = scale))
}
plot(-x_seq+20, uni_ccdf, log = "y")


fit_gpd <- fitpot(a, evi = -shape, threshold = 0)
theo_gpd <- evir::qgpd(uni_ccdf, 
                       xi = -shape, 
                       beta = scale,
                       mu = 0)

plot(x = x_seq+200, y = uni_ccdf, log= "y")
lines(x = theo_gpd, y = rev(uni_ccdf), col = "red")


# Weibull 10-3 ------------------------------------------------------------
p_seq <- lseq(0.001, 1-10^(-12), length.out  = 100000)

threshold <- evir::qgev(p = p_seq, xi = -shape, sigma = scale)



a_thresh <- a[a > threshold]/threshold
fit_gpd <- fitpot(a_thresh, evi = -shape, threshold = 1)
theo_gpd <- evir::qgpd(uni_ccdf, 
                       xi = -shape, 
                       beta = fit_gpd$coeff["psi"],
                       mu = 1)

theo_gpd_re <- theo_gpd*threshold 
plot(x = x_seq, y = uni_ccdf, log= "y")
lines(x = theo_gpd_re, y = 1-uni_ccdf, col = "red")





a <- evir::qgev(p=uni_ccdf, xi = -shape, sigma = scale)
a <- a[a!=Inf]
plot(x = a, y = 1-uni_ccdf, log= "y")


threshold <- qweibull(1-10^-6, shape = 1/shape, scale = scale)
a_thr <- a[a>threshold]/threshold

lambda <- pweibull(threshold, shape = 1/shape, scale = scale)
thrselect(a, threshold = 130)

fit_gpd <- fitpot(a_thr, evi = -shape, threshold = 1)
probs <- uni_ccdf[uni_ccdf < 10^-6]
probs <- probs[probs > 0]
prob_lambda <- ((1-probs) - lambda)/( 1 - lambda)

theo_gpd <- evir::qgpd(prob_lambda, 
                       xi = -shape, 
                       beta = fit_gpd$coeff["psi"],
                       mu = 1)
quants <- (threshold*theo_gpd)

plot(x = a, y = 1-uni_ccdf, log= "y")
lines(x = quants, y = (1-prob_lambda), col = "red")

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


# Beta --------------------------------------------------------------------

# Make mixture
shape1 <- 1/8
shape2 <- 8

x_seq <- seq(from = 0, to = 1, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pbeta(x, shape1 = shape1, shape2 = shape2))
}
plot(x_seq, uni_ccdf, log = "y")


p <- 0.5
thresh <- qbeta(p = p, shape1 = shape1, shape2 = shape2)
lambda <- pbeta(q = thresh, shape1 = shape1, shape2 = shape2)
#shape2 <- 1/8
#shape1 <- 1/8



x_thresh <- x_seq[x_seq > thresh] 
theo_gpd <- evd::pgpd(x_thresh, 
                      shape = -1/shape2, 
                      scale = (1-thresh)/shape2,
                      loc = thresh)

#save plots log and no log
plot(x_seq, uni_ccdf, log = "y")
lines(x_thresh , (1-lambda)*(1-theo_gpd),col = "red")

p_test <- c(1-10^(-6),  1-10^(-9), 1-10^(-12))

q_fit <- evd::qgpd(p_test, 
          shape = -1/shape2, 
          scale = (1-thresh)/shape2,
          loc = thresh)

q_true <- qbeta(p = p_test, shape1 = shape1, shape2 = shape2)

q_fit/q_true


# Weibull -----------------------------------------------------------------



# Make mixture
shape <- 8
scale <- 80

p_seq <- lseq(from = 0.00000001, to = 1-10^(-14), length.out = 500)
x_weib <- c()
for (p in p_seq){
  x_weib  <- c(x_weib , qweibull(p, shape = shape, scale = scale))
}
plot(x_weib, 1-p_seq, log = "y")


p <- 0.5
thresh <- qweibull(p=p, shape = shape, scale = scale)
lambda <- pweibull(q=thresh, shape = shape, scale = scale)

p_test <- c(1-10^(-6),  1-10^(-9), 1-10^(-12))

q_fit <- evd::qgpd(p_test, 
                   shape = -1/shape, 
                   scale = (max(x_weib)-thresh)/shape,
                   loc = thresh)

q_true <- qweibull(p = p_test, shape = shape, scale = scale)

q_fit/q_true

