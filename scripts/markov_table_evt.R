library(truncnorm)
library(tidyverse)
library(data.table)
source("scripts/markov_functions.R")


# Gaussian mixture ---------------------------------------------------------


# Make mixture
mean1 <- 5
mean2 <- 50
mean3 <- 100
sd1 <- 10
sd2 <- 10
sd3 <- 10
weights <- c(0.6, 0.399, 0.001)
# Range of mixture
x_seq <- seq(1, 170, length.out = 1000)

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
thresholds <- c(14.53, 46.67, 56.65)

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
  x_theo <- seq(trunc, max_dist+1000, length.out = 1000)
  
  # For high quantiles, get the truncated normal mean

    mostra_trunc <- mostra[mostra > trunc] - trunc
    # The rate for the exp is 1/mean
    rate_mostra_trunc <- 1 / mean(mostra_trunc)

  
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
lines(x_thresh , (1-lambda)*(1-theo_gpd),col = "red", lwd = 2)

p_test <- c(1-10^(-6),  1-10^(-9), 1-10^(-12))

q_fit <- evd::qgpd(p_test, 
          shape = -1/shape2, 
          scale = (1-thresh)/shape2,
          loc = thresh)

q_true <- qbeta(p = p_test, shape1 = shape1, shape2 = shape2)

q_fit/q_true


# Weibull -----------------------------------------------------------------


library(evd)
library(ercv)

# Make mixture
shape <- 4
scale <- 80

p_seq <- seq(from = 0, to = 1-10^(-13), length.out = 100000)
x_weib <- c()
for (p in p_seq){
  x_weib  <- c(x_weib , qweibull(p, shape = shape, scale = scale))
}
x_rweib <- rev(-x_weib + max(x_weib))

#plot(x_rweib, 1-p_seq, log = "y")
#hist(x_rweib)

rweib_data <- tibble(x = x_rweib, ccdf = 1-p_seq)



# Make mixture
shape <- 4
scale <- 80

x_seq <- seq(from = 0, to = 100, length.out = 1000)
p_weib <- c()
for (x in x_seq){
  p_weib  <- c(p_weib , pweibull(x, shape = shape, scale = scale))
}

a <- -x_seq + max(x_seq)
plot(a, p_weib, log = "y")
rweib_data <- tibble(x = a, ccdf = p_weib)


p_fit <- evd::qgpd(p_weib[-1], 
                   shape = -1/shape, 
                   scale = (max(a)-53.45)/shape,
                   loc = 0)


plot(a, p_weib, log = "y")
lines(a[-1], p_fit, col = "red")


results <- tibble(x = a[-1], fit = p_fit)




#samp <- rweibull(n = 10^5, shape = shape, scale = scale)

#cvplot(-samp +200 , evi = -1/shape)

threshold <- qweibull(0.9, shape = shape, scale = scale)
lambda <- 1-pweibull(threshold , shape = shape, scale = scale)
p_fit <- evd::pgpd(x_rweib, 
                   shape = -1/shape, 
                   scale = (max(x_rweib)-threshold)/shape,
                   loc = threshold)


plot(x_rweib, 1-p_seq)
lines(x_rweib, (1-lambda)*(1-p_fit), col = "red")

p_test <- c(1-10^(-6),  1-10^(-9), 1-10^(-12))
threshold <- qweibull(0.90, shape = shape, scale = scale)
q_fit <- evd::qgpd(p_test, 
                   shape = -1/shape, 
                   scale = (max(x_rweib)-40)/shape,
                   loc = 40)


q_true <- 64

q_fit/q_true


# Rev Weibull sampl -------------------------------------------------------


library(evd)
library(ercv)

# Make mixture
shape <- 4
scale <- 80

sample_w <- rweibull(10^4, shape = shape, scale = scale)
sample_rw <- sort(-sample_w + max(sample_w ))
hist(sample_rw)

fit <- fitpot(sample_rw)
weib_ccdf <- ecdf(sample_rw)


threshold <- 0
lambda <-  1-weib_ccdf(threshold) 
p_fit <- evd::pgpd(sample_rw, 
                   shape = -1/shape, 
                   scale = (max(sample_w)-threshold)/shape,
                   loc = threshold)


plot(sample_rw, 1-weib_ccdf(sample_rw))
lines(sample_rw,(1-p_fit), col = "red")


q_fit <- evd::qgpd(p_test, 
                   shape = -1/shape, 
                   scale = (max(x_weib)-thresh)/shape,
                   loc = thresh)

q_true <- qweibull(p = p_test, shape = shape, scale = scale)

q_fit/q_true

# Gaussian ----------------------------------------------------------------


library(truncnorm)
# Make mixture
mean <- 100
sd <- 50


# Range of mixture
x_seq <- seq(1, 1000, length.out = 1000)

# COmpute mixture ccdf
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_ccdf <- c(mixture_ccdf, 1 - 
   pnorm(x, mean = mean, sd = sd)) 
}

# Get a big sample to compute the expected value^k
mostra <- sort(c(rnorm(10000000, mean = mean, sd = sd)))

# Initiate pdf
plot_data <- tibble(theo =  mixture_ccdf, x = x_seq)

thresholds <- qnorm(p = c(0.5, 0.75, 0.90), mean = mean, sd = sd)
# compute Markov's
for (i in 3:3) {
  max_dist <- max(x_seq)
  trunc <- thresholds[i]
  
  # Get lambda to scale the probability to the threshold
  lambda <- 1 - pnorm(trunc, mean = mean, sd = sd)
  x_theo <- seq(trunc, max_dist+500, length.out = 10000)
  
  mostra_trunc <- mostra[mostra > trunc] - trunc
  # The rate for the exp is 1/mean
  rate_mostra_trunc <- 1 / mean(mostra_trunc)

  # Compute the ccdf from the threshold
  ccdf_mostra_trunc <- as.tibble(lambda * (1 - pexp(x_theo - trunc, rate_mostra_trunc)))
  plot_data<- cbind(plot_data, x_theo, ccdf_mostra_trunc)
}

plot(x = plot_data$x, y = plot_data$theo, log = "y")
lines(x = plot_data$x_theo, y = plot_data$value, col = "red")


p_test <- c(1-10^(-6),  1-10^(-9), 1-10^(-12))

q_fit <- (1-lambda)*qexp(p_test, rate = trunc_rate) 

q_true <-  qnorm(p=p_test, mean = mean, sd = sd)

q_fit/q_true



# Weibull mixture ---------------------------------------------------------


# Make mixture
scale1 <- 5
scale2 <- 50
scale3 <- 100
shape1 <- 8
shape2 <- 8
shape3 <- 8
weights <- c(0.6, 0.399, 0.001)
# Range of mixture
x_seq <- seq(0, 200, length.out = 10000)

# COmpute mixture ccdf
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_ccdf <- c(mixture_ccdf, (
    weights[1] * pweibull(x, shape = shape1, scale = scale1) +
      weights[2] * pweibull(x, shape = shape2, scale = scale2) +
      weights[3] * pweibull(x, shape = shape3, scale = scale3)
  ))
}


x_weib <- -x_seq +max(x_seq) 
plot(x_weib, mixture_ccdf, log = "y")
weib_data <- tibble(x = x_weib, ccdf = mixture_ccdf)

# Threshold from where to fit the exponential

p_test <- c(1-10^(-6),  1-10^(-9), 1-10^(-12))

# Make plot empty to draw the line later
plot(x_weib, mixture_ccdf, log = "y" , xlab = "Time", ylab = "CCDF")

max_dist <-max(x_seq) 
trunc <- 195.35

# Get lambda to scale the probability to the threshold
lambda <- (weights[1] * pweibull(trunc, shape = shape1, scale = scale1) +
               weights[2] * pweibull(trunc, shape = shape2, scale = scale2) +
               weights[3] * pweibull(trunc, shape = shape3, scale = scale3))
x_theo <- seq(trunc, max_dist, length.out = 100000)

mix_data <- tibble(x = x_weib, ccdf = mixture_ccdf) %>% 
  filter(mixture_ccdf != 0, mixture_ccdf != 1) %>% 
  filter(x > trunc)

fit <- evd::pgpd(mix_data$x, 
          shape = -1/shape1, 
          scale = (max_dist-trunc)/shape1,
          loc = trunc)

lines(x = mix_data$x, y =1- fit , col = "red")

mix_data$fit <- 1-fit


fit_data_3 <- fit_data[[3]]


# Beta mixture ------------------------------------------------------------


# Make mixture
shape1 <- 4
shape2 <- 8
scale1 <- 5
scale2 <- 10
scale3 <- 20
weights <- c(0.6, 0.399, 0.001)
# Range of mixture
x_seq <- seq(0, 1, length.out = 1000)

# COmpute mixture ccdf
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_ccdf <- c(mixture_ccdf, (
    weights[1] * pbeta(x, shape1 =  4, shape2 = 4) +
      weights[2] * pbeta(x, shape1 = 80, shape2 = 80) +
      weights[3] * pbeta(x, shape1 = 160, shape2 = 160)
  ))
}

plot(x_seq, 1-mixture_ccdf, log = "y")
weib_data <- tibble(x = x_weib, ccdf = mixture_ccdf)

# Threshold from where to fit the exponential

p_test <- c(1-10^(-6),  1-10^(-9), 1-10^(-12))

# Make plot empty to draw the line later
plot(x_seq, 1-mixture_ccdf, log = "y" , xlab = "Time", ylab = "CCDF")

max_dist <-max(x_seq) 
trunc <- 0.1

mix_data <- tibble(x = x_seq, ccdf = mixture_ccdf) %>% 
  filter(mixture_ccdf != 0, mixture_ccdf != 1) %>% 
  filter(x > trunc)

fit <- evd::qgpd( mix_data$ccdf, 
                  shape = -1/160, 
                  scale = (1-trunc)/160,
                  loc = trunc)

lines(x = mix_data$x, y = 1-mix_data$ccdf, col = "red")


