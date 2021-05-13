source("scripts/markov_functions.R")
library(truncnorm)
library(tidyverse)
library(data.table)

# Gaussian ----------------------------------------------------------------


# Make mixture
mean1 <- 100
mean2 <- 50
mean3 <- 100
sd1 <- 10
sd2 <- 20
sd3 <- 30



x_seq <- seq(from = 1,  to = 200, length.out = 100)
mixture_ccdf <- c()
for (x in x_seq){
  mixture_ccdf <- c(mixture_ccdf, 1 - pnorm(x, mean = mean1, sd = sd1))
}



# pdf("markov_mixture.pdf", height = 6, width = 12)
# par(mfrow=c(1,2))
# plot(x_seq, mixture_density, xlab = "Time", ylab = "Prob Density")
# plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
# lines(x_seq, cota, col = "red")
#  dev.off()


colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(100)

pdf("markov_guassian_k_multi.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- c(10,50,90)
plot_data <- tibble(theo =  mixture_ccdf, x = x_seq)
for (k in ks){
  e_k <- tail(gaussian_expected_value_k(k + 1, mean1, sd1), n = 1)
  cota_k <-  as.tibble(e_k / x_seq^(k))
  color_k <- rep(colors_palette[k], nrow(cota_k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
  #lines(x_seq, cota_k, col = colors_palette[k], lwd = 2)
}


legend(25, 1e-9, legend=c("k = 10", "k = 50", "k = 90"),
       col=c(colors_palette[10], colors_palette[50],colors_palette[90]), lty=1, cex=1, lwd = 2)
dev.off()





pdf("markov_k_200.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, log = "y", xlab = "Time", ylab = "CCDF", type = "n")
ks <- 1:200
for (k in ks){
  e_k <-truncatedNormalMoment(k = k, x_lower = 0, x_upper = Inf, mean = mean1, sd = sd1)
  cota_k <-  e_k / (x_seq)^(k)
  lines(x_seq, cota_k, col = colors_palette[k])
}
points(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF", lwd = 4)


pdf("markov_k_200.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, log = "y", xlab = "Time", ylab = "CCDF", type = "n")
ks <- 1:200
for (k in ks){
  e_k <-  gaussian_abs_expected_value_k(k = k,  mu = mean1, sigma = sd1)
  cota_k <-  e_k / (x_seq)^(k)
  lines(x_seq, cota_k, col = colors_palette[k])
}
points(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF", lwd = 4)



dev.off()
lines(xs, ys, col = "red")



pdf("markov_k_100.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, log = "y", xlab = "Time", ylab = "CCDF", type = "n")
ks <- 1:100
mostra <- c(rnorm(6000000, mean = mean1, sd = sd1), rnorm(3000000, mean = mean1, sd = sd1), rnorm(1000000, mean = mean1, sd = sd1))
for (k in ks) {
  cotak <- mean(mostra^k) / x_seq^k
  lines(x_seq, cotak, col = colors_palette[k])
}
points(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF", lwd = 4)



# Weibull -----------------------------------------------------------------

library(evd)




# Reverse weibull ---------------------------------------------------------


# Make mixture
location <- 40000
scale <- 1
shape <- 1/2


x_seq <- -seq(from = 0, to = 10, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - evd::prweibull(x, shape = shape, scale = scale))
}
plot(-x_seq, rev(uni_ccdf), log = "y")

colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(100)

x_seq_2 <- -x_seq
#pdf("markov_k_weibull_4.pdf", height = 4, width = 6)
plot(x_seq_2, rev(uni_ccdf), type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- 1:10
for (k in ks){
  e_k <- reverse_weibull_k_moment(k = k, shape = -shape, scale = scale)
  cota_k <-  e_k / (x_seq_2)^(k)
  lines(x_seq_2, cota_k, col = colors_palette[k])
}
#dev.off()

# Reverse weibull ---------------------------------------------------------


# Make mixture
location <- 40000
scale <- 100
shape <- 1/8

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

colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(100)


pdf("markov_k_rev_weibull_8.pdf", height = 4, width = 6)
plot(x_seq, uni_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF", xlim = c(0, 200))
ks <- 1:200
for (k in ks){
  e_k <- weibull_k_moment(k = k, shape = 1/shape, scale = scale)
  cota_k <-  e_k / (x_seq)^(k)
  lines(x_seq, cota_k, col = colors_palette[k])
}
dev.off()




# Frechet -----------------------------------------------------------------

 
# Make mixture
location <- 40000
scale <- 100
shape <- 1/8


x_seq <- lseq(from = 1, to = 10000, length.out = 500)
#x_seq <- lseq(from = 1, to = qfrechet(1- 1e-1, shape = shape, scale = scale), length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - evd::pfrechet(x, shape = 1/shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(7)

pdf("markov_k_frechet_8.pdf", height = 4, width = 6)
plot(x_seq, uni_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- 1:(1/shape-1)
for (k in ks){
  e_k <- reverse_weibull_k_moment(k = k, scale = scale, shape = 1/shape)
  cota_k <-  e_k / (x_seq)^(k)
  lines(x_seq,  cota_k, col = colors_palette[k])
}
dev.off()



# Gumbel ------------------------------------------------------------------


# Make mixture
location <- 40000
scale <- 100


x_seq <- seq(from = 0, to = 40, length.out = 500)
#x_seq <- lseq(from = 1, to = qfrechet(1- 1e-1, shape = shape, scale = scale), length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - evd::pgumbel(x, scale = 1))
}
plot(x_seq, uni_ccdf, log = "y")

colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(50)

pdf("markov_k_gumbel.pdf", height = 4, width = 6)
plot(x_seq, uni_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- 1:50
for (k in ks){
  e_k <- gamma(1 + k)
  cota_k <-  e_k / exp(k*x_seq)
  lines(exp(x_seq),  cota_k, col = colors_palette[k])
}
dev.off()
