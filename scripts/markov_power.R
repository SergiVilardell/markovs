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

cota <- (0.6 * 5 + 0.399 * 50 + 0.001 * 100) / x_seq


# pdf("markov_mixture.pdf", height = 6, width = 12)
# par(mfrow=c(1,2))
# plot(x_seq, mixture_density, xlab = "Time", ylab = "Prob Density")
# plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
# lines(x_seq, cota, col = "red")
#  dev.off()


plot(x_seq, mixture_ccdf, log = "y", xlab = "Time", ylab = "CCDF", type = "n")
ks <- 1:100
for (k in ks){
  mostra <- c(rnorm(6000000, mean = mean1, sd = sd1), rnorm(3000000, mean = mean2, sd = sd2), rnorm(1000000, mean = mean3, sd = sd3))
  e_k <- mean(mostra^k)
  cota_k <-  e_k / x_seq^(k)
  lines(x_seq, cota_k, col = colors_palette[k])
}
points(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF", lwd = 4)



colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(200)

pdf("markov_k_vs_evt.pdf", height = 4, width = 10)
par(mfrow=c(1,2))
plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- 1:200
for (k in ks){
  e_k <- 0.6 * (tail(gaussian_expected_value_k(k+1, mean1, sd1), n = 1) )/(1-pnorm(0, mean = mean1, sd = sd1)) + 
                 0.399 * (tail(gaussian_expected_value_k(k+1,mean2, sd2), n = 1))/(1-pnorm(0, mean = mean2, sd = sd2)) + 
                 0.001 * (tail(gaussian_expected_value_k(k+1, mean3, sd3), n = 1))/(1-pnorm(0, mean = mean3, sd = sd3))
  cota_k <-  e_k / x_seq^(k)
  lines(x_seq, cota_k, col = colors_palette[k])
}
dev.off()

plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- 1:100
for (k in ks){
  e_k <- 0.6 * (tail(gaussian_expected_value_k(k+1, mean1, sd1), n = 1)) + 
    0.3 * (tail(gaussian_expected_value_k(k+1, mean2, sd2), n = 1))+ 
    0.1 * (tail(gaussian_expected_value_k(k+1, mean3, sd3), n = 1))
  cota_k <-  e_k / x_seq^(k)
  lines(x_seq, cota_k, col = colors_palette[k])
}


dev.off()
lines(xs, ys, col = "red")

e0_mixture <- (0.6*ptruncnorm(0, mean = mean1, sd = sd1)+
  0.3*pnorm(0, mean = mean2, sd = sd2)+
  0.1*pnorm(0, mean = mean3, sd = sd3))


pdf("markov_trunc_k_200.pdf", height = 6, width = 6)
library(hpa)
plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- 1:100
for (k in ks){
  e_k <- 1/(1-e0_mixture)* (0.6 * (1-pnorm(0, mean = mean1, sd = sd1))* truncatedNormalMoment(k = k, x_lower = 0, x_upper = Inf, mean = mean1, sd = sd1) + 
    0.3 * (1-pnorm(0, mean = mean2, sd = sd2))*truncatedNormalMoment(k = k, x_lower = 0, x_upper = Inf, mean = mean2, sd = sd2)  + 
    0.1 * (1-pnorm(0, mean = mean3, sd = sd3))*truncatedNormalMoment(k = k, x_lower = 0, x_upper = Inf, mean = mean3, sd = sd3) 
  )
  cota_k <-  e_k / (x_seq)^(k)
  lines(x_seq, cota_k, col = colors_palette[k])
}
dev.off()
lines(xs, ys, col = "red")


pdf("markov_abs_k_200.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- 1:200
for (k in ks){
  e_k <- 0.6 * gaussian_abs_expected_value_k(k = k,  mu = mean1, sigma = sd1) + 
    0.3 * gaussian_abs_expected_value_k(k = k,  mu = mean2, sigma = sd2) +
    0.1 * gaussian_abs_expected_value_k(k = k,  mu = mean3, sigma = sd3)
  cota_k <-  e_k / (x_seq)^(k)
  lines(x_seq, cota_k, col = colors_palette[k])
}



dev.off()
lines(xs, ys, col = "red")

# Scale does not matter ---------------------------------------------------


# Make mixture
mean1 <- 5
mean2 <- 50
mean3 <- 100
sd1 <- 20
sd2 <- 20
sd3 <- 20

x_seq <- seq(1, 240, 1)
mixture_density <- c()
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_density <- c(
    mixture_density,
    0.6 * dnorm(x, mean = mean1, sd = sd1) +
      0.3 * dnorm(x, mean = mean2, sd = sd2) +
      0.1 * dnorm(x, mean = mean3, sd = sd3)
  )
  mixture_ccdf <- c(mixture_ccdf, 1 - (
    0.6 * pnorm(x, mean = mean1, sd = sd1) +
      0.3 * pnorm(x, mean = mean2, sd = sd2) +
      0.1 * pnorm(x, mean = mean3, sd = sd3)
  ))
}
#plot(x_seq, mixture_density)

cota <- (0.6 * 5 + 0.3 * 50 + 0.1 * 100) / x_seq

pdf("markov_simple_scale_20.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
lines(x_seq, cota, col = "red")
dev.off()


# Make mixture
mean1 <- 5
mean2 <- 50
mean3 <- 100
sd1 <- 30
sd2 <- 30
sd3 <- 30

x_seq <- seq(1, 310, 1)
mixture_density <- c()
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_density <- c(
    mixture_density,
    0.6 * dnorm(x, mean = mean1, sd = sd1) +
      0.3 * dnorm(x, mean = mean2, sd = sd2) +
      0.1 * dnorm(x, mean = mean3, sd = sd3)
  )
  mixture_ccdf <- c(mixture_ccdf, 1 - (
    0.6 * pnorm(x, mean = mean1, sd = sd1) +
      0.3 * pnorm(x, mean = mean2, sd = sd2) +
      0.1 * pnorm(x, mean = mean3, sd = sd3)
  ))
}
#plot(x_seq, mixture_density)

cota <- (0.6 * 5 + 0.3 * 50 + 0.1 * 100) / x_seq

pdf("markov_simple_scale_30.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
lines(x_seq, cota, col = "red")
dev.off()



########################## 3

n <- 1000
level <- 0.95
mostra <- c(rnorm(0.6 * n, mean = mean1, sd = sd1), rnorm(0.3 * n, mean = mean2, sd = sd2), rnorm(0.1 * n, mean = mean3, sd = sd3))
ecdf_mostra <- ecdf(mostra)

cotak0s <- c()
cotak1s <- c()
cotak2s <- c()
colfunc <- colorRampPalette(c("red", "orange"))
colors_palette <- colfunc(5)

pdf("markov_k_sim_80.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, type = "l", log = "y", col = "grey", xlab = "Time", ylab = "CCDF")
lines(y= 1-sort(ecdf_mostra(mostra)), x=sort(mostra))
for (k in 1:80){
  # cotak<-mean(mostra^k)/xs^k
  x0 <- c()
  x1 <- c()
  x2 <- c()
  for (x in x_seq){
    mostra0 <- mostra / x
    x0 <- c(x0, mean(mostra0^k))
    x1 <- c(x1, t.test(mostra0^k, conf.level = level)$conf.int[1])
    x2 <- c(x2, t.test(mostra0^k, conf.level = level)$conf.int[2])
  }
  x1[x1 < 0] <- NA
  x1[x1 > 1] <- NA
  x2[x2 < 0] <- NA
  x2[x2 > 1] <- NA
  cotak0s <- cbind(cotak0s, x0)
  cotak1s <- cbind(cotak1s, x1)
  cotak2s <- cbind(cotak2s, x2)
}

#lines(x_seq, apply(cotak2s, 1, min, na.rm = T), col = "red")
#lines(x_seq, apply(cotak1s, 1, min, na.rm = T), col = "orange")
lines(x_seq, apply(cotak0s, 1, min, na.rm = T), col = colors_palette[1])
#lines(x_seq, sort(mostra), col = "orange")
abline(h = 1 / n, lty = 2)
dev.off()
