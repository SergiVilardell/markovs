
source("scripts/markov_functions.R")

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
      0.3 * dnorm(x, mean = mean2, sd = sd2) +
      0.1 * dnorm(x, mean = mean3, sd = sd3)
  )
  mixture_ccdf <- c(mixture_ccdf, 1 - (
    0.6 * pnorm(x, mean = mean1, sd = sd1) +
      0.3 * pnorm(x, mean = mean2, sd = sd2) +
      0.1 * pnorm(x, mean = mean3, sd = sd3)
  ))
}


colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(100)

pdf("chernoff_s_1e-4_1.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
s_seq <- seq(0.0001, 1, length.out = 200)
colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(length(s_seq))
gaussian_expected_value_k(k = 4, mu = mean1, sigma = sd1)

iter <- 1
for (s in s_seq){
  chernoff <- (0.6 * mgf_gaussian(t = s, mu = mean1, sigma = sd1) + 
                              0.3 * mgf_gaussian(t = s, mu = mean2, sigma = sd2) + 
                              0.1 * mgf_gaussian(t = s, mu = mean3, sigma = sd3) )/  exp(s * x_seq)
  lines(x_seq, chernoff, col = colors_palette[iter])
  iter <- iter + 1
}
dev.off()

lines(xs, ys, col = "red")
