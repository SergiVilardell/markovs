library(truncnorm)
# Make mixture
mean1 <- 5
mean2 <- 50
mean3 <- 100
sd1 <- 10
sd2 <- 10
sd3 <- 10

# Range of mixture
x_seq <- seq(1, 170, length.out = 1000)

# COmpute mixture ccdf
mixture_ccdf <- c()
for (x in x_seq) {
  mixture_ccdf <- c(mixture_ccdf, 1 - (
    0.6 * pnorm(x, mean = mean1, sd = sd1) +
      0.399 * pnorm(x, mean = mean2, sd = sd2) +
      0.001 * pnorm(x, mean = mean3, sd = sd3)
  ))
}

# Get a big sample to compute the expected value^k
mostra <- sort(c(rnorm(6000000, mean = mean1, sd = sd2), rnorm(3990000, mean = mean2, sd = sd2), rnorm(10000, mean = mean3, sd = sd3)))

# Initiate pdf
pdf("markov_vs_evt_001.pdf", height = 4, width = 6)

# Threshold from where to fit the exponential
thresholds <- 50:240

# Colors for the plot
colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(length(thresholds))

# Make plot empty to draw the line later
plot(x_seq, mixture_ccdf, log = "y", xlab = "Time", ylab = "CCDF", type = "n")

# compute Markov's
for (i in 1:length(thresholds)) {
  max_dist <- max(x_seq)
  trunc <- thresholds[i]
  
  # Get lambda to scale the probability to the threshold
  lambda <- 1 - (0.6 * pnorm(trunc, mean = mean1, sd = sd1) +
                 0.399 * pnorm(trunc, mean = mean2, sd = sd2) +
                 0.001 * pnorm(trunc, mean = mean3, sd = sd3))
  x_theo <- seq(trunc, max_dist, length.out = 1000)
  
  # For high quantiles, get the truncated normal mean
  if (trunc > 130) {
    rate_mostra_trunc <- 1 / (etruncnorm(a = trunc, b = Inf, mean = mean3, sd = sd3) - trunc)
  } else {
    mostra_trunc <- mostra[mostra > trunc] - trunc
    # The rate for the exp is 1/mean
    rate_mostra_trunc <- 1 / mean(mostra_trunc)
  }
  
  # Compute the ccdf from the threshold
  ccdf_mostra_trunc <- lambda * (1 - pexp(x_theo - trunc, rate_mostra_trunc))
  lines(x_theo, ccdf_mostra_trunc, col = colors_palette[i])
}

# Draw the distribution line now so it is on top of the others
points(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF", lwd = 4)
dev.off()

pdf("markov_k_100.pdf", height = 6, width = 6)
plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- 1:100
for (k in ks) {
  mostra <- c(rnorm(6000000, mean = mean1, sd = sd2), rnorm(3990000, mean = mean2, sd = sd2), rnorm(10000, mean = mean3, sd = sd3))
  cotak <- mean(mostra^k) / x_seq^k
  lines(x_seq, cotak, col = colors_palette[k])
}
dev.off()








# Silva2017 distros -------------------------------------------------------

library(extRemes)
loc <- 40000
scale <- 100
shapes <- c(-1 / 2, -1 / 4, -1 / 8, 0, 1 / 8, 1 / 4, 1 / 2)

mean_GEV <- function(loc, scale, shape) {
  loc + gamma(1 - shape) * (scale / shape)
}

mean_GEV(loc = loc, scale = scale, shape = shapes[2])

x_seq <- seq(39800, 40450, length.out = 100)
gev_ccdf <- 1 - extRemes::pevd(q = x_seq, loc = loc, scale = scale, shape = shapes[2], type = "GEV")
mean(extRemes::revd(n = 1000000, loc = loc, scale = scale, shape = shapes[1], type = "GEV"))

plot(x_seq, gev_ccdf, log = "y")

pdf("markov_gev_100_2.pdf", height = 6, width = 6)
plot(x_seq, gev_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(100)

ks <- seq(2, 100, length.out = 100)
cotak <- c()

for (j in 1:length(ks)) {
  # cotak <- mean(mostra^k) / x_seq^k
  cotak <- mean_GEV(loc = loc, scale = scale, shape = shapes[2] / ks[j]) * ks[j] / x_seq^ks[j]
  lines(x_seq, cotak, col = colors_palette[j])
}
dev.off()

lines(xs, ys, col = "red")
