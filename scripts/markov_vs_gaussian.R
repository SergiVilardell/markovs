library(mixtools)
library(tidyverse)

# Here we will compare how tight MI when changing Gaussian parameters

sd_seq <- seq(0.001, 1000, length.out = 1000)
results <- data.frame()

for(j in 1:length(sd_seq)){
  m <- 5
  sd <- sd_seq[j]
  gaussian_narrow <- sort(rnorm(3000, mean = m, sd = sd))
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
  
  res <- c(min(markov), 1-max(p_a), sd)
  results <- rbind(results, res)
}


names(results) <- c("markov", "p", "sd")
plot(y = results$markov -results$p, x = results$sd, xlab = "Standard Deviation", ylab = "Minimum Markov")

colors <- c(
  "Sample data" = "#000000",
  "Markov" = "#e60000"
)
data.frame(x = gaussian_narrow_transf, y = p_a, m = markov) %>%
  ggplot() +
  geom_point(aes(x = x, y = y, color = "Sample data")) +
  geom_line(aes(x = x, y = m, color = "Markov"), size = 1.2) +
  labs(
    y = "CCDF",
    x = paste0("Sample ", i),
    color = ""
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.70, 0.82),
    legend.background = element_rect(fill = alpha("white", 0.)),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  ) +
  guides(colour = guide_legend(override.aes = list(
    shape = c(95, 95),
    linetype = c(1, 1),
    size = c(1, 1)
  ))) +
  ylim(0, 1) +
  scale_color_manual(values = colors) +
  guides(fill = guide_legend(override.aes = list(
    shape = c(0, 0),
    linetype = c(1, 1),
    size = c(1, 1)
  )))

summary(lm(data = results, markov ~ sd))