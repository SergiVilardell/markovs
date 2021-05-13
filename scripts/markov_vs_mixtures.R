library(mixtools)
library(tidyverse)

mixture_fits <- readRDS("best_fit_mixture_data.RDS")

weights <- c(0.8, 0.15, 0.05)
mu <- c(10, 20, 50)
sigma <- c(1, 5, 10)



# data
x <- c(rnorm(1000*weights[1], mean = mu[1], sd = sigma[1]),
       rnorm(1000*weights[2], mean = mu[2], sd = sigma[2]),
       rnorm(1000*weights[3], mean = mu[3], sd = sigma[3])
) %>% sort()

# translate to 0
translation <- min(x)
x_transf <- x - translation

# mixture params

# compute the expected value
ev_mixture <- weights[1] * mu[1] + weights[2] * mu[2] + weights[3] * mu[3] - translation
p_a <- c()
markov <- c()
# get a random high quantile
for (i in 1:length(x_transf)) {
  p_a[i] <- 1  -  (weights[1] * pnorm(x_transf[i], mean = mu[1] - translation, sd = sigma[1]) +
                   weights[2] * pnorm(x_transf[i], mean = mu[2] - translation, sd = sigma[2]) +
                   weights[3] * pnorm(x_transf[i], mean = mu[3] - translation, sd = sigma[3])
  )
  markov[i] <- ev_mixture / x_transf[i]
}
markov_data <- data.frame(x = x_transf, y = p_a, m = markov)
colors <- c(
  "Sample data" = "#000000",
  "Markov" = "#e60000"
)
p <- markov_data %>%
  ggplot() +
  geom_point(aes(x = x, y = y, color = "Sample data")) +
  geom_line(aes(x = x, y = m, color = "Markov"), size = 1.2) +
  labs(
    y = "CCDF",
    x = paste0("Sample "),
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
ggsave(p, file = "markov_mixture.pdf", height = 3, width =6)
