library(tidyverse)
library(ercv)
# Number of Peaks ---------------------------------------------------------

params_list <- readRDS("params_list_1.RDS")

# Constants
n_peaks <- 3
n_samples <- 1000
quant <- c(1-10^-6, 1-10^-9, 1-10^-12)
mus <- c(0, 50, 100)
shapes <- c(1, 2, 5)
scales <- c(10, 50, 100)
weights <- c(0.5, 0.3, 0.2)

# Make the mixture
mixture <- list()
for(i in 1:n_peaks){
    mixture[[i]] <- rweib(n_samples, alpha = shapes[i], lambda = scales[i], mu = mus[i] ) 
}
mixture_full <- unlist(mixture)


# Theoretical Expected value of Weibull Mixture
ev_mixture <- weibull_mixture_ev(weights = weights, lambda = scales, alpha = shapes, mu = mus)
# Compute theoretical EV

ev_mixture/q
markov <- c()
p_a <- c()
# get a random high quantile
for (i in 1:length(mixture_full)) {
  p_a[i] <- 1  -  (weights[1] * pweib(mixture_full[i], alpha = shapes[1], lambda = scales[1], mu = mus[1]) +
                   weights[2] * pweib(mixture_full[i], alpha = shapes[2], lambda = scales[2], mu = mus[2]) +
                   weights[3] * pweib(mixture_full[i], alpha = shapes[3], lambda = scales[3], mu = mus[3])
                   )
  markov[i] <- ev_mixture / mixture_full[i]
}

markov_data <- data.frame(x = mixture_full, y = p_a, m = markov)
colors <- c(
  "Sample data" = "#000000",
  "Markov" = "#e60000"
)
pl <- markov_data %>%
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
pl

ggsave(pl, file = "best_markov_bound_params_1.pdf", height = 3, width = 6)
