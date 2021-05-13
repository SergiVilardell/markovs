library(mixtools)
library(tidyverse)


mixture_fits <- readRDS("best_fit_mixture_data.RDS")
get_markov_inequality_data <- function(fits) {
  markov_data <- list()
  for (j in 1:length(fits)) {
    # data
    x <- sort(read.table(paste0("data/TEST", j - 1, ".txt"))[, 1])
    # translate to 0
    translation <- min(x)
    x_transf <- x - translation
    
    # mixture params
    mix <- fits[[j]]
    weights <- mix$lambda
    means <- mix$mu
    sigmas <- mix$sigma
    
    # compute the expected value
    ev_mixture <- weights[1] * means[1] + weights[2] * means[2] + weights[3] * means[3]
    p_a <- c()
    markov <- c()
    # get a random high quantile
    for (i in 1:length(x_transf)) {
      p_a[i] <- 1 - (weights[1] * pnorm(x_transf[i], mean = means[1], sd = sigmas[1]) +
                       weights[2] * pnorm(x_transf[i], mean = means[2], sd = sigmas[2]) +
                       weights[3] * pnorm(x_transf[i], mean = means[3], sd = sigmas[3])
      )
      markov[i] <- ev_mixture / x_transf[i]
    }
    markov_data[[j]] <- data.frame(x = x_transf, y = p_a, m = markov)
  }
  return(markov_data)
}


markov_data <- get_markov_inequality_data(mixture_fits)

pdf(file = "markovs_train_data.pdf", width = 6, height = 3)
colors <- c(
  "Sample data" = "#000000",
  "Markov" = "#e60000"
)
for (i in 1:length(markov_data)) {
  p <- markov_data[[i]] %>%
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
  plot(p)
}
dev.off()