library(mixtools)
library(tidyverse)

mixture_fits <- readRDS("best_fit_mixture_data.RDS")
x <- sort(read.table(paste("data/TEST",2,".txt",sep=""))[,1])

mix <- mixture_fits[[2]]
weights <- mix$lambda
means <- mix$mu
sigmas <- mix$sigma


#translate to 0
translation <- abs(min(x))

# get a random high quantile
a <- 18000
p_a <- 1- (weights[1]*pnorm(a, mean = means[1], sd = sigmas[1]) + 
           weights[2]*pnorm(a, mean = means[2], sd = sigmas[2]) +
           weights[3]*pnorm(a, mean = means[3], sd = sigmas[3])
           )


# compute the expected value
ev_mixture <- weights[1]*means[1] + weights[2]*means[2] + weights[3]*means[3]


# check markovs inequality
# P(x >= a) = 0.05
(ev_mixture - translation)/(a - translation)
