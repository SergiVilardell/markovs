library(mixtools)

# First check the usefulness of the inequality for a mixture with more weight on the highest values
mean1 <- 0.5
mean2 <- 5

# generate two guassians in different proportions, 0,2 and 0.8
mix1 <- rnorm(200, mean = mean1, sd = 0.2)
mix2 <- rnorm(800, mean = mean2, sd = 2)

# mix them
mixture <- c(mix1, mix2)
fit_mixture <-  normalmixEM(mixture, k=2)
plot(fit_mixture, which=2)


#translate to 0
translation <- abs(min(mixture))
mixture0 <- mixture + translation

# get a random high quantile
a <- 8
p_a <- 1- (0.2*pnorm(8, mean = mean1, sd = 0.2) + 0.8*pnorm(8, mean = mean2, sd = 2))


# compute the expected value
ev_mixture <- 0.2*mean1 + 0.8*mean2 + translation


# check markovs inequality
# P(x >= a) = 0.05
ev_mixture/(a + translation)
