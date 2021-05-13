# Create parameters list for all thepermutations
get_params_list <- function(alpha = NULL, lambda = NULL){
  params_list <- list()
  iter <- 1
  for(i in 1:nrow(alpha)){
    for(k in 1:nrow(lambda)){
      alphas <- alpha[i,]
      lambdas <- lambda[k,]
      params_list[[iter]] <- list(alpha = alphas,
                                  lambda = lambdas)
      iter <- iter + 1
    }
  }
  return(params_list)
}

# Get the cumulative prob of a given quantile of the weibull
pweib <- function(q = NULL, lambda = NULL, alpha = NULL, mu = NULL) {
  p <- 1 - exp(- ((q - mu)/lambda)^(alpha))
  return(p)
}

# Get the quantile of a given cumulative prob of the weibull
qweib <- function(p = NULL, lambda = NULL, alpha = NULL, mu = NULL) {
  quant <- mu + lambda * (log(1 / (1-p)))^(1 / alpha)
  return(quant)
}

# Get random samples of a weibull
rweib <- function(n = NULL, lambda = NULL, alpha = NULL, mu = NULL) {
  probs <- runif(n)
  quant <- qweib(p = probs, lambda = lambda, alpha = alpha, mu = mu)
  return(quant)
}

# Get the expected value of the weibull
weibull_ev <- function(lambda = NULL, alpha = NULL, mu = NULL){
  ev <- lambda * gamma(1 + 1/alpha) + mu
  return(ev)
}

# Get the expected value of the mixture
weibull_mixture_ev <- function(weights = NULL, lambda = NULL, alpha = NULL, mu = NULL){
  mixture_ev <- sum(weights * weibull_ev(lambda = lambda, alpha = alpha, mu = mu))
  return(mixture_ev)
}

# Optimization function for quantile 0.9
mixture_quantile_opt_1 <- function(q) {
  (0.9 -   sum(weights * (1 - exp(- ((q - mu)/lambda)^(alpha)))))^2
}

# Optimization function for quantile 0.99
mixture_quantile_opt_2 <- function(q) {
  (0.99 -   sum(weights * (1 - exp(- ((q - mu)/lambda)^(alpha)))))^2
}

# Optimization function for quantile 0.999
mixture_quantile_opt_3 <- function(q) {
  (0.999 -   sum(weights * (1 - exp(- ((q - mu)/lambda)^(alpha)))))^2
}

# Expected value of Gaussian distribution to the power of K
gaussian_expected_value_k <- function(k, mu, sigma){
  mu_r <- c()
  mu_r[1] <- 1
  mu_r[2] <- mu
  if(k >= 2){
    mu_r[3] <- mu^2 + sigma^2
  }
  if (k >= 3) {
    for (i in 3:(k-1)) {
      mu_r[i+1] <- (i-1) * sigma^2 * mu_r[i - 1] + mu * mu_r[i]
    }
  }
  if(k == 1){
    return(mu_r[2])
  } else {
    return(mu_r[-1])
  }
}

# Expected value of Gaussian distribution to the power of K
gaussian_trunc_expected_value_k <- function(k, sigma){
  if((k%%2) == 0){
    mom <- (factorial(2*k)*sigma^(2*k))/(factorial(k)*2^k)
  } else { mom <- 0}
  return(mom)
}

# Expected absolute value of Gaussian distribution to the power of K
gaussian_abs_expected_value_k <- function(k, mu, sigma){
 (sigma^k)*(2^(k/2))*(gamma((1+k)/2)/sqrt(pi))*Re(hypergeom1F1(-k/2, 1/2, -0.5*(mu/sigma)^2))
}

# Expected value of Absolute Gaussian distribution to the power of K

gaussian_absolute_expected_value_k <- function(k, sigma){
 pi^((k-1)/2)*sigma^(-k)*gamma(0.5*(k+1))
}

# Expected value of moment generating function to the power of K
chernoff_expected_value_k <- function(k, t, mu, sigma){
  e_k <- gaussian_expected_value_k(k = k + 1, mu = mu, sigma = sigma)
  taylor <- t^(1:k) / factorial(1:k)
  chernoff_k <- sum(e_k * taylor + 1)
  return(chernoff_k) 
}


mgf_gaussian <- function(t, mu, sigma){
  exp(t * mu + 0.5 * sigma^2 * t^2)
}

# kth moment of the weibull distribution
weibull_k_moment <- function(k, scale, shape){
  (scale^k)*gamma(1 + k/shape)
}

# kth moment of the weibull distribution
reverse_weibull_k_moment <- function(k, scale, shape){
  (scale^k)*gamma(1 - k/shape)
}

# kth moment of the weibull distribution
frechet_k_moment <- function(k, shape){
  gamma(1 - k*shape)
}

gumbel_mgf <- function(t, scale){
  gamma(1 - t*scale)
}

lseq <- function(from, to, length.out) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}
