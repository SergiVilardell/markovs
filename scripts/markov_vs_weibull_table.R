library(tidyverse)
library(ercv)
library(gtools)
source("markov_functions.R")

# Constants
n_peaks <- 3
weights <- c(w1 = 1/3, w2 = 1/3, w3 = 1/3)
lambda <- c(10,10,10, 50,50,50, 100,100,100) 
mu <- c(0,50,100) 
alpha <- c(1,1,1, 2,2,2, 5,5,5) 
quants <- c(q1 = 1-10^-6, q2 = 1-10^-9, q3 = 1-10^-12)

# All permutations of size 3 for the parameters
alpha_perm <- permutations(length(alpha), n_peaks, alpha, set = F) %>% unique()
lambda_perm <- permutations(length(lambda), n_peaks, lambda, set = F) %>% unique()

# Create the list of all input parameters
params_list <- get_params_list(alpha = alpha_perm, lambda = lambda_perm)

# Compute EV of the mixture and Markov's inequality for the selected quantiles for each param instance of the list
results <- data.frame()
for(i in 1:length(params_list)){
  
  q_1 <- qweib(p = quants[1], lambda = p$lambda[3], alpha = p$alpha[3], mu = mu[3])
  q_2 <- qweib(p = quants[2], lambda = p$lambda[3], alpha = p$alpha[3], mu = mu[3])
  q_3 <- qweib(p = quants[3], lambda = p$lambda[3], alpha = p$alpha[3], mu = mu[3])
  if(min(c(q1, q2, q3)) == q3){
    next
  }
  print(i)
  p <- params_list[[i]]
  
  ev_mixture <- weibull_mixture_ev(weights = weights, lambda = p$lambda, alpha = p$alpha, mu = mu)
  
  seed_1 <- qweib(p = quants[1], lambda = p$lambda[3], alpha = p$alpha[3], mu = mu[3])
  seed_2 <- qweib(p = quants[2], lambda = p$lambda[3], alpha = p$alpha[3], mu = mu[3])
  seed_3 <- qweib(p = quants[3], lambda = p$lambda[3], alpha = p$alpha[3], mu = mu[3])
  
  solfit_1 <- stats::nlm(mixture_quantile_opt_1, seed_1)
  solfit_2 <- stats::nlm(mixture_quantile_opt_2, seed_2)
  solfit_3 <- stats::nlm(mixture_quantile_opt_3, seed_3)
  
  a_1 <- solfit_1$estimate
  a_2 <- solfit_2$estimate
  a_3 <- solfit_3$estimate
  a <- c(a_1 = a_1,
         a_2 = a_2,
         a_3 = a_3)
  
  markov_bound <- ev_mixture / a
  markov_prop <- markov_bound / (1 - quants)
  markov_diff <- markov_bound - (1 - quants)
  names(markov_prop) <- c("prop_q_1", "prop_q_2", "prop_q_3")
  names(markov_diff) <- c("diff_q_1", "diff_q_2", "diff_q_3")
  names(markov_bound) <- c("ev /a_1", "ev /a_2", "ev /a_3")
  
  res <- c(weights, unlist(p), a, markov_bound, 1 - quants, markov_prop, markov_diff) %>% round(3)
  results <- rbind(results, res)
}

colnames(results) <- names(res)

# Save
write.csv(results, "results_markov_theo_param_5.csv", row.names = F)
saveRDS(params_list,"params_list_5.RDS")




plot1 <- results %>% 
  ggplot()+
  geom_boxplot(aes(x = as.factor(lambda3), y = diff_q_3))+
  ylab("MB - q3")+
  scale_y_log10()

plot1 <- results %>% 
  ggplot()+
  geom_boxplot(aes(x = as.factor(alpha3), y = diff_q_3))+
  ylab("MB - q3")+
  scale_y_log10()








res1 <- read.csv("results_markov_theo_param_1.csv")
res2 <- read.csv("results_markov_theo_param_2.csv")
res3 <- read.csv("results_markov_theo_param_3.csv")
res4 <- read.csv("results_markov_theo_param_4.csv")

res_full <- rbind(res1, res3, res4)

plot1 <- res_full %>% 
  mutate(peaks = as.factor(rep(c(3,2,1), each = nrow(res1)))) %>% 
  ggplot()+
  geom_boxplot(aes(x = peaks, y = MB...q_3.1))+
  ylab("MB - q3")+
  scale_y_log10()
ggsave(plot1, file = "peaks_vs_mb.pdf", height = 3, width = 6)