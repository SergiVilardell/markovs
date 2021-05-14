#plots paper
source("scripts/markov_functions.R")
library(truncnorm)
library(tidyverse)
library(data.table)
library(showtext)
library(ggrepel)
font_add("lmroman", regular = "D:/BSC/Latin-Modern-Roman/lmroman12-regular.otf")
font_add("lmroman", regular = "D:/BSC/Latin-Modern-Roman/lmroman10-regular.otf")
showtext_auto()
theme_set(theme_bw()+
            theme(text = element_text(family = "lmroman"),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "black")) 
)
# Gaussian ----------------------------------------------------------------


# Make mixture
mean1 <- 100
mean2 <- 50
mean3 <- 100
sd1 <- 10
sd2 <- 20
sd3 <- 30



x_seq <- seq(from = 1,  to = 600, length.out = 1000)
mixture_ccdf <- c()
for (x in x_seq){
  mixture_ccdf <- c(mixture_ccdf, 1 - pnorm(x, mean = mean1, sd = sd1))
}


colfunc <- colorRampPalette(c("red", "blue"))
colors_palette <- colfunc(100)


# Gaussian markov ---------------------------------------------------------


plot_data <- tibble(theo =  mixture_ccdf, x = x_seq)
for (k in 1:1){
  e_k <- mean1
  cota_k <-  as.tibble(e_k / x_seq^(k))
  color_k <- rep(colors_palette[k], nrow(cota_k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
  #lines(x_seq, cota_k, col = colors_palette[k], lwd = 2)
}

viridis_colors <- viridis::viridis(4)
colors <- c("Gaussian" ="#00000F", 
            "Markov" = viridis_colors[3])

p_gauss <- plot_data %>% 
  ggplot()+
  geom_line(aes(x = x, y = theo, color = "Gaussian"))+
  geom_line(aes(x = x, y = k1, color = "Markov"))+
  scale_y_log10(limits = c(1e-8,2))+
  scale_color_manual(values =  colors,
                     name= "")+
  theme(legend.position = c(0.6, 0.4))+
  guides(color = guide_legend(reverse = TRUE))+
  labs(x = "Time", y = "CCDF")
p_gauss
ggsave(p, width = 3, height = 3, file = "markov_gaussian.pdf") 


# Weibull Markov ----------------------------------------------------------

# Make mixture
location <- 40000
scale <- 100
shape <- 1/4

# kth moment of the weibull distribution
weibull_k_moment <- function(k, scale, shape){
  (scale^k)*gamma(1 + k/shape)
}

x_seq <- seq(from = 0, to = 600, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pweibull(x, shape = 1/shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

plot_data <- tibble(theo =  uni_ccdf, x = x_seq)
ks <- seq(from = 1, to = 100, by = 5)
for (k in ks){
  e_k <- weibull_k_moment(k = k, shape = 1/shape, scale = scale)
  cota_k <-  as.tibble(e_k / (x_seq)^(k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
}

viridis_colors <- viridis::viridis(4)
colors <- c("Weibull" ="#00000F", 
            "Markov" = viridis_colors[3])

p_weib <- plot_data %>% 
  ggplot()+
  geom_line(aes(x = x, y = theo, color = "Weibull"))+
  geom_line(aes(x = x, y = k1, color = "Markov"))+
  scale_y_log10(limits = c(1e-8,2))+
  scale_color_manual(values =  colors,
                     name= "")+
  theme(legend.position = c(0.6, 0.4))+
  labs(x = "Time", y = "CCDF")
p_weib


# Gumbel Markov -----------------------------------------------------------


location <- 40000
scale <- 10


x_seq <- seq(from = 0, to = 600, length.out = 500)
#x_seq <- lseq(from = 1, to = qfrechet(1- 1e-1, shape = shape, scale = scale), length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pexp(x, rate = 1/scale))
}
plot(x_seq, uni_ccdf, log = "y")

plot_data <- tibble(theo =  uni_ccdf, x = x_seq)
ks <- 1:1
for (k in ks){
  e_k <- scale
  cota_k <-  as.tibble(e_k / x_seq)
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
}

viridis_colors <- viridis::viridis(4)
colors <- c("Gumbel" ="#00000F", 
            "Markov" = viridis_colors[3])


p_gumbel <- melt(plot_data, id.vars = c("x", "theo")) %>% 
  ggplot()+
  geom_line(aes(x = x, y = theo, color = "Gumbel"))+
  geom_line(aes(x = x, y = value, color = "Markov"), show.legend = F)+
  scale_y_log10(limits = c(1e-8,5))+
  scale_color_manual(values =  colors,
                     name= "")+
  theme(legend.position = c(0.6, 0.4))+
  xlim(c(0,600))+
  labs(x = "Time", y = "CCDF")
p_gumbel
ggsave(plo_gumbel, width = 3, height = 3, file = "markov_gumbel.pdf") 



# Combination markov ------------------------------------------------------
library(ggpubr)
ggsave(ggarrange(p_gauss , p_weib,nrow = 1), width = 6, height = 3, file = "markov_comb.pdf")

# Multi-example -----------------------------------------------------------

#Gaussian


# Make mixture
mean1 <- 100
sd1 <- 10


x_seq <- seq(from = 1,  to = 600, length.out = 1000)
mixture_ccdf <- c()
for (x in x_seq){
  mixture_ccdf <- c(mixture_ccdf, 1 - pnorm(x, mean = mean1, sd = sd1))
}

ks <- c(10,50,90)
plot_data <- tibble(theo =  mixture_ccdf, x = x_seq)
for (k in ks){
  e_k <- tail(gaussian_expected_value_k(k + 1, mean1, sd1), n = 1)
  cota_k <-  as.tibble(e_k / x_seq^(k))
  color_k <- rep(colors_palette[k], nrow(cota_k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
  #lines(x_seq, cota_k, col = colors_palette[k], lwd = 2)
}

melted_plot_data <- melt(plot_data, id.vars = "x")

p <- melted_plot_data %>% 
  ggplot()+
  geom_line(aes(x = x, y = value, color = variable, size = variable))+
  scale_y_log10(limits = c(1e-8,10))+
  xlim(c(0,300))+
  scale_color_manual(values = c("black", "red", "darkviolet", "blue"), 
                     name= "",
                     labels = c("Gaussian","k = 10","k = 50","k = 90"))+
  scale_size_manual(values = c(1, 1, 1, 1), guide = FALSE)+
  theme(legend.position = c(0.25, 0.3),
        legend.key.size = unit(0.4, 'cm'))+
  labs(x = "Time", y = "CCDF")
p

#Weibull


# Make mixture
location <- 40000
scale <- 100
shape <- 1/4

# kth moment of the weibull distribution
weibull_k_moment <- function(k, scale, shape){
  (scale^k)*gamma(1 + k/shape)
}

x_seq <- seq(from = 0, to = 600, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pweibull(x, shape = 1/shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

plot_data <- tibble(theo =  uni_ccdf, x = x_seq)
ks <- c(10,50,90)
for (k in ks){
  e_k <- weibull_k_moment(k = k, shape = 1/shape, scale = scale)
  cota_k <-  as.tibble(e_k / (x_seq)^(k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
}

melted_plot_data <- melt(plot_data, id.vars = "x")

p_weib <- melted_plot_data %>% 
  ggplot()+
  geom_line(aes(x = x, y = value, color = variable, size = variable))+
  scale_y_log10(limits = c(1e-8,10))+
  xlim(c(0,300))+
  scale_color_manual(values = c("black", "red", "darkviolet", "blue"), 
                     name= "",
                     labels = c("Weibull","k = 10","k = 50","k = 90"))+
  scale_size_manual(values = c(1, 1, 1, 1), guide = FALSE)+
  theme(legend.position = c(0.3, 0.3),
        legend.key.size = unit(0.4, 'cm'))+
  labs(x = "Time", y = "CCDF")
p_weib

ggsave(p, width = 5, height = 3, file = "markov_gaussian_multik_example.pdf") 

ggsave(ggarrange(p , p_weib), width = 6, height = 3, file = "markov_comb_K3.pdf")



# Contour bound -----------------------------------------------------------

plot(x_seq, mixture_ccdf, type = "l", log = "y", xlab = "Time", ylab = "CCDF")
ks <- seq(from = 1, to = 100, by = 5)
plot_data <- tibble(theo =  mixture_ccdf, x = x_seq)
for (k in ks){
  e_k <- tail(gaussian_expected_value_k(k + 1, mean1, sd1), n = 1)
  cota_k <-  as.tibble(e_k / x_seq^(k))
  color_k <- rep(colors_palette[k], nrow(cota_k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
  #lines(x_seq, cota_k, col = colors_palette[k], lwd = 2)
}



# Get best k --------------------------------------------------------------

p_test <- 1-1e-3

get_data_best_k <- function(p_test){
  x_test <- qnorm(p_test, mean = mean1, sd = sd1)
  
  ks <- seq(from = 5, to = 200, by = 2)
  iter <- 1
  data_test_k <- c()
  for (k in ks){
    e_k <- tail(gaussian_expected_value_k(k + 1, mean1, sd1), n = 1)
    cota_k <-  e_k / x_test^(k)
    names(cota_k) <- paste0("k", k)
    data_test_k[iter] <- cota_k
    iter <- iter + 1
    #lines(x_seq, cota_k, col = colors_palette[k], lwd = 2)
  }
  
  data_test_k <- na.omit(data_test_k)
  data_test_k <- data_test_k[data_test_k!=0]
  data_test_k <- data_test_k[data_test_k!=Inf]
  index_min_k <- which(data_test_k==min(data_test_k))
  return(data_test_k)
}

data_test_k_6 <- get_data_best_k(1-1e-6)
data_test_k_9 <- get_data_best_k(1-1e-9)
data_test_k_12 <- get_data_best_k(1-1e-12)
qnorm(p_test, mean = mean1, sd = sd1)

get_min_k <- function(data_plot){
  index_min_k <- which(data_plot$y==min(data_plot$y))
  data_plot$x[index_min_k]
}
get_min_k(data_plot_best_k_12) 

data_plot_best_k_6 <- tibble(y = data_test_k_6, x = ks[1:length(data_test_k_6)], label = rep(get_min_k(data_plot_best_k_6),
                                                                                             length(data_test_k_6)))
data_plot_best_k_9 <- tibble(y = data_test_k_9, x = ks[1:length(data_test_k_9)],label = rep(get_min_k(data_plot_best_k_9),
                                                                                            length(data_test_k_9)))
data_plot_best_k_12 <- tibble(y = data_test_k_12, x = ks[1:length(data_test_k_12)],label = rep(get_min_k(data_plot_best_k_12),
                                                                                               length(data_test_k_12)))




viridis_colors <- viridis::viridis(4)
colors <- c("1e-6" = viridis_colors[3], 
            "1e-9" = viridis_colors[1],
            "1e-12" = viridis_colors[2])

labels <- factor(c("1e-06","1e-09","1e-12"), levels = c("1e-6","1e-9","1e-12"))
colors <- c("1e-06" = viridis_colors[3], 
            "1e-09" = viridis_colors[2],
            "1e-12" = viridis_colors[1])
pl_best_k <- data_plot_best_k %>% 
  ggplot()+
  geom_point(data= data_plot_best_k_6,aes(x = x, y = y, color = "1e-06"))+
  geom_point(data= data_plot_best_k_9,aes(x = x, y = y,color = "1e-09"))+
  geom_point(data= data_plot_best_k_12,aes(x = x, y = y,color = "1e-12"))+
  geom_label_repel(data          = subset(data_plot_best_k_6, x == get_min_k(data_plot_best_k_6) ),
                   aes(x = x, y = y, color = "1e-06", label = label),
                   show.legend = FALSE,
                   size          = 4,
                   vjust = 1,
                   box.padding   = 0.1,
                   nudge_x = 2,
                   point.padding = 0.5,
                   force         = 100,
                   segment.size  = 0.2,
                   segment.color = "1e-06",
                   direction     = "y")+
  geom_label_repel(data          = subset(data_plot_best_k_9, x == get_min_k(data_plot_best_k_9) ),
                   aes(x = x, y = y, color = "1e-09", label = label),
                   show.legend = FALSE,
                   size          = 4,
                   vjust = 1,
                   nudge_x = 3,
                   box.padding   = 0.1,
                   point.padding = 0.5,
                   force         = 100,
                   segment.size  = 0.2,
                   segment.color = "1e-09",
                   direction     = "y")+
geom_label_repel(data          = subset(data_plot_best_k_12, x == get_min_k(data_plot_best_k_12) ),
                 aes(x = x, y = y, color = "1e-12", label = label),
                 show.legend = FALSE,
                 size          = 4,
                 vjust = 1,
                 nudge_x = 3,
                 box.padding   = 0.1,
                 point.padding = 0.5,
                 force         = 100,
                 segment.size  = 0.2,
                 segment.color = "1e-12",
                 direction     = "y")+
  scale_color_manual(values = colors)+
  labs(
    y = "",
    color = "Quantile"
  ) +
  scale_y_log10()+
  theme(legend.position = c(0.1, 0.72))+
  labs(x = "k", y = "CCDF")
pl_best_k

ggsave(pl_best_k, height = 3, width = 5, file = "plot_best_k.pdf")

countour <- plot_data[, -c(1:2)] %>% 
  t() %>% 
  apply(2,min)

countour_data <- cbind(plot_data[, c(1:2)], countour)




colors <- c("Contour Bound" = "#00cfc4", "Distribution" = "#000000")
pl <- countour_data %>% 
  ggplot()+
  geom_line(aes(x = x, y = countour, color = "Contour Bound"), size = 1.5)+
  geom_line(aes(x = x, y = theo, color = "Distribution"))+
  scale_y_log10(limits = c(1e-8,10))+
  scale_color_manual(values = colors)+
  labs(
    y = "",
    color = ""
  ) +
  guides(colour = guide_legend(override.aes = list(
    size = c(1.5, 0.5)
  ))) +
  theme(legend.position = c(0.2, 0.4))+
  labs(x = "Time", y = "CCDF")
  
ggsave(pl, width = 5, height = 3, file = "markov_gaussian_contour.pdf") 




# All bounds --------------------------------------------------------------

plo <- melt(plot_data, id.vars = c("x", "theo")) %>% 
  ggplot()+
  geom_line(aes(x = x, y = theo))+
  geom_line(aes(x = x, y = value, color = variable), show.legend = F)+
  scale_y_log10(limits = c(1e-8,10))+
  scale_color_viridis_d(end = 0.8)+
  labs(x = "Time", y = "CCDF")
ggsave(plo, width = 5, height = 3, file = "markov_gaussian_all_bounds.pdf")   



# Gumbel ------------------------------------------------------------------

# Make mixture
location <- 40000
scale <- 1


x_seq <- seq(from = 0, to = 60, length.out = 500)
#x_seq <- lseq(from = 1, to = qfrechet(1- 1e-1, shape = shape, scale = scale), length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - evd::pgumbel(x, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

plot_data <- tibble(theo =  uni_ccdf, x = x_seq)
ks <- 1:50
for (k in ks){
  e_k <- gamma(1 + k)
  cota_k <-  as.tibble(e_k / exp(k*x_seq))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
}

plo_gumbel <- melt(plot_data, id.vars = c("x", "theo")) %>% 
  ggplot()+
  geom_line(aes(x = x, y = theo))+
  geom_line(aes(x = exp(x), y = value, color = variable), show.legend = F)+
  scale_y_log10(limits = c(1e-15,10))+
  xlim(c(0,40))+
  scale_color_viridis_d(end = 0.8)+
  labs(x = "Time", y = "CCDF")
plo_gumbel
ggsave(plo_gumbel, width = 5, height = 3, file = "markov_gumbel.pdf")



# Frechet -----------------------------------------------------------------

# Make mixture
location <- 40000
scale <- 100
shape <- 1/2


x_seq <- lseq(from = 1, to = 1e8, length.out = 500)
#x_seq <- lseq(from = 1, to = qfrechet(1- 1e-1, shape = shape, scale = scale), length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - evd::pfrechet(x, shape = 1/shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")

plot_data <- tibble(theo =  uni_ccdf, x = x_seq)
ks <- 1:(1/shape-1)
for (k in ks){
  e_k <- reverse_weibull_k_moment(k = k, scale = scale, shape = 1/shape)
  cota_k <-  as.tibble(e_k / (x_seq)^(k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
}

plo_frech_2 <- melt(plot_data, id.vars = c("x", "theo")) %>% 
  ggplot()+
  geom_line(aes(x = x, y = theo))+
  geom_line(aes(x = x, y = value, color = variable), show.legend = F)+
  scale_y_log10(limits = c(1e-13,10))+
  scale_color_viridis_d(begin = 0.8, end = 1)+
  labs(x = "Time", y = "CCDF")
plo_frech_2
ggsave(plo_frech_2, width = 5, height = 3, file = "markov_frechet_2.pdf")



# Reverse Weibull ---------------------------------------------------------

# Make mixture
location <- 40000
scale <- 100
shape <- 1/2

# kth moment of the weibull distribution
weibull_k_moment <- function(k, scale, shape){
  (scale^k)*gamma(1 + k/shape)
}

x_seq <- seq(from = 0, to = 600, length.out = 500)
uni_ccdf <- c()
for (x in x_seq){
  uni_ccdf <- c(uni_ccdf, 1 - pweibull(x, shape = 1/shape, scale = scale))
}
plot(x_seq, uni_ccdf, log = "y")



plot_data <- tibble(theo =  uni_ccdf, x = x_seq)
ks <- seq(from = 1, to = 100, by = 5)
for (k in ks){
  e_k <- weibull_k_moment(k = k, shape = 1/shape, scale = scale)
  cota_k <-  as.tibble(e_k / (x_seq)^(k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)
}

plo_weib_2 <- melt(plot_data, id.vars = c("x", "theo")) %>% 
  ggplot()+
  geom_line(aes(x = x, y = theo))+
  geom_line(aes(x = x, y = value, color = variable), show.legend = F)+
  scale_y_log10(limits = c(1e-13,10))+
  scale_color_viridis_d(end = 0.8)+
  labs(x = "Time", y = "CCDF")
plo_weib_2
ggsave(plo_weib_2, width = 5, height = 3, file = "markov_weibull_2.pdf")


# Markov vs EVT -----------------------------------------------------------

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
plot_data <- tibble(theo =  mixture_ccdf, x = x_seq)

# Threshold from where to fit the exponential
thresholds <- seq(from=50, to = 160, by = 2)

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
  ccdf_mostra_trunc <- as.tibble(lambda * (1 - pexp(x_theo - trunc, rate_mostra_trunc)))
  x_theo <- as.tibble(x_theo)
  colnames(ccdf_mostra_trunc) <- paste0("threshold", thresholds[i])
  colnames(x_theo) <- paste0("x_theo_", thresholds[i])
  plot_data<- cbind(plot_data, x_theo, ccdf_mostra_trunc)
}

a <- plot_data[colnames(plot_data[grep("x_theo_", colnames(plot_data))])] %>% melt()
colnames(a) <- c("x_theo", "value_x")
b <- plot_data[colnames(plot_data[grep("threshold", colnames(plot_data))])] %>% melt()
colnames(b) <- c("threshold", "value_t")

d <- cbind(a,b)
plo_evt<- d %>% 
  ggplot()+
  geom_line(aes(x = value_x, y = value_t, color = threshold), show.legend = F)+
  geom_line(data = plot_data, aes(x = x, y = theo), show.legend = F, size = 1.5)+
  scale_y_log10(limits = c(1e-15,10))+
  scale_color_viridis_d(end = 0.8)+
  labs(x = "Time", y = "CCDF")
plo_evt 

ggsave(plo_evt, width = 5, height = 3, file = "markov_evt.pdf")




# Markov evt gaussian -----------------------------------------------------

source("scripts/markov_functions.R")
library(CharFun)
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
      0.399 * dnorm(x, mean = mean2, sd = sd2) +
      0.001 * dnorm(x, mean = mean3, sd = sd3)
  )
  mixture_ccdf <- c(mixture_ccdf, 1 - (
    0.6 * pnorm(x, mean = mean1, sd = sd1) +
      0.399 * pnorm(x, mean = mean2, sd = sd2) +
      0.001 * pnorm(x, mean = mean3, sd = sd3)
  ))
}

plot_data <- tibble(theo =  mixture_ccdf, x = x_seq)
ks <- seq(from = 5, to = 120, by = 3)
for (k in ks){
  e_k <- 0.6 * (tail(gaussian_expected_value_k(k+1, mean1, sd1), n = 1) )/(1-pnorm(0, mean = mean1, sd = sd1)) + 
    0.399 * (tail(gaussian_expected_value_k(k+1,mean2, sd2), n = 1))/(1-pnorm(0, mean = mean2, sd = sd2)) + 
    0.001 * (tail(gaussian_expected_value_k(k+1, mean3, sd3), n = 1))/(1-pnorm(0, mean = mean3, sd = sd3))
  cota_k <-  as.tibble(e_k / (x_seq)^(k))
  colnames(cota_k) <- paste0("k", k)
  plot_data <- cbind(plot_data, cota_k)

}

plo_evt_gauss <- melt(plot_data, id.vars = c("x", "theo")) %>% 
  ggplot()+
  geom_line(aes(x = x, y = theo))+
  geom_line(aes(x = x, y = value, color = variable), show.legend = F)+
  scale_y_log10(limits = c(1e-15,10))+
  scale_color_viridis_d(end = 0.8)+
  labs(x = "Time", y = "CCDF")
plo_evt_gauss
ggsave(plo_evt_gauss, width = 5, height = 3, file = "markov_evt_gauss.pdf")

ggsave(ggarrange(plo_evt , plo_evt_gauss), width = 10, height = 3, file = "markov_evt_gauss_comb.pdf")
