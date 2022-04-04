library(tidyverse)

data_captures <- read_csv("captures_copy.csv")

# rOIZTP <- function(n, p, lambda) { #sample from one-inflated zero-trunkated poisson
#   sample <- rpois(n, lambda*(1 - p)) #sample from poisson
#   sample = sample + 1 #adjust for zero trunkated
#   return(sample)
# }
# 
# rOIZTP_w_ghosts <- function(n, lambda, p) {
#   sample_n_ghosts <- rpois(n, lambda) %>% #sample observations of each bear
#     map_df(~data.frame(n_obs = .x, ghosts = rbinom(1, .x, p))) %>% #calculate n of ghosts per observed individual
#     filter(n_obs != 0) %>% #remove unobserved bears
#     mutate(correct_obs = n_obs - ghosts) #add column with amount of correct observations
#   return(sample_n_ghosts)
# }
# 
# rOIZTNB_w_ghosts <- function(n, size, prob, p) { #sample from one-inflated zero-trunkated negative binomial
#   sample_n_ghosts <- rnbinom(n, size, prob) %>% #sample observations of each bear
#     map_df(~data.frame(n_obs = .x, ghosts = rbinom(1, .x, p))) %>% #calculate n of ghosts per observed individual
#     filter(n_obs != 0) %>% #remove unobserved bears
#     mutate(correct_obs = n_obs - ghosts) #add column with amount of correct observations
#   return(sample_n_ghosts)
# }
# 
# data_t <- data_captures %>% 
#   mutate(simulated = FALSE) %>%
#   group_by(id, simulated) %>%
#   summarise(n_obs = n()) %>%
#   full_join(rOIZTNB_count_ghosts(4000, 6, 0.6, 0.001) %>% mutate(id = str_c("id", as.character(row_number())), simulated = TRUE),
#             by = c("id", "simulated", "n_obs")) 
# 
# data_t %>%
#   ggplot(aes(x = n_obs)) +
#   geom_histogram(binwidth = 1) +
#   facet_wrap(~simulated)
# 
# lloipp_nocov_TEST <- function(la, theta, data){
#   n0 <- length(data)
#   n1 <- sum(data==1)
#   w = 1/(1+exp(-theta))
#   n1*log(w + (1-w)*(la/(exp(la)-1))) + (n0-n1)*log(1-w) + (sum(data)-n1)*log(la) - 
#     (n0-n1)*log(exp(la)-1) - sum(log(factorial(data)))
# }
# 
# data_antal_osb <- data_captures %>% 
#   group_by(id) %>%
#   summarise(n_obs = n())
# 
# opt_lloipp_nocov <- function(param){
#   la = param[1]
#   theta = param[2]
#   n0 <- length(y)
#   n1 <- sum(y==1)
#   w = 1/(1+exp(-theta))
#   -(n1*log(w + (1-w)*(la/(exp(la)-1))) + (n0-n1)*log(1-w) + (sum(y)-n1)*log(la) - 
#     (n0-n1)*log(exp(la)-1) - sum(log(factorial(y))))
# }
# 
# opt_lloiztnb_nocov <- function(param) {
#   li <- param[1]
#   a  <- param[2]
#   
#   # Use the logistic link (Section 2.1 of the main paper)
#   phi <- param[3]
#   w <- 1 / (1 + exp(phi))
#   
#   ymax <- max(y)
#   terms <- weights <- rep(0, ymax)
#   for(ii in 1:ymax) {
#     terms[ii] <- log(a + ii - 1)
#   }
#   weights[1] <- sum(y > 1)
#   for(ii in 2:ymax) {
#     weights[ii] <- sum(y > (ii - 1))
#   }
#   gterm <- terms * weights
#   return(-(sum(log(1 - w) + (y == 1) * (log(w / (1 - w) + a * ((a / (a + li)) ^ a) 
#                                           * (li / (a + li - a * (1 + (li / a)) ^ (1 - a))))) + (1 - (y == 1)) * 
#                (a * log(a) - log(factorial(y)) + y * log(li) - (a + y) * log(a + li) - 
#                   log(1 - (a / (a + li)) ^ a))) + sum(gterm)))
# }

# EGNA FÖRDELNINGAR

rpp <- function(n, lambda) {
  y <- rpois(n, lambda)
  y <- y[y!=0]
  return(y)
}

roipp <- function(n, lambda, p) {
  y <- rpp(n, lambda) %>% # generate sample without inflation
    map_df(~data.frame(y = .x, 
                       real = rbinom(1, .x, (1 - p)))) %>%
    mutate(ghosts = y - real)
  ghosts = sum(y["ghosts"])
  y <- y %>% 
    filter(real != 0) %>%
    select(real) %>%
    data.matrix %>%
    array()
  y = c(y, rep.int(1, ghosts))
  return(y)
}

# log likelihood computing functions

lloipp <- function(param, data = y) {
  lambda <- param[1]
  p <- param[2]
  theta <- lambda*(1 - p)
  oemga <- 1/(1 + lambda*p)
  n1 <- sum(data==1)
  nm <- length(data) - n1
  sum_data_m <- sum(data) - n1
  l <- n1*log((1 - omega) + omega*(theta/(exp(theta) - 1))) + #log likelihood of "1" observations
    nm*log(omega) + sum_data_m*log(theta) - nm*log(exp(theta) - 1) # log likelihood of non "1"s
  return(l)
}

#   -(n1*log(w + (1-w)*(la/(exp(la)-1))) + (n0-n1)*log(1-w) + (sum(y)-n1)*log(la) - 
#     (n0-n1)*log(exp(la)-1) - sum(log(factorial(y))))

llpp <- function(lambda, data = y) {
  l <- sum(data*log(lambda) - log(exp(lambda) - 1) - log(factorial(data)))
  return(l)
}

llotpp <- function(theta, data = y) {
  data <- data[data != 1]
  l <- sum(data*log(theta) - log(exp(theta) - theta - 1) - log(factorial(data)))
  return(l)
}

#OPT functions for optim/numerical maximization

opt_lloipp <- function(param, data = y) {
  return(-lloipp(param, data))
} 

opt_llpp <- function(lambda, data = y) {
  return(-llpp(lambda, data))
} 

opt_llotpp <- function(theta, data = y) {
  return(-llotpp(theta, data))
} 

#ml fitting

oipp_ml_fit <- function(data = y, param = c(2, 0.5), hess = FALSE) {
  ml <- optim(param, opt_lloipp, data = data, lower = c(1e-3, 0.01), upper = c(Inf, 0.99), hessian = hess, method = "L-BFGS-B")
  return(ml)
}

pp_ml_fit <- function(data = y, param = 2, hess = FALSE) {
  ml <- optim(param, opt_llpp, data = data, lower = 1e-3, hessian = hess, method = "L-BFGS-B")
  return(ml)
}

otpp_ml_fit <- function(data = y, param = 2) {
  ml <- optim(2, opt_llotpp, data = data, lower = 1e-3, hessian = TRUE, method = "L-BFGS-B")
  return(list(par = ml$par, se = as.numeric(sqrt(1 / ml$hessian))))
}

# Population estimation

N_oipp <- function(param_ml_optim, data = y) {
  theta_ml <- param_ml_optim
  p_ml <- param_ml_optim
  n <- length(data)
  N_ml <- (n/((1 + p_ml)*(1 - dpois(0, theta_ml)))) %>%
    round()
  return(N_ml)
}

N_pp <- function(param_ml_optim, data = y) {
  n <- length(data)
  N_ml <- (n/(1 - dpois(0, param_ml_optim))) %>%
    round()
  return(N_ml)
}

N_otpp <- function(param_ml_optim, data = y) {
  n <- length(data) - sum(data==1)
  N_ml <- (n/(1 - dpois(0, param_ml_optim) - dpois(1, param_ml_optim))) %>%
    round()
  return(N_ml)
}

# Simulation

sim_oipp <- function(lambda, p, m, n) {
  simulations <- data.frame(N = rep(n, m)) %>%
    rowwise() %>% 
    mutate(sim = list(roipp(N, lambda, p)),
           fit = list(oipp_ml_fit(data = sim)),
           N_ml = N_oipp(fit$par))
  return(simulations)
}

sim_pp <- function(lambda, p, m, n) {
  simulations <- data.frame(N = rep(n, m)) %>%
    rowwise() %>% 
    mutate(sim = list(roipp(N, lambda, p)),
           fit = pp_ml_fit(data = sim),
           N_ml = N_pp(fit$par))
  return(simulations)
}

sim_otpp <- function(lambda, p, m, n) {
  simulations <- data.frame(N = rep(n, m)) %>%
    rowwise() %>% 
    mutate(sim = list(tibble(y = roipp(N, lambda, p))),
           fit = list(otpp_ml_fit(data = sim$y)),
           N_ml = N_otpp(fit[[1]][1]))
  return(simulations)
}

sim_otpp <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = roipp(N, lambda, p))), # Simulate data
           fit = list(otpp_ml_fit(data$y)), # Fit the MLE
           N_est = N_otpp(fit$par, data$y)
           )
  return(simulations)
}






#TESTER

y <- roipp(1903, 8, 0.6)
oipp_ml_fit()$par


y <- roipp(1903, 14, 0.001)
par <- pp_ml_fit()$par
N_pp(par)
par <- otpp_ml_fit()$par
N_otpp(par)

y <- roipp(1903, 14, 0.9)
par <- otpp_ml_fit()[[1]][1]
N_otpp(par)

test <- sim_otpp(n = 1000, lambda = 1:5, p = seq(from = 0, to = 0.5, by = 0.25), n_sim = 10)




