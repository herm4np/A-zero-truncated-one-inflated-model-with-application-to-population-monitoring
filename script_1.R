library(tidyverse)

# fördelningar

rpp <- function(n, lambda) {
  y <- rpois(n, lambda)
  y <- y[y!=0]
  return(y)
}

roipp <- function(n, lambda, p) {
  y <- rpp(n, lambda) %>% # generate sample without inflation
    map_df(~data.frame(y = .x, 
                       real = rbinom(1, .x, (1 - p)))) %>% # generate inflation
    mutate(ghosts = y - real)
  ghosts = sum(y["ghosts"]) # count total ghosts
  y <- y %>% 
    filter(real != 0) %>%
    select(real) %>%
    data.matrix %>%
    array()
  y = c(y, rep.int(1, ghosts))
  return(y)
}

rztnb <- function(n, lambda, k) {
  y <- rnbinom(n, size = lambda*k, prob = 1 - 1/(1 + k))
  y <- y[y!=0]
  return(y)
}

rztoinb <- function(n, lambda, k, p) {
  y <- rztnb(n, lambda, k) %>% # generate sample without inflation
    map_df(~data.frame(y = .x, 
                       real = rbinom(1, .x, (1 - p)))) %>% # generate inflation
    mutate(ghosts = y - real)
  ghosts = sum(y["ghosts"]) # count total ghosts
  y <- y %>% 
    filter(real != 0) %>%
    select(real) %>%
    data.matrix %>%
    array()
  y = c(y, rep.int(1, ghosts))
  return(y)
}

# log likelihood computing functions

lloipp <- function(param, data = y){
  lambda <- param[1]
  p <- param[2]
  l <- sum(log((dpois(data, lambda * (1 - p)) + lambda * p * (data == 1)) / (lambda * p + 1 - dpois(0, lambda * (1 - p)))))
  return(l)
}

llpp <- function(lambda, data = y) {
  l <- sum(data*log(lambda) - log(exp(lambda) - 1) - log(factorial(data)))
  return(l)
}

llotpp <- function(theta, data = y) {
  data <- data[data != 1]
  l <- sum(data*log(theta) - log(exp(theta) - theta - 1) - log(factorial(data)))
  return(l)
}

llztoinb <- function(param, data = y) {
  lambda <- param[1]
  p <- param[3]
  size <- lambda*param[2] #lambda * k
  prob <- 1 - 1/(1 + param[2]/(1 - p)) # 1/(1 + k/(1 - p))
  omega <- 1/(1 + lambda*p/(1 - dnbinom(0, size, prob)))
  l = sum(log((1 - omega) + omega*dnbinom(data[data==1], size, prob)/(1 - prob^size))) + # ll of ones
    sum(log(omega*dnbinom(data[data!=1], size, prob)/(1 - prob^size))) # ll of non ones
  return(l)
}

llztnb <- function(param, data = y) {
  size <- param[1]*param[2] #lambda * k
  prob <- 1 - 1/(1 + param[2]) # 1/(1 + k/(1 - p)), p = 0
  l <- sum(log(dnbinom(data, size, prob)/(1 - prob^size)))
  return(l)
}

llzotnb <- function(param, data = y) {
  data_sub <- data[data != 1]
  size <- param[1]*param[2] #lambda * k
  prob <- 1 - 1/(1 + param[2]) # 1/(1 + k/(1 - p)), p = 0
  l <- sum(log(dnbinom(data_sub, size, prob)/(1 - prob^size - size*(1 - prob)*prob^size)))
  return(l)
}

#OPT functions for optim/numerical maximization

opt_lloipp <- function(param, data = y) {
  if (0 < param[1] & 0 <= param[2] & param[2] <= 1) { #likes to try values which creates error in likelihood even thought it is not supposed to see lower in oipp_ml_fit
    fit <- -lloipp(param, data)
    if (abs(fit) == Inf) {
      return(1e+100)
    } else {
      return(fit)
    }
  } else {
    return(1e+100) # optim avoids non allowed values
  }
} 

opt_llpp <- function(lambda, data = y) {
  if (0 < lambda) { 
    fit <- -llpp(lambda, data)
    if (abs(fit) == Inf) {
      return(1e+100)
    } else {
      return(fit)
    }
  } else {
    return(1e+100)
  }
} 

opt_llotpp <- function(theta, data = y) {
  if (0 < theta) { 
    fit <- -llotpp(theta, data)
    if (abs(fit) == Inf) {
      return(1e+100)
    } else {return(fit)}
  } else {
    return(1e+100)
  }
} 

opt_llztoinb <- function(param, data = y) {
  if (0 < param[1] & 0 < param[2] & 0 <= param[3] & param[3] <= 1) { 
    fit <- -llztoinb(param, data)
    if (is.nan(fit) | abs(fit) == Inf) {
      return(1e+100)
    } else {return(fit)}
  } else {
    return(1e+100)
  }
}

opt_llztnb <- function(param, data = y) {
  if (0 < param[1] & 0 < param[2]) { 
    fit <- -llztnb(param, data)
    if (is.nan(fit) | abs(fit) == Inf) {
      return(1e+100)
    } else {return(fit)}
  } else {
    return(1e+100)
  }
}

opt_llzotnb <- function(param, data = y) {
  if (0 < param[1] & 0 < param[2]) { 
    fit <- -llzotnb(param, data)
    if (is.nan(fit) | abs(fit) == Inf) {
      return(1e+100)
      } else {return(fit)}
    }
  else {
    return(1e+100)
  }
}

#ml fitting

oipp_ml_fit <- function(data = y, param = c(2, 0.5)) {
  ml <- optim(par = param, fn = opt_lloipp, data = data, lower = c(0, 0), upper = c(Inf, 1), hessian = TRUE, method = "L-BFGS-B") # use method to allow limits
  return(list(par = ml$par, hess = ml$hessian))
}

pp_ml_fit <- function(data = y, param = 2) {
  ml <- optim(par = param, fn = opt_llpp, data = data, hessian = TRUE, lower = 0, method = "L-BFGS-B") 
  return(list(par = ml$par, se = as.numeric(sqrt(1 / ml$hessian))))
}

otpp_ml_fit <- function(data = y, param = 2) {
  ml <- optim(par = param, fn = opt_llotpp, data = data, lower = 0, hessian = TRUE, method = "L-BFGS-B")
  return(list(par = ml$par, se = as.numeric(sqrt(1 / ml$hessian))))
}

ztoinb_ml_fit <- function(data = y, param = c(2, 2, 0.5)) {
  ml <- optim(par = param, fn = opt_llztoinb, data = data, hessian = TRUE, lower = c(0, 0, 0), upper = c(Inf, Inf, 1), method = "L-BFGS-B")
  return(list(par = ml$par, hess = ml$hessian))
}

ztnb_ml_fit <- function(data = y, param = c(2, 2)) {
  ml <- optim(par = param, fn = opt_llztnb, data = data, hessian = TRUE, lower = c(0, 0), upper = c(Inf, Inf), method = "L-BFGS-B") # use method to allow limits
  return(list(par = ml$par, hess = ml$hessian))
}

zotnb_ml_fit <- function(data = y, param = c(2, 2)) {
  ml <- optim(par = param, fn = opt_llzotnb, data = data, hessian = TRUE, lower = c(0, 0), upper = c(Inf, Inf), method = "L-BFGS-B") # use method to allow limits
  return(list(par = ml$par, hess = ml$hessian))
}

# Population estimation

N_oipp <- function(param, data = y) {
  lambda <- param[1]
  p <- param[2]
  theta <- lambda*(1 - p)
  n <- length(data) - sum(data==1)*(p*lambda/(p*lambda + (dpois(1, theta)))) #subtract estimated proportion of "ghosts"
  N <- (n/(1 - dpois(0, theta))) %>%
    round()
  return(N)
}

N_pp <- function(param, data = y) {
  n <- length(data)
  N <- (n/(1 - dpois(0, param))) %>%
    round()
  return(N)
}

N_otpp <- function(param, data = y) {
  n <- length(data) - sum(data==1)
  N <- (n/(1 - dpois(0, param) - dpois(1, param))) %>%
    round()
  return(N)
}

N_ztoinb <- function(param, data = y) {
  lambda <- param[1]
  p <- param[3]
  size <- lambda*param[2] #lambda * k
  prob <- 1 - 1/(1 + param[2]/(1 - p)) # 1/(1 + k/(1 - p))
  n <- length(data) - sum(data==1)*p*lambda/(p*lambda + dnbinom(1, size, prob)) #subtract estimated proportion of "ghosts"
  N <- (n/(1 - dnbinom(0, size, prob))) %>%
    round()
  return(N)
}

N_ztnb <- function(param, data = y) {
  size <- param[1]*param[2] #lambda * k
  prob <- 1 - 1/(1 + param[2]) # 1/(1 + k/(1 - p)), p = 0
  n <- length(data)
  N <- (n/(1 - dnbinom(0, size, prob))) %>%
    round()
  return(N)
}

N_zotnb <- function(param, data = y) {
  size <- param[1]*param[2] #lambda * k
  prob <- 1 - 1/(1 + param[2]) # 1/(1 + k/(1 - p)), p = 0
  n <- length(data) - sum(data==1)
  N <- (n/(1 - dnbinom(0, size, prob) - dnbinom(1, size, prob))) %>%
    round()
  return(N)
}

# Simulation

sim_oipp <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = roipp(N, lambda, p))), # Simulate data
           fit = list(oipp_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_oipp(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_pp <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = roipp(N, lambda, p))), # Simulate data
           fit = list(pp_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_pp(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_otpp <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = roipp(N, lambda, p))), # Simulate data
           fit = list(otpp_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_otpp(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_p <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = roipp(N, lambda, p))), # Simulate data
           fit_oipp = list(oipp_ml_fit(data = data$y)), # Fit the MLE
           fit_pp = list(pp_ml_fit(data = data$y)),
           fit_otpp = list(otpp_ml_fit(data = data$y)),
           N_est_oipp = N_oipp(par = fit_oipp$par, data = data$y), # Estimate population
           N_est_pp = N_pp(par = fit_pp$par, data = data$y),
           N_est_otpp = N_otpp(par = fit_otpp$par, data = data$y))
  return(simulations)
}

sim_ztoinb <- function(n, lambda, k, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, k = k, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoinb(N, lambda, k, p))), # Simulate data
           fit = list(ztoinb_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_ztoinb(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_ztnb <- function(n, lambda, k, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, k = k, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoinb(N, lambda, k, p))), # Simulate data
           fit = list(ztnb_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_ztnb(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_zotnb <- function(n, lambda, k, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, k = k, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoinb(N, lambda, k, p))), # Simulate data
           fit = list(zotnb_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_zotnb(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_nb <- function(n, lambda, k, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, k = k, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoinb(N, lambda, k, p))), # Simulate data
           fit_ztoinb = list(ztoinb_ml_fit(data = data$y)), # Fit the MLE
           fit_ztnb = list(ztnb_ml_fit(data = data$y)),
           fit_zotnb = list(zotnb_ml_fit(data = data$y)),
           N_est_ztoinb = N_ztoinb(par = fit_ztoinb$par, data = data$y), # Estimate population
           N_est_ztnb = N_ztnb(par = fit_ztnb$par, data = data$y),
           N_est_zotnb = N_zotnb(par = fit_zotnb$par, data = data$y))
  return(simulations)
}

#

filter_obs <- function(data, n) {
  return(data %>% filter(n_obs < n))
}

qqplot_gen <- function(data, type = 0) { # 0 for oiztnb estimate, 1 for zotnb
  data_obs <- data$n_obs
  par <- ztoinb_ml_fit(data_obs)$par
  lambda <- par[1]
  p <- par[3]
  size <- lambda*par[2] #lambda * k
  prob <- 1/(1 + par[2]/(1 - p)) # 1/(1 + k/(1 - p))
  data_w_removed_ghosts <- c(data_obs[data_obs!=1], # remove all ones
                             rep.int(1, round(sum(data_obs==1)*(1 - (p*lambda/(p*lambda + dnbinom(1, size, prob))))))) #add estimated true value of 1s
  nbin_q <- qnbinom(ppoints(length(data_w_removed_ghosts)), size, prob)
  qqplot <- qqplot(data_w_removed_ghosts, nbin_q[nbin_q != 0])
  return(qqplot)
}

qqplot_gen_2 <- function(data, type = 0) { # 0 for oiztnb estimate, 1 for zotnb
  data_obs <- data$n_obs
  par <- zotnb_ml_fit(data_obs)$par
  lambda <- par[1]
  size <- lambda*par[2] #lambda * k
  prob <- 1 - 1/(1 + par[2]) # 1/(1 + k/(1 - p))
  data_w_removed_ghosts <- data_obs[data_obs!=1] # remove all ones
  nbin_q <- qnbinom(ppoints(1000000), size, prob)
  return(nbin_q[nbin_q >= 1])
}

# code for EX.Rmd

data_captures <- read_csv("captures_copy.csv")

data_summary <- data_captures %>%
  group_by(id) %>%
  summarise(n_obs = n())

data_summary_2015 <- data_captures %>%
  filter(year == 2015) %>%
  group_by(id) %>%
  summarise(n_obs = n())

data_summary_2016 <- data_captures %>%
  filter(year == 2016) %>%
  group_by(id) %>%
  summarise(n_obs = n())

data_summary_2017 <- data_captures %>%
  filter(year == 2017) %>%
  group_by(id) %>%
  summarise(n_obs = n())

data_summary_2019 <- data_captures %>%
  filter(year == 2019) %>%
  group_by(id) %>%
  summarise(n_obs = n())

data_summary_2020 <- data_captures %>%
  filter(year == 2020) %>%
  group_by(id) %>%
  summarise(n_obs = n())

# set.seed(5118761)
# P_SIM <- sim_p(c(500, 900, 1300), c(1,3,5), c(0, 0.1, 0.2), 1000)
# set.seed(5118761)
# P_SIM <- sim_p(c(500, 1500, 2500), c(1,2, 3), c(0, 0.05, 0.1), 1000)
P_SIM <- readRDS("P_SIM2")[[1]]
# set.seed(2844707)
# PP_SIM <- sim_pp(c(500,900,1300), c(1,3,5), c(0, 0.025, 0.05), 1000)
PP_SIM <- readRDS("PP_SIM")[[1]]
# set.seed(6556485)
# NB_SIM <- sim_nb(c(500,1500,2500), c(1,2,3), c(0.2,0.6,1),  c(0, 0.05, 0.1), 1000)
NB_SIM <- readRDS("NB_SIM2")[[1]]





NB_SIM %>%
  filter() %>%
  group_by(lambda, p, N, k) %>%
  summarize(avg_N_est_percent_error = mean(N_est_ztoinb)/N) %>%
  filter(N == 1500) %>%
  ggplot(aes(x = p, y = avg_N_est_percent_error, group = as.factor(lambda), color = as.factor(lambda))) +
  geom_line(alpha = 0.75) +
  scale_color_manual(values=c("red", "blue", "green")) +
  scale_x_continuous(breaks = c(0, 0.1)) +
  facet_wrap(~ k, scales = "free_y")

sim_nb <- function(n, lambda, k, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, k = k, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoinb(N, lambda, k, p))))
  return(simulations)
}
















### UPPDATERAD


# 
# 
# 
# 
# NB_SIM %>%
#   filter(lambda != 0.5) %>%
#   group_by(lambda, p, k) %>%
#   summarize(avg_N_est_percent_error = mean(N_est_ztoinb)/N) %>%
#   ggplot(aes(x = p, y = avg_N_est_percent_error, group = as.factor(lambda), color = as.factor(lambda))) +
#   geom_line() +
#   scale_color_manual(values=c("red", "blue", "green")) +
#   facet_grid(~k)
# 
# t2t2 <- sim_zotnb(1000, c(3,6), c(0.5,1),  c(0, 0.1, 0.2, 0.3), 100)
# 
# 
# t2t2%>%
#   group_by(lambda, p, k) %>%
#   summarize(avg_N_est_percent_error = mean(N_est)/N) %>%
#   ggplot(aes(x = p, y = avg_N_est_percent_error, group = as.factor(lambda), color = as.factor(lambda))) +
#   geom_line() +
#   scale_color_manual(values=c("red", "blue", "green")) +
#   facet_grid(~k)
# 
# 
# t11t1t1 <- sim_oipp(900, 1, seq(from = 0, to = 0.9, by = 0.1), 100)
# 
# 
# t11t1t1 %>%
#   group_by(lambda, p, N) %>%
#   filter(p < 0.8) %>%
#   summarize(avg_N_est_percent_error = mean(N_est)/N) %>%
#   ggplot(aes(x = p, y = avg_N_est_percent_error, group = as.factor(lambda), color = as.factor(lambda))) +
#   geom_line() +
#   scale_color_manual(values=c("red", "blue", "green"))
# 
# ztoinb_ml_fit(data_summary_2020$n_obs)
# 
# 
# NB_SIM %>%
#   filter(k == 1, lambda != 0.5, p != 0.2) %>%
#   group_by(lambda, p, N) %>%
#   summarize(avg_N_est_percent_error = mean(N_est_zotnb)/N) %>%
#   ggplot(aes(x = p, y = avg_N_est_percent_error, group = as.factor(lambda), color = as.factor(lambda))) +
#   geom_line() +
#   scale_color_manual(values=c("red", "blue", "green")) +
#   facet_grid(~N) 
# 
# NB_SIM
# 
# saveRDS(list(NB_SIM), "NB_SIM")




# test

# y <- roipp(2000, 8, 0)
# oipp_ml_fit()
# 
# 
# y <- roipp(1903, 14, 0.001)
# par <- pp_ml_fit()$par
# N_pp(par)
# par <- otpp_ml_fit()$par
# N_otpp(par)
# 
# y <- roipp(1903, 14, 0.9)
# par <- otpp_ml_fit()[[1]][1]
# N_otpp(par)
# 
# test <- sim_oipp(n = 1000, lambda = seq(from = 1, to = 9, by = 4), p = seq(from = 0, to = 0.1, by = 0.05), n_sim = 100)
# 
# test <- sim_pp(n = 100, lambda = seq(from = 1, to = 10, by = 5), p = seq(from = 0, to = 0.1, by = 0.05), n_sim = 3)
# 
# test <- sim_otpp(n = 1000, lambda = 1:5, p = seq(from = 0, to = 0.5, by = 0.25), n_sim = 10)
# 
# k = 10
# lambda = 20
# p = 0.1
# data.frame(dat = rbinom(10000, k*lambda, 1/(1 + k/(1 - p)))) %>%
#   ggplot(aes(x = dat)) +
#   geom_histogram(binwidth = 1)
# 
# 
# dat <- rztoinb(10000, 2, 1, 0.05)
# N_ztoinb(ztoinb_ml_fit(dat)$par, dat)
# 
# test <- data_captures %>% 
#   select(year, id) %>%
#   distinct() %>%
#   group_by(id) %>%
#   summarise(n())
# 
# 
# 
# 
# 
# 
# fit <- ztoinb_ml_fit(data_summary_2015$n_obs)
# lambda <- fit$par[1]
# p <- fit$par[3]
# size <- lambda*fit$par[2] #lambda * k
# prob <- 1/(1 + fit$par[2]/(1 - p)) # 1/(1 + k/(1 - p))
# omega <- 1/(1 + lambda*p)

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
#
#
# doipp <- function(x, lambda, p){
#   return((dpois(x, lambda * (1 - p)) + lambda * p * (x == 1)) / (lambda * p + 1 - dpois(0, lambda * (1 - p))))
# }
# 
# lloipp <- function(param, data = y){
#   return(sum(log(doipp(data, param[1], param[2]))))
# }
# 
# lloipp <- function(param, data = y) {
#   lambda <- param[1]
#   p <- param[2]
#   theta <- lambda*(1 - p)
#   omega <- 1/(1 + lambda*p)
#   n1 <- sum(data==1)
#   nm <- length(data) - n1
#   sum_data_m <- sum(data) - n1
#   l <- n1*log((1 - omega) + omega*(theta/(exp(theta) - 1))) + #log likelihood of "1" observations
#     nm*log(omega) + sum_data_m*log(theta) - nm*log(exp(theta) - 1) # log likelihood of non "1"s
#   return(l)
# }

