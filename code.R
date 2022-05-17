library(tidyverse)

# Data generationg functions, parametrisations as found in report.

rztp <- function(n, lambda) {
  y <- rpois(n, lambda)
  y <- y[y!=0]
  return(y)
}

rztoip <- function(n, lambda, p) {
  y <- rztp(n, lambda) %>% # generate sample without inflation
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

llztoip <- function(param, data = y){
  lambda <- param[1]
  p <- param[2]
  l <- sum(log((dpois(data, lambda * (1 - p)) + lambda * p * (data == 1)) / (lambda * p + 1 - dpois(0, lambda * (1 - p)))))
  return(l)
}

llztp <- function(lambda, data = y) {
  l <- sum(data*log(lambda) - log(exp(lambda) - 1) - log(factorial(data)))
  return(l)
}

llzotp <- function(theta, data = y) {
  data <- data[data != 1]
  l <- sum(data*log(theta) - log(exp(theta) - theta - 1) - log(factorial(data)))
  return(l)
}

llztoinb <- function(param, data = y) {
  size <- param[1]*param[2] #lambda * k
  prob <- param[2]/(1 + param[2] - param[3]) # k/(1 - p + k)
  omega <- 1/(1 + param[1]*param[3]/(1 - dnbinom(0, size, prob)))
  l = sum(log((1 - omega) + omega*dnbinom(data[data==1], size, prob)/(1 - prob^size))) + # ll of ones
    sum(log(omega*dnbinom(data[data!=1], size, prob)/(1 - prob^size))) # ll of non ones
  return(l)
}

llztnb <- function(param, data = y) {
  size <- param[1]*param[2] #lambda * k
  prob <- param[2]/(1 + param[2]) # k/(1 - p + k), p = 0
  l <- sum(log(dnbinom(data, size, prob)/(1 - prob^size)))
  return(l)
}

llzotnb <- function(param, data = y) {
  data_sub <- data[data != 1]
  size <- param[1]*param[2] #lambda * k
  prob <- param[2]/(1 + param[2]) # k/(1 - p + k), p = 0
  l <- sum(log(dnbinom(data_sub, size, prob)/(1 - prob^size - size*(1 - prob)*prob^size)))
  return(l)
}

#OPT functions for optim/numerical maximization

opt_llztoip <- function(param, data = y) {
  if (0 < param[1] & 0 <= param[2] & param[2] <= 1) { #likes to try values which creates error in likelihood even thought it is not supposed to see lower in ztoip_ml_fit
    fit <- -llztoip(param, data)
    if (abs(fit) == Inf) {
      return(1e+100)
    } else {
      return(fit)
    }
  } else {
    return(1e+100) # optim avoids non allowed values
  }
} 

opt_llztp <- function(lambda, data = y) {
  if (0 < lambda) { 
    fit <- -llztp(lambda, data)
    if (abs(fit) == Inf) {
      return(1e+100)
    } else {
      return(fit)
    }
  } else {
    return(1e+100)
  }
} 

opt_llzotp <- function(theta, data = y) {
  if (0 < theta) { 
    fit <- -llzotp(theta, data)
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

ztoip_ml_fit <- function(data = y, param = c(2, 0.5)) {
  ml <- optim(par = param, fn = opt_llztoip, data = data, lower = c(0, 0), upper = c(Inf, 1), hessian = TRUE, method = "L-BFGS-B") # use method to allow limits
  return(list(par = ml$par, hess = ml$hessian))
}

ztp_ml_fit <- function(data = y, param = 2) {
  ml <- optim(par = param, fn = opt_llztp, data = data, hessian = TRUE, lower = 0, method = "L-BFGS-B") 
  return(list(par = ml$par, se = as.numeric(sqrt(1 / ml$hessian))))
}

zotp_ml_fit <- function(data = y, param = 2) {
  ml <- optim(par = param, fn = opt_llzotp, data = data, lower = 0, hessian = TRUE, method = "L-BFGS-B")
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

N_otp <- function(param, data = y) {
  n <- length(data)
  N <- (n/(1 - dpois(0, param))) %>%
    round()
  return(N)
}

N_p <- function(param, data = y) {
  param <- c(param, 0)
  theta <- param[1]*(1 - param[2])
  n <- length(data) - sum(data==1)
  N <- (n/(1 - dpois(0, theta) - dpois(1, theta))) %>%
    round()
  return(N)
}

N_ztnb <- function(param, data = y) {
  size <- param[1]*param[2] #lambda * k
  prob <- param[2]/(1 + param[2]) # k/(1 - p + k), p = 0
  n <- length(data)
  N <- (n/(1 - dnbinom(0, size, prob))) %>%
    round()
  return(N)
}

N_nb <- function(param, data = y) {
  param <- c(param, 0)
  size <- param[1]*param[2] #lambda * k
  prob <- param[2]/(1 + param[2] - param[3]) # k/(1 - p + k)
  n <- length(data) - sum(data==1)
  N <- (n/(1 - dnbinom(0, size, prob) - dnbinom(1, size, prob))) %>%
    round()
  return(N)
}

# Simulation

sim_oiztp <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoip(N, lambda, p))), # Simulate data
           fit = list(ztoip_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_ztoip(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_ztp <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoip(N, lambda, p))), # Simulate data
           fit = list(pp_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_ztp(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_zotp <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoip(N, lambda, p))), # Simulate data
           fit = list(zotp_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_p(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_p <- function(n, lambda, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoip(N, lambda, p))), # Simulate data
           fit_ztoip = list(ztoip_ml_fit(data = data$y)), # Fit the MLE
           fit_ztp = list(ztp_ml_fit(data = data$y)),
           fit_zotp = list(zotp_ml_fit(data = data$y)),
           N_est_ztoip = N_p(par = fit_ztoip$par, data = data$y), # Estimate population
           N_est_ztp = N_ztp(par = fit_ztp$par, data = data$y),
           N_est_zotp = N_p(par = fit_zotp$par, data = data$y))
  return(simulations)
}

sim_ztoinb <- function(n, lambda, k, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, k = k, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoinb(N, lambda, k, p))), # Simulate data
           fit = list(ztoinb_ml_fit(data = data$y)), # Fit the MLE
           N_est = N_nb(par = fit$par, data = data$y)) # Estimate population
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
           N_est = N_nb(par = fit$par, data = data$y)) # Estimate population
  return(simulations)
}

sim_nb <- function(n, lambda, k, p, n_sim) {
  simulations <- expand_grid(lambda = lambda, k = k, p = p, N = n, sim = 1:n_sim) %>% 
    rowwise() %>% 
    mutate(data = list(tibble(y = rztoinb(N, lambda, k, p))), # Simulate data
           fit_ztoinb = list(ztoinb_ml_fit(data = data$y)), # Fit the MLE
           fit_ztnb = list(ztnb_ml_fit(data = data$y)),
           fit_zotnb = list(zotnb_ml_fit(data = data$y)),
           N_est_ztoinb = N_nb(par = fit_ztoinb$par, data = data$y), # Estimate population
           N_est_ztnb = N_ztnb(par = fit_ztnb$par, data = data$y),
           N_est_zotnb = N_nb(par = fit_zotnb$par, data = data$y))
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

data_captures <- read_csv("captures.csv")

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
# P_SIM <- sim_p(c(500, 1500, 2500), c(1,2, 3), c(0, 0.05, 0.1), 1000)
P_SIM <- readRDS("P_SIM")[[1]]

# set.seed(2844707)
# PP_SIM <- sim_pp(c(500, 2500), c(1,2,3), seq(from = 0, to = 0.4, by = 0.1), 100)
PP_SIM <- readRDS("PP_SIM")[[1]]

# set.seed(6556485)
# NB_SIM <- sim_nb(c(500,1500,2500), c(1,2,3), c(0.2,0.6,1),  c(0, 0.05, 0.1), 1000)
NB_SIM <- rbind(readRDS("NB_SIM_part_1")[[1]], readRDS("NB_SIM_part_2")[[1]])

# set.seed(0911764)
# NB_SIM_ext <- sim_zotnb(500,c(1,2,3), 2, 0.1, 1000)
NB_SIM_ext <- readRDS("NB_SIM_ext")[[1]]
