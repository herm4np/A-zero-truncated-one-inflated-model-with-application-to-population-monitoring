library(tidyverse)
library(purrr)

data_captures <- read_csv("captures_copy.csv")

rOIZTP <- function(n, p, lambda) { #sample from one-inflated zero-trunkated poisson
  sample <- rpois(n, lambda*(1 - p)) #sample from poisson
  sample = sample + 1 #adjust for zero trunkated
  return(sample)
}

rOIZTP_w_ghosts <- function(n, lambda, p) {
  sample_n_ghosts <- rpois(n, lambda) %>% #sample observations of each bear
    map_df(~data.frame(n_obs = .x, ghosts = rbinom(1, .x, p))) %>% #calculate n of ghosts per observed individual
    filter(n_obs != 0) %>% #remove unobserved bears
    mutate(correct_obs = n_obs - ghosts) #add column with amount of correct observations
  return(sample_n_ghosts)
}

rOIZTNB_w_ghosts <- function(n, size, prob, p) { #sample from one-inflated zero-trunkated negative binomial
  sample_n_ghosts <- rnbinom(n, size, prob) %>% #sample observations of each bear
    map_df(~data.frame(n_obs = .x, ghosts = rbinom(1, .x, p))) %>% #calculate n of ghosts per observed individual
    filter(n_obs != 0) %>% #remove unobserved bears
    mutate(correct_obs = n_obs - ghosts) #add column with amount of correct observations
  return(sample_n_ghosts)
}

data_t <- data_captures %>% 
  mutate(simulated = FALSE) %>%
  group_by(id, simulated) %>%
  summarise(n_obs = n()) %>%
  full_join(rOIZTNB_count_ghosts(4000, 6, 0.6, 0.001) %>% mutate(id = str_c("id", as.character(row_number())), simulated = TRUE),
            by = c("id", "simulated", "n_obs")) 

data_t %>%
  ggplot(aes(x = n_obs)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~simulated)

lloipp_nocov_TEST <- function(la, theta, data){
  n0 <- length(data)
  n1 <- sum(data==1)
  w = 1/(1+exp(-theta))
  n1*log(w + (1-w)*(la/(exp(la)-1))) + (n0-n1)*log(1-w) + (sum(data)-n1)*log(la) - 
    (n0-n1)*log(exp(la)-1) - sum(log(factorial(data)))
}

data_antal_osb <- data_captures %>% 
  group_by(id) %>%
  summarise(n_obs = n())

opt_lloipp_nocov <- function(param){
  la = param[1]
  theta = param[2]
  n0 <- length(y)
  n1 <- sum(y==1)
  w = 1/(1+exp(-theta))
  -(n1*log(w + (1-w)*(la/(exp(la)-1))) + (n0-n1)*log(1-w) + (sum(y)-n1)*log(la) - 
    (n0-n1)*log(exp(la)-1) - sum(log(factorial(y))))
}

opt_lloiztnb_nocov <- function(param) {
  li <- param[1]
  a  <- param[2]
  
  # Use the logistic link (Section 2.1 of the main paper)
  phi <- param[3]
  w <- 1 / (1 + exp(phi))
  
  ymax <- max(y)
  terms <- weights <- rep(0, ymax)
  for(ii in 1:ymax) {
    terms[ii] <- log(a + ii - 1)
  }
  weights[1] <- sum(y > 1)
  for(ii in 2:ymax) {
    weights[ii] <- sum(y > (ii - 1))
  }
  gterm <- terms * weights
  return(-(sum(log(1 - w) + (y == 1) * (log(w / (1 - w) + a * ((a / (a + li)) ^ a) 
                                          * (li / (a + li - a * (1 + (li / a)) ^ (1 - a))))) + (1 - (y == 1)) * 
               (a * log(a) - log(factorial(y)) + y * log(li) - (a + y) * log(a + li) - 
                  log(1 - (a / (a + li)) ^ a))) + sum(gterm)))
}

