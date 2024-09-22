
# simulate one big dataset and fit the correct model

post_path <- fs::path(getwd(), "posts", "simulating-survival-data-with-non-PH-and-cure") %>%
  print()

library(tidyverse)
library(flexsurvcure)

n_sim <- 50
n_sample <- 1000
n_total <- n_sim * n_sample

set.seed(43)

simulated_data <- tibble(
  #
  # simulate id variables and covariates
  #
  sim_id = rep(1:n_sim, each = n_sample),
  id = rep(1:n_sample, times = n_sim),
  treatment = sample(c(0, 1), size = n_total,
                     replace = TRUE, prob = c(0.5, 0.5)),
  parameters = case_when(
    treatment == 1 ~ "Treatment group (T ~ LL(4, 7))",
    treatment == 0 ~ "Control group (T ~ LL(5, 7))"
  ),
  age = {
    rgamma(n = n_total, shape = 40, rate = 0.7) %>%
      pmin(., 90) %>%
      pmax(20, .) %>%
      round(., 1)
  },
  #
  # simulate latent event times as a function of covariates
  #
  linear_predictor_alpha = log(5) - log(1.25) * treatment,
  linear_predictor_beta =log(7),
  u = runif(n_total, min = 0, max = 1),
  uncured_time = (exp(linear_predictor_alpha) *
                    (((1 / u) - 1) ^ (-1 / exp(linear_predictor_beta)))),
  #
  # simulate the cure / incidence part from a logistic model
  #
  cured = rbinom(n_total, 1, plogis(4 - 0.07 * age)),
  latent_time = ifelse(cured == 1, 10000, uncured_time),
  #
  # simulate censoring times and censor some individuals
  #
  random_censoring = rgamma(n_total, shape = 6, scale = 1),
  # keep the smallest of the random censoring
  # and administrative censoring (at t = 6) times
  censoring_time = pmin(random_censoring, 6),
  time = pmin(latent_time, censoring_time),
  cens = as.numeric(latent_time > censoring_time),
  event = 1 - cens
)

simulated_data %>% glimpse(width = 80)

# fit the correct model
model <- simulated_data %>%
  # rename the treatment variable to make the printed summary nicer looking
  rename(trt = treatment) %>%
  flexsurvcure(
    formula = Surv(time, event, type = "right") ~ age,
    data = .,
    dist = "llogis",
    link = "logistic",
    anc = list(
      scale = ~ trt,
      shape = ~ trt
    ),
    mixture = TRUE
  )

model

# persist model to disk
write_rds(x = model, file = fs::path(post_path, "sanity-check-model.rds"))

# read it back to see if it works correctly
read_rds(file = fs::path(post_path, "sanity-check-model.rds"))
