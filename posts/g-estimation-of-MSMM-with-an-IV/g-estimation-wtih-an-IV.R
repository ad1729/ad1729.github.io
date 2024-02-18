
# This file contains all the code written for the
# MSMM IV post. It contains extra code that was not
# added to the post for reasons of brevity.

library(tidyverse)
library(broom, include.only = "tidy")
# functions from the following packages are called via package::fun()
# to reduce namespace conflicts
# library(geeM)
# library(mgcv)
# library(ivtools)
# library(marginaleffects)
# library(glue)

theme_set(theme_bw())

# simulation parameters ----
confounder_values <- seq(0, 1e5, 50) / 1e4

# exposure distribution
pmax(round(rgamma(5e4, shape = 5, scale = 4)), 1) %>% summary()
hist(pmax(round(rgamma(5e4, shape = 5, scale = 4)), 1), breaks = 30)

# # baseline P[Y = 1] when U = 0 and A = 0 : 0.002 = 2 / 1000
# # RD for U: 0.005
# 0.02 + 0.005 * seq(0, 1e5, 10000) / 1e4
# # RD for A: 0.0007 (effect of treatment on outcome)
# 0.02 + 0.005 * seq(0, 1e5, 10000) / 1e4 + 0.0007 * seq(0, 100, 10)

# change so that a single risk ratio is specified
# instead of risk difference
# RR for U: 1.4
round(exp(log(0.002) + log(1.4) * (seq(0, 1e5, 10000) / 1e4)), 4)
# RR for A: 1.02 (effect of treatment on outcome)
round(exp(log(0.002) + log(1.4) * (seq(0, 1e5, 10000) / 1e4) + log(1.02) * seq(0, 80, 8)), 4)
# generating from a logit-binomial model rather than a log-binomial model
round(plogis(log(0.002) + log(1.4) * (seq(0, 1e5, 10000) / 1e4) + log(1.02) * seq(0, 80, 8)), 4)

# simulating a single dataset ----
n <- 1e4
set.seed(24)
df <- tibble(
  id = 1:n,
  # randomization indicator / instrument, indicates group assignment
  Z = rbinom(n, 1, prob = 0.7),
  # simulate continuous confounder value
  U = sample(x = confounder_values, size = n, replace = TRUE),
  # simulate (conditional) compliance indicator
  prob_C = plogis(2 - 0.5 * U),
  C = rbinom(n, 1, prob_C),
  # simulate observed treatment
  # potential continuous exposure value is simulated
  # for each individual on [1, 100]
  # but for those that have Z = 0 (control) or C = 0 (never-takers), A = 0
  A = Z * C * pmin(pmax(round(rgamma(n, shape = 5, scale = 4)), 1), 100),
  # response variable simulated with risk ratios specified
  # baseline log probabilities
  logpY0 = log(0.002) + log(1.4) * U,
  # impact of treatment
  logpY1 = logpY0 + log(1.02) * A,
  # simulated binary responses from log-binomial model
  Y = rbinom(n, 1, exp(logpY1))
) %>%
  glimpse() %>%
  print(n = 5)

# exploring the simulated dataset ----
df %>%
  filter(A > 0) %>%
  ggplot(aes(x = A)) +
  geom_bar() +
  scale_x_continuous(breaks = seq(0, 100, 5), labels = seq(0, 100, 5))

# probability of compliance as a function of the confounder
df %>%
  ggplot(aes(x = U, y = prob_C)) +
  geom_line() + theme_bw() + coord_cartesian(y = c(0, 1))

# plotting the impact of treatment A on the log probability of outcome
# at various levels of U
df %>%
  ggplot(aes(x = A, y = logpY1, color = logpY0)) +
  geom_point() + theme_bw()

df %>%
  ggplot(aes(x = logpY0, y = logpY1, color = A)) +
  geom_point() + theme_bw()

# for binned values of logit(p) and the exposure
# lines should be approximately parallel
df %>%
  mutate(logpY0_bin = cut_interval(logpY0, n = 5)) %>%
  group_by(A, logpY0_bin) %>%
  summarize(m = mean(logpY1), .groups = "drop") %>%
  ungroup() %>%
  ggplot(aes(x = A, y = m, color = logpY0_bin)) +
  geom_point() + geom_smooth(se = FALSE, method = "lm", formula = 'y ~ x') +
  xlab("A") + ylab("log(P[Y|A])") +
  theme_bw()

# proportion of compliers overall independent of group assignment
df %>% count(C) %>% mutate(p = n / sum(n))

# proportion in the prinicipal strata (never-takers, always-takers, compliers, defiers)
# A is dichotomized here to see how many people would comply with treatment if assigned
df %>%
  mutate(A_bin = as.numeric(A > 0)) %>%
  count(Z, A_bin) %>%
  tidyr::complete(Z, A_bin, fill = list(n = 0)) %>%
  mutate(p = n / sum(n))

# check balance of U across Z
# balanced as expected, because Z is a randomization indicator
df %>%
  mutate(Z = factor(Z, levels = c(1, 0), labels = c("Treatment", "Control"))) %>%
  ggplot(aes(x = U, color = Z)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "bottom")

# check balance across (dichotomized version of) A - not balanced because U and A are associated
df %>%
  mutate(`A (binary)` = factor(as.numeric(A > 0), levels = c(1, 0), labels = c("Treated", "Control"))) %>%
  ggplot(aes(x = U, color = `A (binary)`)) +
  geom_density() +
  theme_bw() +
  theme(legend.position = "bottom")

df %>%
  mutate(`A (categorized)` = cut(A, breaks = c(-Inf, 0, 10, 20, 30, 40, 100),
                    right = TRUE, ordered_result = FALSE)) %>%
  ggplot(aes(x = U, color = `A (categorized)`)) +
  geom_density() +
  theme_bw()

# simulate data function ----
simulate_data <- function(seed = 24, n = 1e5) {
  set.seed(seed)
  tibble(
    id = 1:n,
    # randomization indicator / instrument, indicates group assignment
    Z = rbinom(n, 1, prob = 0.7),
    # simulate continuous confounder value
    U = sample(x = confounder_values, size = n, replace = TRUE),
    # simulate (conditional) compliance indicator
    prob_C = plogis(2 - 0.5 * U),
    C = rbinom(n, 1, prob_C),
    # simulate observed treatment
    # potential continuous exposure value is simulated
    # for each individual on [1, 100]
    # but for those that have Z = 0 (control) or C = 0 (never-takers), A = 0
    A = Z * C * pmin(pmax(round(rgamma(n, shape = 5, scale = 4)), 1), 100),
    # response variable simulated with risk ratios specified
    # baseline log probabilities
    logpY0 = log(0.002) + log(1.4) * U,
    # impact of treatment
    logpY1 = logpY0 + log(1.02) * A,
    # simulated binary responses from log-binomial model
    Y = rbinom(n, 1, exp(logpY1))
  )
}

df_large <- simulate_data() %>%
  glimpse()

# grid search ----
grid_search <- seq(-1, 1, 0.05) %>%
  # get the coefficient from the correct and incorrect propensity models
  map_dfr(.f = ~ {

    temp_df <- df_large %>%
      mutate(H = Y * exp(-.x * A))

    list(Correct = "A ~ H + U", Incorrect = "A ~ H") %>%
      map_dfc(.f = ~ {
        temp_df %>%
          lm(formula(.x), data = .) %>%
          tidy() %>%
          filter(term == "H") %>%
          pull(estimate)
      }) %>%
      mutate(psi = .x, .before = Correct)
  }) %>%
  print()

grid_search %>%
  pivot_longer(-psi) %>%
  ggplot(aes(x = exp(psi), y = value, color = name)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1.02, linetype = "dotdash") +
  geom_hline(yintercept = 0, linetype = "dotdash") +
  xlab("Risk ratio (true value: 1.02)") +
  ylab(
    glue::glue("Estimated coefficient for H(\U1D713) \nin ",
               "the propensity model")
  ) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  coord_cartesian(xlim = c(0, 3), ylim = c(-6, 4))

# look at the p-values
seq(-0.05, 0.05, 0.001) %>%
  # get the coefficient from the correct and incorrect propensity models
  map_dfr(.f = ~ {

    temp_df <- df_large %>%
      mutate(H = Y * exp(-.x * A))

    list(Correct = "A ~ H + U", Incorrect = "A ~ H") %>%
      map_dfc(.f = ~ {
        temp_df %>%
          lm(formula(.x), data = .) %>%
          tidy() %>%
          filter(term == "H") %>%
          pull(p.value)
      }) %>%
      mutate(psi = .x, .before = Correct)
  }) %>%
  mutate(RR = exp(psi), .after = psi) %>%
  print(n = Inf) %>%
  pivot_longer(-c(psi, RR)) %>%
  ggplot(aes(x = exp(psi), y = value, color = name)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1.02, linetype = "dotdash") +
  geom_hline(yintercept = 0, linetype = "dotdash") +
  xlab("Risk ratio (true value: 1.02)") +
  ylab("P-value that alpha1 = 0") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())

# use optimize with p-value
gest_function <- function(psi, sign = 1) {
  df_large %>%
    mutate(H = Y * exp(-psi * A)) %>%
    glm(A ~ H + U, family = gaussian(), data = .) %>%
    tidy() %>%
    filter(term == "H") %>%
    pull(p.value) %>%
    prod(., sign)
}

gest_function(0.01)
gest_function(0.01, sign = -1)

optimize(gest_function, interval = c(-5, 5), maximum = TRUE)
# $maximum
# [1] -4.985535
#
# $objective
# [1] 2.003414e-09

optim(par = 0, fn = gest_function, method = "L-BFGS-B", lower = -5, upper = 5, sign = -1)
# $par
# [1] 5.884365e-09
#
# $value
# [1] -3.440981e-13
#
# $counts
# function gradient
#        2        2
#
# $convergence
# [1] 0
#
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

optim(par = 0.01, fn = gest_function, method = "L-BFGS-B", lower = 0, upper = 0.5, sign = -1)
# $par
# [1] 0.02345147
#
# $value
# [1] -0.997653
#
# $counts
# function gradient
#       55       55
#
# $convergence
# [1] 52
#
# $message
# [1] "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"


optimize(gest_function, interval = c(0, 2), maximum = TRUE)
# $maximum
# [1] 0.02343829
#
# $objective
# [1] 0.9999589

# lm runs faster than glm
bench::mark(lm(A ~ H + U, data = mutate(df_large, H = Y * exp(-0.2 * A))), iterations = 10)
bench::mark(glm(A ~ H + U, data = mutate(df_large, H = Y * exp(-0.2 * A))), iterations = 10)

# give the same estimates and z-scores
identical(
  coef(lm(A ~ H + U, data = mutate(df_large, H = Y * exp(-0.2 * A)))),
  coef(glm(A ~ H + U, data = mutate(df_large, H = Y * exp(-0.2 * A))))
)

# g-computation + marginal effects ----
# the goal is to estimate the population CRR
# G-computation when the confounder is measured / available
glm_model <- glm(Y ~ U + A, data = df_large, family = binomial(link = "log"))

summary(glm_model)

glm_model %>%
  tidy(exponentiate = TRUE, conf.int = TRUE)

# compare with the unadjusted estimate
df_large %>%
  glm(Y ~ A, data = ., family = binomial(link = "log")) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE)

# run with smaller sample size and more simulations to asses the magnitude of the bias
1:100 %>%
  map_dbl(.x = ., .f = ~ {
    simulate_data(n = 1e4, seed = .x) %>%
      glm(Y ~ A, data = ., family = binomial(link = "log")) %>%
      tidy(exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "A") %>%
      pull(estimate)
  }) %>%
  mean() # [1] 0.9809837 instead of 1.02

# conditional risk ratio
glm_model %>%
  marginaleffects::avg_comparisons(variables = "A", comparison = "ratio")

# marginal risk ratio
glm_model %>%
  marginaleffects::avg_comparisons(variables = "A", comparison = "ratioavg")
# marginal estimate is the same as the conditional estimate due to the use of a
# collapsible link function (log instead of logit)

# Slope estimate on the response (probability) scale
# which is why the function is increasing as a function of U
glm_model %>%
  marginaleffects::plot_slopes(variables = "A", condition = "U")

# previous one by default is on the response scale,
# this one is on the link scale
glm_model %>%
  marginaleffects::plot_slopes(variables = "A", condition = "U", type = "link")

df_large %>% pull(U) %>% mean()

# P[Y = 1 | A = a] while conditioning on mean of U (= 5)
glm_model %>%
  marginaleffects::plot_predictions(condition = "A") +
  theme_bw() + ylab("Predicted probability P(Y = 1 | A = a, U = 5)")

# these two are identical
glm_model %>%
  marginaleffects::plot_predictions(condition = "A", type = "response") +
  theme_bw() + ylab("Predicted probability P(Y = 1 | A = a, U = 5)")

glm_model %>%
  marginaleffects::plot_predictions(by = "A", type = "response",
                                    newdata = expand_grid(A = 0:73, U = 4.99813)) +
  theme_bw() + ylab("Predicted probability P(Y = 1 | A = a, U = 5)")

# this differs from the conditional dose response curve above because
# it averages over the observed levels of U
# rather than interpolating over the grid defined by the ranges of A and U
glm_model %>%
  marginaleffects::plot_predictions(by = "A") +
  theme_bw() + ylab("Predicted probability of Y = 1")

glm_model %>%
  marginaleffects::plot_predictions(by = "A", type = "response") +
  theme_bw() + ylab("Predicted probability of Y = 1")

glm_model %>%
  marginaleffects::plot_predictions(by = "A", type = "response", draw = FALSE) %>%
  arrange(A) %>%
  ggplot(aes(x = A, y = estimate)) + geom_line()

# when the smooth grid is provided as newdata, the plot looks (nearly) identical to
# the one with condition = list("A"))
# the differences between this and the ones without newdata are because the one without newdata
# uses ranges of the variables from the observed dataset
glm_model %>%
  marginaleffects::plot_predictions(
    by = "A", newdata = expand_grid(A = 0:73, U = 0:10)
  ) +
  theme_bw() + ylab("Predicted probability of Y = 1")

# get predicted probabilities from the model at quartiles of U
glm_model %>%
  marginaleffects::plot_predictions(condition = list("A", "U" = "fivenum")) +
  theme_bw() + ylab("Predicted probability of Y = 1")

# G-estimation ----
# fit Gamma GEE model when U is observed

df_large %>%
  ggplot(aes(x = U, y = A)) +
  geom_point(size = 0.2, alpha = 0.2) +
  geom_smooth(method = "lm", formula = 'y ~ x') +
  geom_smooth(method = "gam", color = "orange") +
  theme_bw()

propensity_scores <- mgcv::gam(A ~ s(U, k = 10), data = df_large) %>%
  fitted()

propensity_scores %>% summary()

gest_confounder <- df_large %>%
  mutate(PS = propensity_scores) %>%
  geeM::geem(Y ~ PS + A, data = ., family = Gamma(link = "log"),
             corstr = "independence", sandwich = TRUE)

summary(gest_confounder)

tidy.geem <- function(x, conf.int = FALSE, exponentiate = FALSE, conf.level = 0.95, robust = TRUE, ...) {
  # works only when broom is attached
  # as it's added as an S3 method for broom::tidy(x, ...)
  # could also do library(broom, include.only = "tidy")
  #
  # arguments are the same as broom::tidy.geeglm()
  stopifnot(is.logical(robust))
  s <- summary(x)
  ret <- c("term" = "coefnames", "estimate" = "beta",
           "std.error" = if (robust) "se.robust" else "se.model",
           "statistic" = "wald.test", "p.value" = "p") %>%
    map(.x = ., .f = ~ pluck(s, .x)) %>%
    as_tibble()

  if (conf.int) {
    p <- (1 + conf.level) / 2
    ret <- ret %>%
      mutate(
        conf.low = estimate - (qnorm(p = p) * std.error),
        conf.high = estimate + (qnorm(p = p) * std.error),
      )
  }

  if (exponentiate) {
    # exponentiate the point estimate and the CI endpoints (if present)
    ret <- ret %>%
      mutate(
        across(
          .cols = -c(term, std.error, statistic, p.value),
          .fns = exp
        )
      )
  }

  return(ret)
}

# use the tidy method to get point estimate and CI
gest_confounder %>%
  tidy(exponentiate = TRUE, conf.int = TRUE, robust = TRUE) %>%
  filter(term == "A")

# using the estimating equation from Palmer et al
MSNMM_function <- function(psi, data = df_large, positive = TRUE) {
  result <- data %>%
    mutate(
      Zbar = mean(Z),
      s = (Z - Zbar) * Y * exp(-psi * A)
    ) %>%
    pull(s) %>%
    sum()

  if (positive) result <- abs(result)

  return(result)
}

seq(-0.05, 0.5, 0.01) %>%
  tibble(psi = .) %>%
  rowwise() %>%
  mutate(
    `Raw value` = MSNMM_function(psi, positive = FALSE),
    `Absolute value` = abs(`Raw value`)
  ) %>%
  ungroup() %>%
  pivot_longer(-psi) %>%
  mutate(name = fct_rev(factor(name))) %>%
  ggplot(aes(x = exp(psi), y = value, color = name)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0, linetype = "dotdash") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_x_continuous(n.breaks = 10) +
  xlab("Risk ratio") +
  ylab("Covariance between Z and H(\U1D713)")

optimize(f = MSNMM_function, interval = c(-2, 2), maximum = FALSE)

optimize(f = MSNMM_function, interval = c(-2, 2), maximum = FALSE, data = simulate_data(seed = 26))

# MSMM with A instead of Z
MSNMM_function_incorrect <- function(psi, data = df_large, positive = TRUE) {
  result <- data %>%
    mutate(
      Abar = mean(A),
      s = (A - Abar) * Y * exp(-psi * A)
    ) %>%
    pull(s) %>%
    sum()

  if (positive) result <- abs(result)

  return(result)
}

seq(-0.05, 0.5, 0.01) %>%
  tibble(psi = .) %>%
  rowwise() %>%
  mutate(
    `Raw value` = MSNMM_function_incorrect(psi, positive = FALSE),
    `Absolute value` = abs(`Raw value`)
  ) %>%
  ungroup() %>%
  pivot_longer(-psi) %>%
  mutate(name = fct_rev(factor(name))) %>%
  ggplot(aes(x = exp(psi), y = value, color = name)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0, linetype = "dotdash") +
  geom_vline(xintercept = 0.99, linetype = "dotdash") +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_x_continuous(n.breaks = 10) +
  xlab("Risk ratio") +
  ylab("Covariance between A and H(\U1D713)")

# average IV SMM estimates across simulations
IV_estimates <- map_dbl(
  .x = 1:100,
  .f = function(seed) {
    seed %>%
      simulate_data(seed = .) %>%
      optimize(
        f = MSNMM_function, interval = c(-2, 2), maximum = FALSE, data = .
      ) %>%
      pluck("minimum")
  })

exp(mean(IV_estimates) + (var(IV_estimates) / 2))
mean(exp(IV_estimates))
# the last two should be identical when the sample size is large
median(exp(IV_estimates))
summary(exp(IV_estimates))

# using ivtools instead
IV_model <- glm(Z ~ 1, family = binomial(link = "logit"), data = df_large)
summary(IV_model)

# should be identical
df_large %>% pull(Z) %>% mean()
plogis(coef(IV_model))

IV_gest <- df_large %>%
  as.data.frame() %>%
  ivtools::ivglm(estmethod = "g",
                 X = "A", Y = "Y",
                 fitZ.L = IV_model,
                 link = "log", data = .)

IV_gest %>% summary()

# get the coefficient
exp(summary(IV_gest)$coefficients[1, 1])
exp(confint(IV_gest))

# bootstrap to get the predicted marginal RR curves as a function of A
# plot individual curves and mean curve based on the point estimate
# plot the true RR curve on top
bootstrap_function <- function(proportion, replacement) {
  map_dfr(
    .x = 1:100,
    .f = function(seed) {
      set.seed(seed)
      boot_data <- df_large %>%
        slice_sample(prop = proportion, replace = replacement)

      est_min <- boot_data %>%
        optimize(
          f = MSNMM_function, interval = c(-2, 2), maximum = FALSE, data = .
        ) %>%
        pluck("minimum")  # estimate on the logRR scale

      Y0 <- boot_data %>%
        mutate(H = Y * exp(-est_min * A)) %>%
        pull(H) %>%
        mean()

      # dataset to return
      tibble(id = seed, estimate = est_min, Amax = max(boot_data$A), Y0 = Y0)
    })
}

boot_est <- bootstrap_function(proportion = 1, replacement = TRUE)

boot_est %>% summary()

10000 / nrow(df_large)

subsample_est <- bootstrap_function(proportion = 0.1, replacement = FALSE)

subsample_est %>% summary()

# assemble data for plotting
marginal_predictions <- bind_rows(
  boot_est %>%
    mutate(category = "Bootstrap"),
  subsample_est %>%
    mutate(category = "Subsampling (b = 10000)")
) %>%
  expand_grid(., A = 0:100) %>%
  mutate(probability = Y0 * exp(estimate * A),
         extrapolation = A > Amax) %>%
  print(n = 10)

# plot with the extrapolated predictions
#
# some extreme plots, but that can probably be improved by increasing the
# number of simulations to 5,000-10,000 and using the results between the
# 2.5% and 97.5% percentiles
marginal_predictions %>%
  ggplot(aes(x = A, y = probability, group = interaction(id, category),
             color = extrapolation)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray20") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray20") +
  theme_bw() +
  facet_wrap(~ category) +
  coord_cartesian(ylim = c(0, 1.2)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2))

# plot with only the observed values
marginal_predictions %>%
  filter(A < Amax) %>%
  mutate(`Risk ratio` = exp(estimate)) %>%
  ggplot(aes(x = A, y = probability, group = interaction(id, category),
             color = `Risk ratio`)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray20") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray20") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1.2)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2)) +
  facet_wrap(~ category)

# plot all the marginal curves with their CIs ----

# predictions from g-computation ----
df_large %>% summarize(across(.cols = c(A, U), .fns = c(min, max)))

marginaleffects::plot_predictions(
  glm_model, by = "A", type = "response", newdata = expand_grid(A = 0:73, U = 0:10)
)

marginaleffects::plot_predictions(
  glm_model, by = "A", type = "response"
)

marginaleffects::predictions(
  glm_model, by = "A", type = "response"
) %>%
  as_tibble() %>%
  arrange(A) %>%
  print(n = Inf)

df_large %>%
  mutate(pY = exp(logpY1),
         U = cut_width(U, width = 0.5, center = 0.25)) %>%
  group_by(A, U) %>%
  summarize(n = n(), p = mean(pY)) %>%
  ggplot(aes(U, A, color = cut_width(p, width = 0.05, center = 0.025), size = n)) +
  geom_point()

# compare g-computation from marginaleffects vs hand-rolled implementation
glm_model %>%
  marginaleffects::predictions(
    by = "A", type = "response",
    newdata = marginaleffects::datagrid(A = c(0, 15, 69), grid_type = "counterfactual")
  )

# gives identical estimates
map_dbl(
  .x = c(0, 15, 69),
  .f = ~ {
    df_large %>%
      mutate(A = .x) %>%
      predict(glm_model, newdata = ., type = "response") %>%
      mean()
  }
) %>%
  round(., digits = 4)

g_computation_dose_response_extrapolated_levels <- glm_model %>%
  marginaleffects::predictions(
    by = "A", type = "response",
    newdata = marginaleffects::datagrid(A = c(0, seq(5, 70, 5), 73), grid_type = "counterfactual")
  ) %>%
  select(A, estimate, conf.low, conf.high) %>%
  as_tibble() %>%
  arrange(A) %>%
  mutate(type = "G-computation") %>%
  print()

g_computation_dose_response_observed_levels <- glm_model %>%
  # marginaleffects::avg_predictions(
  #   by = "A", type = "response", newdata = expand_grid(A = 0:73, U = 0:10)
  # ) %>%
  marginaleffects::plot_predictions(
    by = "A", type = "response", draw = FALSE
  ) %>%
  select(A, estimate, conf.low, conf.high) %>%
  as_tibble() %>%
  arrange(A) %>%
  mutate(type = "G-computation (observed (A, U) levels)")

g_computation_dose_response_observed_levels

# get estimates of Pr[Y^{a = 0} = 1] ----
# use g-computation from the model adjusting for confounders
df_large %>%
  select(A, U) %>%
  mutate(A = 0) %>%
  predict(glm_model, newdata = ., type = "response") %>%
  mean()

df_large %>%
  filter(Z == 0) %>%
  summarize(mean(exp(logpY0)), mean(Y))

df_large %>%
  summarize(mean(exp(logpY0)))

# get this from the baseline log probability for the main dataset
map_dbl(.x = 1:100,
        .f = ~ {
          set.seed(.x)
          df_large %>%
            pull(logpY0) %>%
            exp() %>%
            rbinom(1e5, 1, .) %>%
            mean()
        }) %>%
  mean()

# simulate 100 datasets and take the mean of probabilities under A = 0
map_dbl(.x = 1:100,
        .f = ~ {
          simulate_data(seed = .x, n = 1e5) %>%
            pull(logpY0) %>%
            exp() %>%
            mean()
        }) %>%
  mean()

true_dose_response_curve <- tibble(
  A = 0:73,
  estimate = 0.01660574 * exp(log(1.02) * A),
  type = "True curve"
) %>%
  print(n = Inf)

# from the dataset after shaving off the effect of true A
glm_model %>%
  tidy(exponentiate = TRUE) %>%
  filter(term == "A") %>%
  pull(estimate)

df_large %>%
  mutate(H = Y * exp(-log(1.0201458) * A)) %>%
  pull(H) %>%
  mean()

# do this for 100 simulated datasets
map_dbl(.x = 1:100,
        .f = ~ {
          simulate_data(seed = .x) %>%
            mutate(H = Y * exp(-log(1.02) * A)) %>%
            pull(H) %>%
            mean()
        }) %>%
  mean()

# g-estimation curves ----
IV_gest %>% summary()
IV_gest %>% pluck("est") #%>% exp()

IV_gest %>% confint() #%>% exp()

# bootstrap percentile CIs slightly wider
boot_est %>%
  pull(estimate) %>%
  quantile(., probs = c(0.025, 0.975))

plot(density(boot_est$estimate))

g_estimation_pY0 <- c(
  IV_gest %>% pluck("est") %>% unname(),
  IV_gest %>% confint() %>% unname()
) %>%
  set_names(nm = c("estimate", "conf.low", "conf.high")) %>%
  map_dfr(
    .f = ~ {
      # get baseline P[Y = 1 | A = 0]
      df_large %>%
        mutate(H = Y * exp(-.x * A)) %>%
        pull(H) %>%
        mean() %>%
        tibble(logRR = .x, pY0 = .)
    }, .id = "name"
  ) %>%
  print()

df_large %>%
  lm(A ~ Z, data = .) %>%
  summary()

g_estimation_dose_response <- g_estimation_pY0 %>%
  expand_grid(., A = 0:73) %>%
  # get P[Y = 1 | A = a]
  mutate(pY = pY0 * exp(logRR * A)) %>%
  select(-pY0, -logRR) %>%
  pivot_wider(id_cols = A, names_from = "name", values_from = "pY") %>%
  mutate(type  = "G-estimation (using an IV)") %>%
  print()

g_estimation_dose_response

g_computation_dose_response_observed_levels

# plot the curves
bind_rows(
  true_dose_response_curve,
  g_computation_dose_response_extrapolated_levels,
  g_computation_dose_response_observed_levels,
  g_estimation_dose_response
  ) %>%
  ggplot(aes(x = A, y = estimate, color = type)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              #fill = "gray70",
              alpha = 0.1) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 0.25, 0.05)) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ylab("Predicted probability of Y = 1")

title_plot <- bind_rows(
  g_computation_dose_response_extrapolated_levels,
  g_estimation_dose_response
) %>%
  ggplot(aes(x = A, y = estimate, color = type)) +
  geom_line() +
  # plot the true curve
  geom_line(
    data = true_dose_response_curve,
    aes(x = A, y = estimate),
    color = "gray20", linewidth = 1.2,, linetype = "dotdash",
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = type), alpha = 0.2
  ) +
  scale_y_continuous(breaks = seq(0, 0.25, 0.05)) +
  scale_x_continuous(breaks = seq(0, 80, 10)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ylab("Predicted probability of Y = 1") +
  xlab("Number of impressions (true curve in gray)")

title_plot

ggsave(plot = title_plot, filename = "image-MSMM-IV.png",
       device = "png", units = "px", width = 900, height = 800, dpi = 180,
       path = fs::path("posts", "g-estimation-of-MSMM-with-an-IV"))
