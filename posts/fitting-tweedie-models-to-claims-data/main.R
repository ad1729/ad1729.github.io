
# Load packages and read data -----------------

post_dir <- fs::path("posts", "fitting-tweedie-models-to-claims-data")

get_path <- function(fname) fs::path(post_dir, fname)

library(tidyverse)
library(magrittr, include.only = "%$%") # importing the exposition pipe
library(tweedie)

theme_set(theme_bw())

claims_full <- read_rds(get_path("claims.rds")) %>%
  glimpse()

claims <- read_csv(get_path("claims_subset.csv"),
                   show_col_types = FALSE) %>%
  # select(IDpol, Exposure, ClaimAmount) %>%
  glimpse()

totals <- claims %>% select(Exposure, ClaimAmount) %>% map(.f = sum) %>% print()

expected_claim_amount <- round(totals$ClaimAmount / totals$Exposure, 2) %>% print()

expected_claim_amount * 60000

# number of levels in each variable
claims %>%
  summarise(across(.cols = everything(), .fns = ~ length(unique(.x)))) %>%
  glimpse()

# Bootstrapping the weighted mean -----------------

expected_claim_cost_fun <- function(data, indx, ...) {
  data <- data[indx, ]

  expected_value <- sum(data$ClaimAmount) / sum(data$Exposure)

  # this is the single / largest outlier in the data
  outlier_counts <- nrow(data[data$ClaimAmount == 1404185.52, ])

  return(c(expected_value, outlier_counts))
}

# test run with small number of resamples
boot_fun <- function(data, R = 100, parallel = "snow") {
  stopifnot(parallel %in% c("no", "snow"))

  # TRUE if using parallelization, otherwise FALSE
  simple <- parallel == "snow"

  boot::boot(
    data = data,
    statistic = expected_claim_cost_fun,
    R = R,
    sim = "ordinary",
    stype = "i",
    simple = simple,
    parallel = parallel,
    ncpus = 18
  )
}

# boot_fun(data = claims, R = 100, parallel = "no")
# boot_fun(data = claims, R = 100, parallel = "snow")

# full bootstrap run with 10k resamples
# boot_fit <- boot_fun(data = claims, R = 10000)
#
# boot_fit
#
# save the results to disk
# saveRDS(boot_fit, file = get_path("bootstrap_expected_claim_cost.rds"))

# read the bootstrapped results back into R
boot_fit <- readRDS(get_path("bootstrap_expected_claim_cost.rds"))

# compute intervals, these are quick to compute
boot::boot.ci(boot_fit, type = c("norm", "basic", "perc"))

# convert the results into a data frame
boot_dist <- tibble(
  expected_claim_cost = boot_fit$t[, 1],
  outlier_counts = boot_fit$t[, 2],
  # used for coloring / facetting plots
  `Outlier counts` = paste0(boot_fit$t[, 2], " replicates")
)

glimpse(boot_dist)

# Plotting the multimodal bootstrap distribution -----------------

sampling_dist_plot <- boot_dist %>%
  ggplot(aes(x = expected_claim_cost)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = expected_claim_amount, color = "orange",
             linewidth = 1.2, linetype = "dashed") +
  xlab("Bootstrap distribution of expected exposure-adjusted claim amount per person per year") +
  ylab("Count")

sampling_dist_plot

sampling_dist_plot_colored <- boot_dist %>%
  ggplot(aes(x = expected_claim_cost,
             fill = `Outlier counts`)) +
  geom_histogram(bins = 100) +
  xlab("Bootstrap distribution of expected exposure-adjusted claim amount per person per year") +
  ylab("Count") +
  theme(legend.position = "bottom")

sampling_dist_plot_colored +
  geom_vline(xintercept = expected_claim_amount, color = "gray40",
             linewidth = 1.2, linetype = "dashed")

sampling_dist_plot_colored +
  facet_wrap(~ `Outlier counts`, scales = "free") +
  theme(legend.position = "none")

nonparametric_bootstrap_totals <- boot_dist %>%
  mutate(
    total_claim_amount = 60000 * expected_claim_cost,
    total_claim_amount_in_millions = total_claim_amount / 1e6,
    method = "Nonparametric Bootstrap (Weighted mean)"
  ) %>%
  select(method, total_claim_amount_in_millions) %>%
  print(n = 10)

nonparametric_bootstrap_totals %>%
  pull(total_claim_amount_in_millions) %>%
  quantile(probs = sort(c(seq(0, 1, 0.25), 0.95, 0.99, 0.999))) %>%
  round(., 2)

next_year_total_plot <- nonparametric_bootstrap_totals %>%
  ggplot() +
  stat_ecdf(aes(x = total_claim_amount_in_millions), pad = FALSE) +
  geom_vline(xintercept = (6e4 * expected_claim_amount) / 1e6,
             color = "orange", linewidth = 1.2, linetype = "dashed") +
  scale_x_continuous(breaks = seq(8, 30, 2)) +
  scale_y_continuous(labels = scales::label_percent()) +
  xlab("Plausible values for the following year's total claim amount assuming unit exposure (in millions)") +
  ylab("Empirical distribution function")

next_year_total_plot +
  theme(legend.position = "none")

# Tweedie GLM with intercept -----------------

# switch to tweedie models for using covariates to predict
# fit an intercept only tweedie model
#
# for regression models, either model offsets, or
# exposure weights (modelling of pure premium = ClaimAmount / Exposure) can be used
# for a tweedie distribution, these lead to different estimates
# (for poisson they give identical estimates)#
weighted_mean_claim_amount <- weighted.mean(
  x = claims$ClaimAmount / claims$Exposure,
  w = claims$Exposure
) %>% print()

sum(claims$ClaimAmount) / sum(claims$Exposure)

# using link.power = 1 implies identity link for the mean parameter mu
tweedie_intercept_only <- glm(
  I(ClaimAmount / Exposure) ~ 1,
  weights = Exposure,
  data = claims,
  family = statmod::tweedie(var.power = 1.6, link.power = 1)
)

summary(tweedie_intercept_only)

# using link.power = 0 implies log link for the mean parameter mu
summary(glm(I(ClaimAmount / Exposure) ~ 1, weights = Exposure, data = claims,
            family = statmod::tweedie(var.power = 1.6, link.power = 0)))

# log link
# tweedie::tweedie.profile(I(ClaimAmount / Exposure) ~ 1, weights = Exposure, data = claims, link.power = 0, p.vec = seq(1.2, 1.7, by = 0.1)) %>%
#   summary()
# 1.2 1.3 1.4 1.5 1.6 1.7
# ......Done.
#   No valid values of the likelihood computed: smooth aborted
#    Consider trying another value for the input  method.
# Error in if ((xi.max > 0) & (xi.max < 1)) { :
#   missing value where TRUE/FALSE needed
# In addition: Warning messages:
# 1: In model.matrix.default(mt, mf, contrasts) :
#   non-list contrasts argument ignored
# 2: In max(y, na.rm = TRUE) :
#   no non-missing arguments to max; returning -Inf
# prof <- claims %>%
#   mutate(pp = ClaimAmount / Exposure) %>%
#   tweedie::tweedie.profile(pp ~ 1, weight = Exposure, data = ., link.power = 0,
#                            p.vec = seq(1.4, 1.63, length.out = 8),
#                            do.smooth = TRUE, method = "series", verbose = 2)
# print(prof)
#
# no matter what values of p or method I try, the log-likelihood is -Inf
#
# switching to the exposure offset formulation makes the function run without errors
# prof <- claims %>%
#   mutate(pp = ClaimAmount / Exposure) %>%
#   tweedie::tweedie.profile(
#   ClaimAmount ~ 1, offset = log(Exposure), data = ., link.power = 0,
#   #p.vec = seq(1.35, 1.65, by = 0.05),
#   p.vec = c(1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7),
#   do.smooth = TRUE, method = "series", verbose = 2)
# #plot(prof)
# print(prof) # best p is about 1.630612
#
# use alternate parameterization to cast the weighted pure premium
# model into an unweighted with exposure as offset to see if
# tweedie.profile() doesn't give infinite LL estimates
# claims %>%
#   mutate(
#     response = ClaimAmount * (Exposure ^ ((1.6 - 1) / (2 - 1.6))),
#     off = Exposure ^ (1 / (2 - 1.6))
#   ) %>%
#   select(response, off) %>%
#   summary()
#
# temp_mod <- claims %>%
#   mutate(
#     response = ClaimAmount * (Exposure ^ ((1.6 - 1) / (2 - 1.6))),
#     off = Exposure ^ (1 / (2 - 1.6))
#   ) %>%
#   select(response, off) %>%
#   glm(response ~ 1, offset = log(off), data = .,
#       family = statmod::tweedie(var.power = 1.6, link.power = 0))
#
# temp_mod
# temp_mod %>% coef() %>% exp()
#
# summary(temp_mod)
#
# tweedie::logLiktweedie(temp_mod, dispersion = summary(temp_mod)$dispersion)
# tweedie::AICtweedie(temp_mod, dispersion = summary(temp_mod)$dispersion)
# tweedie::AICtweedie(tweedie_intercept_only,
#                     dispersion = summary(temp_mod)$dispersion)
#
# ll_vs_p <- map_dfr(
#   .x = seq(1.2, 1.8, length.out = 20),
#   .f = ~ {
#     p <- .x
#
#     temp_mod <- glm(
#       I(ClaimAmount / Exposure) ~ 1, weights = Exposure, data = claims,
#       family = statmod::tweedie(var.power = p, link.power = 0)
#     )
#     temp_mu <- exp(coef(temp_mod))
#     temp_phi <- summary(temp_mod)$dispersion
#
#     tibble(
#       p = p, mu = temp_mu, phi = temp_phi,
#       ll = tweedie::logLiktweedie(temp_mod, dispersion = temp_phi),
#       AIC = tweedie::AICtweedie(temp_mod, dispersion = temp_phi),
#       converged = temp_mod$converged
#     )
#   })
#
# ll_vs_p
#
# prof <- claims %>%
#   mutate(
#     response = ClaimAmount
#   ) %>%
#   tweedie::tweedie.profile(
#     response ~ 1, offset = off, data = ., link.power = 0,
#     p.vec = c(1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7),
#     do.smooth = TRUE, method = "series", verbose = 2)
# print(prof) # best p is about 1.630612
#
# those parameters above are for the population with unit exposure
# so claims for next year can be sampled as Y_i ~ rtweedie(60000, mu_i, phi_i, p = 1.6)
# doing the same with varying exposures (to compare actuals vs total for example)
# the phi_i gets replaced with phi_i / w_i
#
# here we're assuming constant mu and phi, so mu_i = mu, and phi_i = phi for all i
#
# to compare, we can sample Z_i = w_i * Y_i, where Y_i ~ tweedie(mu, phi / w_i, p)
#
# quick sanity check on simulated data
set.seed(43)
tibble(w = runif(2e5, min = 0.1, max = 1)) %>%
  # sample pure premium values
  mutate(y = tweedie::rtweedie(n(), mu = 200, phi = 350 / w, power = 1.3)) %>%
  # model pure premium
  glm(y ~ 1, weights = w, data = ., family = statmod::tweedie(var.power = 1.3, link.power = 1)) %>%
  summary()

set.seed(43)
tibble(w = runif(2e5, min = 0.1, max = 1)) %>%
  mutate(
    # sample latent pure premium values under exposure w_i
    y = tweedie::rtweedie(n(), mu = 200, phi = 350 / w, power = 1.3),
    # create obs claim amount z for exposure w
    z = w * y
  ) %>%
  glm(I(z / w) ~ 1, weights = w, data = .,
      family = statmod::tweedie(var.power = 1.3, link.power = 1)) %>%
  summary()

# draw values from tweedie with exposures and compute the claim totals across different resamples
mu <- coef(tweedie_intercept_only)
mu %>% as.character()

phi <- summary(tweedie_intercept_only) %>% pluck("dispersion")
phi %>% as.character()

# Posterior predictive checking for observed exposure -----------------

# how well does the model fit the data?
# a form of posterior predictive checking
# ppc_obs_exposure <- map(
#   .x = 1:10000,
#   .f = ~ {
#     if(.x %% 100 == 0) {
#       print(.x)
#     }
#     set.seed(.x)
#     draw <- claims$Exposure * tweedie::rtweedie(60000, mu = mu,
#                                                 phi = phi / claims$Exposure,
#                                                 power = 1.6)
#     tibble(prop_zero = mean(draw == 0), sample_total = sum(draw),
#            sample_max = max(draw), n_nonzero = sum(draw > 0))
#   }) %>%
#   list_rbind()
#
# saveRDS(ppc_obs_exposure, file = get_path("ppc_obs_exposure.rds"))

ppc_obs_exposure <- readRDS(file = get_path("ppc_obs_exposure.rds"))

sample_statistics_obs_exposure <- claims %>%
  summarise(
    prop_zero = mean(ClaimAmount == 0),
    sample_total = sum(ClaimAmount),
    sample_max = max(ClaimAmount),
    n_nonzero = sum(ClaimAmount > 0)
  )

ppc_obs_exposure %>% summary()
sample_statistics_obs_exposure

# can make this into a tidy function of sorts to clean up summary output
tidy_df_summary <- function(data) {
  data %>%
    summary() %>%
    as.data.frame() %>%
    separate(col = Freq, into = c("Statistic", "Value"), sep = ":") %>%
    select(-Var1) %>%
    rename(Column = Var2) %>%
    pivot_wider(names_from = Statistic, values_from = Value) %>%
    mutate(
      Column = str_trim(Column, side = "left"),
      across(
        .cols = where(is.character),
        .fns = ~ str_trim(.x, side = "both")
      )
    )
}

ppc_obs_exposure %>%
  tidy_df_summary()

plot_data_obs_exposure <- bind_rows(
  ppc_obs_exposure %>% mutate(group = "sampled", .before = 0),
  sample_statistics_obs_exposure %>% mutate(group = "observed", .before = 0)
) %>%
  mutate(
    across(c(sample_total, sample_max), ~ .x / 1e6),
    prop_zero = 100 * prop_zero
  ) %>%
  rename(
    `% of policies with zero claims` = prop_zero,
    `Number of policies with non zero claim amounts` = n_nonzero,
    `Maximum claim amount (in millions)` = sample_max,
    `Total claim amount (in millions)` = sample_total
  ) %>%
  pivot_longer(cols = -group, names_to = "statistic", values_to = "values")

# compare these visually
ppc_mean_obs_exposure <- plot_data_obs_exposure %>%
  filter(group == "sampled") %>%
  summarise(values = mean(values), .by = statistic) %>%
  print(n = Inf)

plot_data_obs_exposure %>%
  filter(group == "sampled") %>%
  ggplot(aes(x = values, group = statistic)) +
  stat_ecdf(pad = FALSE) +
  # plot the sample statistic
  geom_vline(data = filter(plot_data_obs_exposure,
                           group == "observed"),
             aes(xintercept = values, group = statistic),
             color = "orange", linewidth = 1.2, linetype = "dashed") +
  # plot the distribution means
  geom_vline(data = ppc_mean_obs_exposure,
             aes(xintercept = values, group = statistic),
             color = "darkred", linewidth = 1.2, linetype = "dashed") +
  facet_wrap(~ statistic, scales = "free") +
  scale_y_continuous(labels = scales::label_percent()) +
  xlab(glue::glue("Sample statistics for the observed exposures (in orange);",
                  "\nmean of the sampling distributions (in red)")) +
  ylab("Empirical distribution function")

sample_statistics_obs_exposure
ppc_mean_obs_exposure

# plot eCDFs of sampled datasets with the observed data
map(
  .x = 1:10,
  .f = ~ {
    set.seed(.x)
    draw <- claims$Exposure * tweedie::rtweedie(60000, mu = mu,
                                                phi = phi / claims$Exposure,
                                                power = 1.6)
    tibble(sim_id = .x, y = draw, grp = "Simulated")
  }) %>%
  list_rbind() %>%
  bind_rows(
    .,
    claims %>%
      mutate(sim_id = 100, grp = "True") %>%
      select(sim_id, y = ClaimAmount, grp)
  ) %>%
  filter(y > 0) %>%
  ggplot(aes(x = y, group = sim_id, color = grp)) +
  stat_ecdf() +
  scale_color_manual(values = c("Simulated" = "gray70", "True" = "Black")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks() +
  xlab("Policies with non-zero claim amounts (log10 scale)") +
  ylab("Empirical distribution function") +
  theme(legend.title = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.9, 0.2))

# Parametric bootstrap for next year's totals -----------------

# # do the previous for 10k steps assuming unit exposure,
# # and compare with the totals
# # here we can do it for a smaller rep B
# predicted_totals_for_unit_exposure <- map_dbl(.x = 1:10000, .f = ~ {
#   if(.x %% 100 == 0) {
#     print(.x)
#   }
#   set.seed(.x)
#   sum(tweedie::rtweedie(60000, mu = mu, phi = phi, power = 1.6))
# })
#
# saveRDS(predicted_totals_for_unit_exposure,
#         file = get_path("predicted_totals_for_unit_exposure.rds"))

predicted_totals_for_unit_exposure <- readRDS(file = get_path("predicted_totals_for_unit_exposure.rds"))

summary(predicted_totals_for_unit_exposure)

next_year_total_plot +
  stat_ecdf(aes(x = predicted_totals_for_unit_exposure / 1e6),
            linetype = "dotdash", pad = FALSE) +
  scale_x_continuous(breaks = seq(3, 30, 2))

# Nonparametric bootstrap for joint distr. of (mu, phi) -----------------

# bootstrap and get the joint distribution of (mu, phi)
# bootstrap_mu_phi <- map_dfr(
#   .x = 1:10000,
#   .f = ~ {
#     if(.x %% 100 == 0) {
#       print(.x)
#     }
#     set.seed(.x)
#     data <- claims %>%
#       slice_sample(n = 60000, replace = TRUE)
#     outlier_counts <- nrow(data[data$ClaimAmount == 1404185.52, ])
#     mod <- data %>%
#       glm(
#         I(ClaimAmount / Exposure) ~ 1,
#         weights = Exposure,
#         data = .,
#         family = statmod::tweedie(var.power = 1.6, link.power = 0)
#       )
#
#     tibble(
#       mu = exp(coef(mod)),
#       phi = summary(mod)$dispersion,
#       outlier_counts = outlier_counts
#     )
#   }
# )
#
# saveRDS(bootstrap_mu_phi,
#         file = get_path("bootstrap_mu_phi.rds"))

bootstrap_mu_phi <- readRDS(get_path("bootstrap_mu_phi.rds"))

summary(bootstrap_mu_phi)

bootstrap_mu_phi

# how many of the (mu, dispersion) pairs are unique upto 0 significant digits?
bootstrap_mu_phi %>% nrow()
bootstrap_mu_phi %>% map(.f = ~ length(unique(.x)))
bootstrap_mu_phi %>% map(.f = ~ length(unique(round(.x))))

mu_phi_plot <- bootstrap_mu_phi %>%
  mutate(`Outlier counts` = factor(outlier_counts)) %>%
  ggplot(aes(x = mu, y = phi, color = `Outlier counts`)) +
  geom_point() +
  geom_point(data = tibble(mu = mu, phi = phi),
             aes(x = mu, y = phi),
             color = "gray20", size = 5, inherit.aes = FALSE) +
  labs(x = "Mean (\u03bc)", y = "Dispersion (\u03d5)") +
  theme(
    legend.position = "inside",
    legend.background = element_blank(),
    legend.position.inside = c(0.8, 0.18)
  ) +
  guides(color = guide_legend(nrow = 2))

mu_phi_plot

# ggExtra::ggMarginal(mu_phi_plot, type = "density")
suppressWarnings(ggExtra::ggMarginal(mu_phi_plot, type = "densigram", xparams = list(bins = 100), yparams = list(bins = 100)))

ggsave(plot = mu_phi_plot,
       filename = "tweedie_bootstrap_plot.png",
       device = "png", units = "px", width = 900, height = 800, dpi = 160,
       path = fs::path("posts", "fitting-tweedie-models-to-claims-data"))

# compare with running a bootstrap on tweedie simulated data
# set.seed(43)
# foo <- tweedie::rtweedie(n = 60000, mu = 250, phi = 3000, p = 1.6) %>% round(., 2)
# #foo %>% sort() %>% tail(n = 30)
# foo_boot <- map_dfr(
#   .x = 1:500,
#   .f = ~ {
#     if(.x %% 50 == 0) {
#       print(.x)
#     }
#     set.seed(.x)
#     samp <- foo %>%
#       sample(x = ., size = 60000, replace = TRUE)
#
#     mod <- glm(samp ~ 1, family = statmod::tweedie(var.power = 1.6, link.power = 0))
#
#     tibble(mu = exp(coef(mod)), phi = summary(mod)$dispersion)
#   }
# )
#
# foo_boot %>%
#   ggplot(aes(x = mu, y = phi)) +
#   geom_point() +
#   geom_point(inherit.aes = FALSE,
#              data = tibble(mu = 250, phi = 3000),
#              aes(x = mu, y = phi), color = "orange", size = 5)


# Use posterior distribution for (mu, phi) for next year's totals -----------------

# sample Y_i from the posterior predictive distribution and sum them
# do this 10k times to get the posterior distribution
# posterior_distribution_samples_for_total_claims <- map(
#   .x = 1:10000,
#   .f = ~ {
#     if(.x %% 100 == 0) {
#       print(.x)
#     }
#     set.seed(.x)
#     draw <- bootstrap_mu_phi %>%
#       slice_sample(n = 1)
#
#     set.seed(.x)
#     sum_values <- draw %$%
#       tweedie::rtweedie(n = 60000, mu = mu, phi = phi, power = 1.6) %>%
#       sum()
#
#     draw %>%
#       mutate(total = sum_values)
#   }
# ) %>%
#   list_rbind()
#
# saveRDS(posterior_distribution_samples_for_total_claims,
#         file = get_path("posterior_distribution_samples_for_total_claims.rds"))

posterior_distribution_samples_for_total_claims <- readRDS(get_path("posterior_distribution_samples_for_total_claims.rds"))

posterior_distribution_samples_for_total_claims %>% summary()

posterior_distribution_samples_for_total_claims %>%
  mutate(
    totals_in_millions = total / 1e6,
    totals_bins = cut_width(totals_in_millions, width = 3)
  ) %>%
  ggplot(aes(x = mu, y = phi, color = totals_bins)) +
  geom_point() +
  geom_point(data = tibble(mu = mu, phi = phi),
             aes(x = mu, y = phi),
             color = "orange", size = 10, inherit.aes = FALSE) +
  labs(x = "Mu", y = "Dispersion (Phi)")

# Use closed-form tweedie distribution for next year's totals -----------------

# distribution of means with total exposure w = sum(w_i)
# is distributed as T_i ~ Tweedie(mu_i, phi / w, p = 1.6)
# add eCDF for 10k samples from this to the plot as well
#
# this is the reproductive property from this page:
# https://en.wikipedia.org/wiki/Tweedie_distribution#Reproductive_exponential_dispersion_models
#
# use the scale invariance property too: https://en.wikipedia.org/wiki/Tweedie_distribution#Scale_invariance
set.seed(43)
predicted_totals_tweedie_sampling_dist <- tweedie::rtweedie(
  n = 10000,
  mu = 60000 * mu,
  phi = ((60000 ^ (2 - 1.6)) * phi) / totals$Exposure,
  power = 1.6
)

predicted_totals_tweedie_sampling_dist %>% summary()

# 95% CIs for the mean -----
broom::tidy(tweedie_intercept_only, conf.int = TRUE, exponentiate = FALSE)
boot::boot.ci(boot_fit, type = c("norm", "basic", "perc"))
# using the boot package for BCa intervals doesn't
# work because of the very large sample sizes
# so we can sue the bcaboot package
# bca_res <- bcaboot::bcajack(
#   x = claims,
#   B = 500,
#   verbose = TRUE,
#   func = function(data) {
#     sum(data$ClaimAmount) / sum(data$Exposure)
#   },
#   alpha = 0.025
# )
#
# saveRDS(bca_res, file = get_path("bca_confint.rds"))

bca_res <- readRDS(file = get_path("bca_confint.rds"))
bca_res

bca_res$lims[c(1,3), 1]

tweedie_CI <- tweedie::rtweedie(
  n = 10000, mu = mu, phi = phi / totals$Exposure, power = 1.6
) %>%
  quantile(probs = c(0.025, 0.975)) %>%
  print()

param_boot_CI <- quantile(predicted_totals_for_unit_exposure / 60000, probs = c(0.025, 0.975))

np_boot_glm_CI <- posterior_distribution_samples_for_total_claims %>%
  mutate(mean = total / 60000) %>%
  pull(mean) %>%
  quantile(probs = c(0.025, 0.975)) %>%
  print()

# some skewness in the sampling distribution (as expected)
hist(predicted_totals_tweedie_sampling_dist / 60000, breaks = 100)

bind_rows(
  broom::tidy(tweedie_intercept_only, conf.int = TRUE, exponentiate = FALSE) %>%
    mutate(method = "Model based") %>%
    select(method, conf.low, conf.high),
  broom::tidy(boot_fit, conf.method = "perc", conf.int = TRUE) %>%
    slice_head(n = 1) %>%
    mutate(method = "NP bootstrap (Percentile)") %>%
    select(method, conf.low, conf.high),
  tibble(method = "NP bootstrap (BCa)",
         conf.low = bca_res$lims[c(1,3), 1][1],
         conf.high = bca_res$lims[c(1,3), 1][2]),
  tibble(method = "Tweedie sampling distr.",
         conf.low = tweedie_CI[1],
         conf.high = tweedie_CI[2]),
  tibble(method = "Parametric bootstrap",
         conf.low = param_boot_CI[1],
         conf.high = param_boot_CI[2]),
  tibble(method = "PPC mean",
         conf.low = np_boot_glm_CI[1],
         conf.high = np_boot_glm_CI[2]),
) %>%
  arrange(desc(conf.high))

# Fit a GLM with a single covariate -----------------
summary(glm(I(ClaimAmount / Exposure) ~ 1 + factor(VehGas),
            weights = Exposure, data = claims,
            family = statmod::tweedie(var.power = 1.6, link.power = 1)))

glm(I(ClaimAmount / Exposure) ~ 1 + factor(VehGas),
    weights = Exposure, data = claims,
    family = statmod::tweedie(var.power = 1.6, link.power = 1)) %>%
  coef() %>% cumsum()

glm(I(ClaimAmount / Exposure) ~ 1 + factor(VehGas),
    weights = Exposure, data = claims,
    family = statmod::tweedie(var.power = 1.6, link.power = 0)) %>%
  coef() %>% cumsum() %>% exp()

claims %>%
  summarise(
    ratio = as.character(round(sum(ClaimAmount) / sum(Exposure), 4)),
    .by = VehGas
  )

# fit main effects model with main variables
claims %>% glimpse()

# number of distinct values for each variable
claims %>%
  summarise(across(.cols = everything(), .fns = ~ length(unique(.x)))) %>%
  glimpse()

# which variables are categorical
claims %>%
  mutate(
    across(c(Area, VehPower, VehBrand, VehGas, Region), ~ as.factor(.x))
  ) %>%
  summary()

# Fit a GLM with main effects -----------------

claims_modelling <-  claims %>%
  mutate(
    across(c(Area, VehPower, VehBrand, VehGas, Region), ~ as.factor(.x)),
    pure_premium = ClaimAmount / Exposure
  )

first_model <- glm(
  pure_premium ~ 1 + Area + VehPower + VehAge + DrivAge
  + BonusMalus + VehBrand + VehGas + Density + Region,
  weights = Exposure,
  data = claims_modelling,
  family = statmod::tweedie(var.power = 1.6, link.power = 0)
)

# get dispersion parameter estimates
first_model %>% summary()

first_model %>%
  broom::tidy(exponentiate = TRUE) %>%
  arrange(p.value)

# 52 parameters + 1 intercept
first_model %>% pluck(coef) %>% length() %>% {. - 1}

fitted_values <- first_model %>%
  broom::augment(newdata = claims_modelling, type.predict = "response") %>%
  glimpse()

fitted_values %>% pull(.fitted) %>% summary()

mean_fitted <- fitted_values %>%
  pull(".fitted") %>%
  mean() %>%
  print()

fitted_values %>%
  pull(.fitted) %>%
  quantile(., probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1.0)) %>%
  c(., "Mean" = mean_fitted) %>%
  sort() %>%
  round(., 2)

fitted_values %>%
  ggplot(aes(x = .fitted)) +
  stat_ecdf() +
  geom_vline(xintercept = mu, color = "orange", linetype = "dashed") +
  annotate(geom = "text", x = 200, y = 0.8, color = "orange",
           label = "Overall mean: 217.8", hjust = "right") +
  geom_vline(xintercept = mean_fitted, color = "red4", linetype = "dashed") +
  annotate(geom = "text", x = 275, y = 0.5, color = "red4",
           label = "Fitted values mean: 256.3", hjust = "left") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks() +
  scale_y_continuous(labels = scales::percent) +
  xlab("Fitted pure premium values") +
  ylab("Empirical distribution function")

# fitted_values %>%
#   mutate(
#     y_act = log(pure_premium + 1),
#     y_pred = log(.fitted)
#   ) %>%
#   ggplot(aes(x = y_pred, y = y_act)) +
#   geom_point()

fitted_values %>%
  arrange(desc(.fitted)) %>%
  print(n = 20)

fitted_values %>%
  arrange(.fitted) %>%
  print(n = 20)

first_model %>% pluck("df.residual")
first_model %>% pluck("df.residual") %>% {. + 53}

# Pearson estimator of dispersion
# see the previous blog post for the formula
fitted_values %>%
  select(pure_premium, fitted = .fitted, weights = Exposure) %>%
  mutate(
    scaled_res = (weights * ((pure_premium - fitted) ^ 2)) / (fitted) ^ 1.6
  ) %>%
  summarise(
    phi = (1 / pluck(first_model, "df.residual")) * sum(scaled_res)
  )

phi_glm <- first_model %>% summary() %>% pluck("dispersion")
phi_glm

# what feature(s) are contributing to the largest k pure premiums
fitted_values %>%
  arrange(desc(.fitted)) %>%
  slice_head(n = 10) %>%
  glimpse()

# Which predictors contribute to the largest predictions? ----

# get the term-wise contribution to the prediction on the
# link scale for the top k largest predictions
top_10_largest_policy_predictions <- fitted_values %>%
  arrange(desc(.fitted)) %>%
  slice_head(n = 10) %>%
  predict(first_model, newdata = ., type = "terms")

attr(top_10_largest_policy_predictions, "constant")

# convert the 'terms' data frame to tibble and add constant value
# then pivot to get the exp(cumsum()) for each policy (i.e. row)
# then pivot back
# there's probably some neat function that does
# this in fewer lines of code
top_10_largest_policy_predictions %>%
  as_tibble() %>%
  mutate(
    id = row_number(),
    constant = attr(top_10_largest_policy_predictions, "constant"),
    .before = 0
  ) %>%
  pivot_longer(cols = -id) %>%
  group_by(id) %>%
  mutate(value2 = cumsum(value),
         value3 = exp(value2)) %>%
  ungroup() %>%
  print(n = 12) %>%
  pivot_wider(id_cols = id, names_from = name, values_from = value3) %>%
  select(-id) %>%
  arrange(desc(Region)) %>%
  print()
  #glimpse() %>%
  #print.data.frame()

claims_modelling %>% count(BonusMalus) %>% arrange(desc(BonusMalus))

first_model %>%
  broom::tidy() %>%
  filter(term == "BonusMalus") # est of 0.0430 on log scale

first_model %>% update(. ~ BonusMalus) # est of 0.0395 on log scale in univ. model

# Predicted totals from GLM main effects -----
# simulate 10,000 times the total claim amounts for the full year of exposure
# using the fitted means
# predicted_totals_from_glm <- map_dbl(.x = 1:10000, .f = ~ {
#   if(.x %% 100 == 0) {
#     print(.x)
#   }
#   set.seed(.x)
#   sum(tweedie::rtweedie(60000, mu = fitted_values$.fitted, phi = phi_glm, power = 1.6))
# })
#
# saveRDS(predicted_totals_from_glm,
#         file = get_path("predicted_totals_from_glm.rds"))

predicted_totals_from_glm <- readRDS(file = get_path("predicted_totals_from_glm.rds"))

summary(predicted_totals_from_glm)
# why is the mean of the conditional means
#  higher here than the one from the intercept only model?
#  The law of iterated expectation indicates these should be the same
#
#  test this on simulated data
set.seed(43)
foo <- tibble(
  x = sample(1:1000, size = 500000, replace = TRUE),
  y = tweedie::rtweedie(n = 500000, p = 1.6, mu = 250 + (4 * x), phi = 1000)
)
foo
summary(foo)
mean(foo$y)
glm(y ~ x,
    data = foo,
    family = statmod::tweedie(var.power = 1.6, link.power = (1 - 1.6))) %>%
  broom::augment(type.predict = "response") %>%
  pull(.fitted) %>%
  mean()
glm(y ~ x,
    data = foo,
    family = statmod::tweedie(var.power = 1.6, link.power = 0)) %>%
  broom::augment(type.predict = "response") %>%
  pull(.fitted) %>%
  mean()
glm(y ~ x,
    data = foo,
    family = statmod::tweedie(var.power = 1.6, link.power = 1)) %>%
  broom::augment(type.predict = "response") %>%
  pull(.fitted) %>%
  mean()
# so the law of iterated expectation (LIE) here holds only for
# the canonical link function
# where E[Y] matches the E[E[Y|X]] for the canonical link
#
# see brady 2019 slides - slide 13
# see Denuit et al 2021 - Autocalibration and Tweedie-dominance for insurance pricing with machine learning
# can add a calibration step before sampling

# try to fit the model with the canonical link
first_model %>%
  update(
    # the default link function is the canonical link function
    family = statmod::tweedie(var.power = 1.6)
  )

# PPC for obs exposure from main effects model -----
# quick check of sampling dist of statistics
# ppc_glm_obs_exposure <- map(
#   .x = 1:1000,
#   .f = ~ {
#     if(.x %% 50 == 0) {
#       print(.x)
#     }
#     set.seed(.x)
#     draw <- (fitted_values$Exposure *
#                tweedie::rtweedie(60000,
#                                  mu = fitted_values$.fitted,
#                                  phi = phi_glm, power = 1.6))
#
#     tibble(prop_zero = mean(draw == 0), sample_total = sum(draw),
#            sample_max = max(draw), n_nonzero = sum(draw > 0))
#   }) %>%
#   list_rbind()
#
# saveRDS(ppc_glm_obs_exposure, file = get_path("ppc_glm_obs_exposure.rds"))

ppc_glm_obs_exposure <- readRDS(file = get_path("ppc_glm_obs_exposure.rds"))

ppc_glm_obs_exposure %>% summary()

# compare against the sample stats
sample_statistics_obs_exposure

# compare with the intercept only GLM
ppc_obs_exposure %>% summary()

# combine the three summaries into a single data frame
list(
  "Intercept" = ppc_obs_exposure,
  "Main effects" = ppc_glm_obs_exposure,
  "Observed" = sample_statistics_obs_exposure
) %>%
  imap(
    .f = ~ {
      .x %>%
        tidy_df_summary() %>%
        mutate(Method = .y, .before = 1)
    }
  ) %>%
  list_rbind() %>%
  mutate(
    Method = factor(Method, levels = c("Observed", "Intercept", "Main effects"))
  ) %>%
  arrange(Column, Method) %>%
  print(n = Inf)

# eCDF plots
map(
  .x = 1:10,
  .f = ~ {
    set.seed(.x)
    draw <- (fitted_values$Exposure *
               tweedie::rtweedie(60000,
                                 mu = fitted_values$.fitted,
                                 phi = phi_glm, power = 1.6))
    tibble(sim_id = .x, y = draw, grp = "Simulated")
  }) %>%
  list_rbind() %>%
  bind_rows(
    .,
    claims %>%
      mutate(sim_id = 100, grp = "Observed") %>%
      select(sim_id, y = ClaimAmount, grp)
  ) %>%
  filter(y > 0) %>%
  ggplot(aes(x = y, group = sim_id, color = grp)) +
  stat_ecdf() +
  xlab("Policies with non-zero claim amounts (log10 scale)") +
  ylab("Empirical distribution function") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks() +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(
    values = c("Simulated" = "gray70", "Observed" = "Black")
  ) +
  theme(legend.title = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.9, 0.2))

# Compare the different distributions for next year's totals -----------------

distributions_of_total_claims <- bind_rows(
  tibble(
    method = "Parametric Bootstrap (Tweedie MLE)",
    total_claim_amount_in_millions = predicted_totals_for_unit_exposure
  ),
  tibble(
    method = "Nonparametric Bootstrap (Tweedie GLM)",
    total_claim_amount_in_millions = posterior_distribution_samples_for_total_claims$total
  ),
  tibble(
    method = "Closed-form Tweedie Distribution for Weighted Means",
    total_claim_amount_in_millions = predicted_totals_tweedie_sampling_dist
  ),
  tibble(
    method = "Main Effects Tweedie GLM",
    total_claim_amount_in_millions = predicted_totals_from_glm
  )
) %>%
  mutate(
    total_claim_amount_in_millions = total_claim_amount_in_millions / 1e6
  ) %>%
  bind_rows(nonparametric_bootstrap_totals, .)

distributions_of_total_claims

# order the methods by increasing values of the ranges
method_order <- distributions_of_total_claims %>%
  group_by(method) %>%
  summarise(range = diff(range(total_claim_amount_in_millions))) %>%
  arrange(range) %>%
  print() %>%
  pull(method)

distributions_of_total_claims <- distributions_of_total_claims %>%
  mutate(
    method = factor(method, levels = method_order)
  )

# compare the curves
distributions_of_total_claims %>%
  split(.$method) %>%
  map(.f = ~ {
    .x %>%
      pull(total_claim_amount_in_millions) %>%
      quantile(probs = sort(c(seq(0, 1, 0.25), 0.95, 0.99, 0.999))) %>%
      round(., 2)
  })

totals_plot <- distributions_of_total_claims %>%
  ggplot() +
  stat_ecdf(aes(x = total_claim_amount_in_millions,
                linetype = method, group = method), pad = FALSE) +
  geom_vline(xintercept = (6e4 * expected_claim_amount) / 1e6,
             color = "orange", linewidth = 1.2, linetype = "dashed") +
  scale_y_continuous(labels = scales::label_percent()) +
  xlab("Plausible values for the following year's total claim amount in millions, assuming unit exposure") +
  ylab("Empirical distribution function")+
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

totals_plot +
  scale_x_continuous(breaks = seq(3, 54, 3)) +
  guides(linetype = guide_legend(nrow = 3))

totals_plot +
  facet_wrap(~factor(method, levels = method_order)) +
  theme(legend.position = "none") +
  geom_vline(
    data = distributions_of_total_claims %>%
      summarize(m = mean(total_claim_amount_in_millions), .by = "method"),
    aes(xintercept = m, group = method, linetype = method)
  ) +
  scale_x_continuous(breaks = seq(3, 54, 6))

# Use deviance to pick the value of power p -----
seq(1.1, 1.9, 0.1) %>%
  set_names(x = ., nm = ~ .x) %>%
  map_dbl(
    .x = .,
    .f = ~ sum(tweedie::tweedie.dev(claims$ClaimAmount, mu, power = .x))
  ) %>%
  sort() # picks 1.8

set.seed(43)
temp <- tweedie::rtweedie(n = 1e5, mu = 250, phi = 20000, power = 1.6)
mean(temp)
glm(temp ~ 1, family = statmod::tweedie(var.power = 1.3, link.power = 1)) %>%
  fitted() %>%
  tweedie::tweedie.dev(temp, ., power = 1.3) %>%
  sum()

seq(1.1, 1.9, 0.1) %>%
  set_names(x = ., nm = ~ .x) %>%
  map_dbl(
    .x = .,
    .f = ~ {
      temp_glm <- glm(temp ~ 1, family = statmod::tweedie(var.power = .x, link.power = 1))
      sum(tweedie::tweedie.dev(temp, fitted(temp_glm), power = .x))
    }
  ) %>%
  sort() # also picks p = 1.8 despite the true value being 1.6
# 1.8       1.9       1.7       1.6       1.5       1.4       1.3
# 3720120   3833189   4882134   7333243  12000058  20989691  38958037
# 1.2       1.1
# 76584135 159523583

true_prof <- tweedie::tweedie.profile(temp ~ 1, link.power = 1)
true_prof
plot(true_prof)
true_prof$xi.max

smallest_p <- map_chr(
  .x = 1:10,
  .f = ~ {
    set.seed(.x)
    temp_y <- tweedie::rtweedie(n = 1e4, mu = 250, phi = 20000, power = 1.6)
    pow <- seq(1.1, 1.9, 0.1) %>%
      set_names(x = ., nm = ~ .x) %>%
      map_dbl(
        .x = .,
        .f = ~ {
          temp_glm <- glm(temp_y ~ 1,
                          family = statmod::tweedie(var.power = .x, link.power = 1))
          sum(tweedie::tweedie.dev(temp, fitted(temp_glm), power = .x))
        }
      ) %>%
      which.min() %>%
      names() %>%
      pluck(1)

    return(pow)
  }
)

smallest_p %>% as.numeric() %>% table()

# Compare the sampled with the non-sampled -----------------
# sample 10-20 draws from the different models and make a plot for each model
# with the raw sample values overlaid
claims %>%
  filter(ClaimAmount > 0) %>%
  ggplot(aes(x = ClaimAmount)) +
  stat_ecdf() +
  scale_x_log10()

claims_full %>%
  filter(ClaimAmount > 0) %>%
  ggplot(aes(x = ClaimAmount)) +
  stat_ecdf() +
  scale_x_log10()

claims_full_sample_indicator <- claims %>%
  mutate(in_sample = 1) %>%
  select(IDpol, in_sample) %>%
  left_join(x = claims_full, y = ., by = "IDpol") %>%
  mutate(in_sample = replace_na(in_sample, 0))

claims_full_sample_indicator %>% count(in_sample)

# cobalt::bal.plot()
