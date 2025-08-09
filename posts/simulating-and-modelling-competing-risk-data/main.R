
library(tidyverse)
library(survival)
library(flexsurv)
library(mirai)

theme_set(theme_bw())

base_path <- "posts/simulating-and-modelling-competing-risk-data"
get_path <- function(fname, ext = ".rds") fs::path(base_path, "saved_objects", paste0(fname, ext))
get_path("foo")

# pick simulation parameters -----------
n <- 1e4
set.seed(25)
random_draws <- tibble(
  t1 = rweibull(n, shape = 1.6, scale = 5),
  # scale = 7 when x1 = 1
  x1 = rweibull(n, shape = 1.6, scale = exp(log(5) + log(1.4))),
  t2 = rllogis(n, shape = 5, scale = 5),
  # exponential distribution, same as rexp(n, 1 / 5)
  t3 = rweibull(n, shape = 1, scale = 5)
)

set.seed(25)
latent_times <- tibble(
  # generate group membership
  x1 = rbinom(n, 1, prob = 0.5),
  # generate counterfactual for x1 = 0
  t1_0 = rweibull(n, shape = 1.6, scale = 5),
  # and counterfactual for x1 = 1 (using scale = 7)
  t1_1 = rweibull(n, shape = 1.6, scale = exp(log(5) + log(1.4))),
  # use the consistency equation to realize the latent t1
  t1 =  (x1 * t1_1) + ((1 - x1) * t1_0),
  t2 = rllogis(n, shape = 5, scale = 5),
  # exponential distribution, same as rexp(n, 1 / 5)
  t3 = rweibull(n, shape = 1, scale = 5),
  # administrative censoring
  cens_time = 4,
  # take the minimum of the latent event and censoring times
  time = pmin(t1, t2, t3, cens_time),
  # record which event / censor time is the smallest
  event = case_when(
    cens_time < t1 & cens_time < t2 & cens_time < t3 ~ 0,
    t1 < t2 & t2 < t3 ~ 1,
    t2 < t1 & t2 < t3 ~ 2,
    t3 < t1 & t3 < t2 ~ 3
  )
) %>%
  mutate(across(.cols = everything(), .fns = ~ round(.x, 2)))

true_params <- tribble(
  ~ cause, ~ term, ~ estimate,
  1L, "shape", 1.6,
  1L, "scale", 5,
  1L, "x1", log(1.4),
  2L, "shape", 5,
  2L, "scale", 5,
  2L, "x1", 0,
  3L, "shape", 1,
  3L, "scale", 5,
  3L, "x1", 0
)

true_params

summary(random_draws)

# check the parameterization
# flexsurvreg(Surv(rweibull(3e4, shape = 1, scale = 5), rep(1, 3e4)) ~ 1, dist = "weibull")

# plot the CDF for each distribution
random_draws %>%
  pivot_longer(everything()) %>%
  mutate(
    Distribution = case_when(
      name == "t1" ~ "Event 1: Weibull(shape = 1.6, scale = 5)",
      name == "x1" ~ "Event 1: Weibull(shape = 1.6, scale = 7)",
      name == "t2" ~ "Event 2: Log-logistic(shape = 5, scale = 5)",
      name == "t3" ~ "Event 3: Exponential(rate = 5)"
    )
  ) %>%
  ggplot(aes(x = value, group = Distribution, color = Distribution)) +
  stat_ecdf() +
  xlab("Time (t)") +
  ylab("CDF (Pr[T <= t])") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

# plot the hazards for each setup
time_seq <- seq(0.01, 50, 0.01)
str(time_seq)
hazards <- bind_rows(
  tibble(
    Distribution = "Event 1: Weibull(shape = 1.6, scale = 5)",
    time = time_seq,
    Hazard = hweibull(time_seq, shape = 1.6, scale = 5)
  ),
  tibble(
    Distribution = "Event 1: Weibull(shape = 1.6, scale = 7)",
    time = time_seq,
    Hazard = hweibull(time_seq, shape = 1.6, scale = 7)
  ),
  tibble(
    Distribution = "Event 2: Log-logistic(shape = 5, scale = 5)",
    time = time_seq,
    Hazard = hllogis(time_seq, shape = 5, scale = 5)
  ),
  tibble(
    Distribution = "Event 3: Exponential(rate = 5)",
    time = time_seq,
    Hazard = hweibull(time_seq, shape = 1, scale = 5)
  )
) %>%
  glimpse()

hazards %>%
  ggplot(aes(x = time, y = Hazard, group = Distribution, color = Distribution)) +
  geom_line() +
  xlab("Time") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

# prob(cause | event = 1)
cause_specific_hazards_by_group <- bind_rows(
  hazards %>%
    # for x = 0, drop the cause 1 for x = 1
    filter(Distribution != "Event 1: Weibull(shape = 1.6, scale = 7)") %>%
    mutate(Group = "x1 = 0"),
  hazards %>%
    # for x = 1, drop the cause 1 for x = 0
    filter(Distribution != "Event 1: Weibull(shape = 1.6, scale = 5)") %>%
    mutate(Group = "x1 = 1")
)

relative_contribution_to_overall_hazard <-cause_specific_hazards_by_group %>%
  arrange(Group, time, Distribution) %>%
  group_by(Group, time) %>%
  mutate(overall_hazard = sum(Hazard), p = Hazard / overall_hazard) %>%
  ungroup() %>%
  print()

relative_contribution_to_overall_hazard %>%
  ggplot(aes(x = time, y = p, group = interaction(Group, Distribution), color = Distribution)) +
  geom_line() +
  xlab("Time") +
  ylab("P[Cause | Event]") +
  facet_wrap(~ Group) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

# calculate the overall survival from the overall hazard
overall_survival <- expand_grid(
  x1 = c(0, 1), time = time_seq
) %>%
  mutate(surv = exp(-(
    Hweibull(x = time, shape = 1.6, scale = exp(log(5) + log(1.4) * x1), log = FALSE) +
    Hllogis(x = time, shape = 5, scale = 5, log = FALSE) +
    Hweibull(x = time, shape = 1, scale = 5, log = FALSE)
  )),
  Group = paste0("x1 = ", x1)) %>%
  select(-x1) %>%
  print()

time_step <- time_seq[2] - time_seq[1]

# plot the integrand (hazard * overall survival) per CIF per group
cause_specific_hazards_by_group %>%
  inner_join(y = overall_survival, by = c("Group", "time")) %>%
  #filter(Group == "x1 = 0", Distribution == "Weibull(shape = 1.6, scale = 5)") %>%
  mutate(integrand = Hazard * surv) %>%
  ggplot(aes(x = time, y = integrand, color = Distribution)) +
  geom_line() +
  facet_wrap(~ Group)

# calculate CIFs
cumulative_incidence_curves <- cause_specific_hazards_by_group %>%
  inner_join(y = overall_survival, by = c("Group", "time")) %>%
  group_by(Group, Distribution) %>%
  mutate(CIF = order_by(time, cumsum(Hazard * surv * time_step))) %>%
  ungroup() %>%
  print()

# sanity check on whether values (nearly) sum to 1 at each covariate level and time point
cumulative_incidence_curves %>%
  summarize(s = sum(CIF), .by = c(Group, time)) %>%
  inner_join(y = overall_survival, by = c("Group", "time")) %>%
  mutate(total = s + surv) %>%
  pull(total) %>%
  summary()

# plot the CIFs
# P[T <= t, cause = k | Event]
cumulative_incidence_curves_plot_data <- cumulative_incidence_curves %>%
  select(-Hazard) %>%
  pivot_wider(id_cols = c(time, Group, surv), names_from = Distribution, values_from = CIF) %>%
  rename(`Event-free` = surv) %>%
  pivot_longer(cols = -c(time, Group), names_to = "Event", values_to = "CIF") %>%
  filter(complete.cases(.)) %>%
  print()

cumulative_incidence_curves_plot_data %>%
  ggplot(aes(x = time, y = CIF, linetype = Event)) +
  geom_line() +
  xlab("Time") +
  ylab("P[T <= t, cause = k]") +
  facet_wrap(~ Group) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  guides(linetype = guide_legend(nrow = 3))

# plot the normalized CIFs conditional on event
# so ignores the no even state
# P[T <= t, cause = k | Event]
cumulative_incidence_curves %>%
  mutate(normalized_CIF = CIF / sum(CIF), .by = c(Group, time)) %>%
  # arrange(Group, time) %>%
  # print()
  ggplot(aes(x = time, y = normalized_CIF, color = Distribution)) +
  geom_line() +
  xlab("Time") +
  ylab("P[T <= t, cause = k | Event]") +
  facet_wrap(~ Group) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

# simulate data from latent failure time model --------

# implement the simulation approach from Beyersmann et al -------
all_cause_cumulative_hazards <- function(t, x1) {
  sum(
    flexsurv::Hweibull(x = t, shape = 1.6, scale = exp(log(5) + log(1.4) * x1), log = FALSE),
    flexsurv::Hllogis(x = t, shape = 5, scale = 5, log = FALSE),
    flexsurv::Hweibull(x = t, shape = 1, scale = 5, log = FALSE)
  )
}

all_cause_cumulative_hazards(t = 0, x1 = 0)
all_cause_cumulative_hazards(t = 10, x1 = 0)
all_cause_cumulative_hazards(t = 10, x1 = 1)

cause_specific_probabilities <- function(t, x1) {
  hazards <- c(
    flexsurv::hweibull(x = t, shape = 1.6, scale = exp(log(5) + log(1.4) * x1), log = FALSE),
    flexsurv::hllogis(x = t, shape = 5, scale = 5, log = FALSE),
    flexsurv::hweibull(x = t, shape = 1, scale = 5, log = FALSE)
  )

  hazards / sum(hazards)
}

cause_specific_probabilities(t = 0, x1 = 0)
cause_specific_probabilities(t = 10, x1 = 0)
cause_specific_probabilities(t = 10, x1 = 1)

obj_fun <- function(t, u, x1) {
  log(1 - u) + all_cause_cumulative_hazards(t, x1)
}

obj_fun(t = 10, u = 0.9, x1 = 0)

survsim::crisk.ncens.sim

# to pick administrative censoring
map(random_draws, .f = ~ quantile(.x, probs = c(0.9, 0.95, 0.99, 1.0)))

simulate_single_outcome <- function(u, x1, admin_censor_time = 5, lower_bound_time = 1e-10) {

  # check whether simulated time would lie outside the time interval we want to solve for
  # in this case the function would have the same sign at the interval endpoints
  # if the function has different signs at the endpoints then a root is within the interval
  same_sign <- (obj_fun(u = u, t = 0.00001, x1 = x1) * obj_fun(u = u, t = admin_censor_time, x1 = x1)) > 0

  if (same_sign) {
    time <- admin_censor_time
    cause <- 0
    iter <- NA_real_
    fn_value_at_root <- NA_real_
  } else {
    # for a handful of cases (out of 100,000s) we can end up with t = 0 from a distribution that doesn't support it
    # https://stats.stackexchange.com/questions/176376/invalid-survival-times-for-this-distribution
    uniroot_obj <- stats::uniroot(f = obj_fun, interval = c(lower_bound_time, admin_censor_time), u = u, x1 = x1)
    time <- uniroot_obj$root
    cause <- which.max(stats::rmultinom(1, 1, cause_specific_probabilities(t = time, x1 = x1)))
    iter <- uniroot_obj$iter
    fn_value_at_root <- uniroot_obj$f.root
  }

  dplyr::tibble(u, x1, admin_censor_time, iter, fn_value_at_root, time, cause)
}

simulate_single_outcome(u = 0.2, x1 = 1, admin_censor_time = 12)
simulate_single_outcome(u = 0.9, x1 = 1, admin_censor_time = 12)
simulate_single_outcome(u = 0.9, x1 = 1, admin_censor_time = 12)

# check if setting seed work correctly
map_dfr(.x = 1:10, .f = ~ simulate_single_outcome(u = 0.3, x1 = 1))
map_dfr(.x = 1:10, .f = ~ {set.seed(43); simulate_single_outcome(u = 0.3, x1 = 1)})
map_dfr(.x = 1:10, .f = ~ {set.seed(.x); simulate_single_outcome(u = 0.3, x1 = 1)})
set.seed(43)
map_dfr(.x = 1:10, .f = ~ simulate_single_outcome(u = 0.3, x1 = 1))
identical(
  map_dfr(.x = 1:10, .f = ~ {set.seed(.x); simulate_single_outcome(u = 0.3, x1 = 1)}),
  map_dfr(.x = 1:10, .f = ~ {set.seed(.x); simulate_single_outcome(u = 0.3, x1 = 1)})
)

simulate_single_dataset <- function(n_rows = 1000, prop_x1 = 0.5, seed = 12345) {
  if(!is.null(seed)) set.seed(seed)
  x1 <- stats::rbinom(n = n_rows, size = 1, prob = prop_x1)
  u <- stats::runif(n_rows)
  purrr::map(
    .x = 1:n_rows,
    .f = ~ simulate_single_outcome(u = u[.x], x1 = x1[.x])
  ) %>%
    purrr::list_rbind()
}

identical(simulate_single_dataset(n_rows = 100, seed = 1), simulate_single_dataset(n_rows = 100, seed = 1)) # should be true due to same seed
identical(simulate_single_dataset(n_rows = 100, seed = 1), simulate_single_dataset(n_rows = 100, seed = 2)) # should be false due to different seeds

simulate_single_dataset(n_rows = 100, seed = 1)
simulate_single_dataset(n_rows = 100, seed = 2)

# seed works correctly
# cor(simulate_single_dataset(n_rows = 1e4, seed = 1)$time, simulate_single_dataset(n_rows = 1e4, seed = 2)$time)
# [1] 0.003499799

# p(censored)
# simulate_single_dataset(n_rows = 1000, seed = 2) %>%
#   count(cause) %>%
#   mutate(p = n / sum(n))
# # # A tibble: 4 × 3
# #   cause     n     p
# #   <dbl> <int> <dbl>
# # 1     0    94 0.094
# # 2     1   279 0.279
# # 3     2   148 0.148
# # 4     3   479 0.479
# with another seed
# simulate_single_dataset(n_rows = 1000, seed = 1) %>%
#   count(cause) %>%
#   mutate(p = n / sum(n))
# # # A tibble: 4 × 3
# #   cause     n     p
# #   <dbl> <int> <dbl>
# # 1     0    87 0.087
# # 2     1   300 0.3
# # 3     2   128 0.128
# # 4     3   485 0.485

create_cause_specific_datasets <- function(data) {
  purrr::map(
    .x = 1:3,
    .f = ~ {
      data %>%
        dplyr::select(x1, time, cause) %>%
        dplyr::mutate(event = as.numeric(cause == .x))
    }
  )
}

simulate_single_dataset() %>%
  create_cause_specific_datasets()

set.seed(5)
flexsurvreg(Surv(rweibull(3e4, shape = 1, scale = 5), rep(1, 3e4)) ~ rbinom(3e4, 1, 0.5), dist = "weibull") %>%
  flexsurv:::tidy.flexsurvreg(transform = "baseline.real") %>%
  mutate(exp = exp(estimate)) # compare these with the ones from print.flexsurvreg

fit_cause_specific_models <- function(cause_specific_datasets, distributions = c("weibull", "llogis", "weibull")) {
  purrr::map2(
    .x = cause_specific_datasets,
    .y = distributions,
    .f = ~ flexsurv::flexsurvreg(survival::Surv(time = time, event = event) ~ x1, data = .x, dist = .y))
}

simulate_single_dataset(seed = 1) %>%
  create_cause_specific_datasets() %>%
  fit_cause_specific_models()

extract_coefs <- function(list_models) {
  purrr::map2_dfr(.x = list_models, .y = 1:length(list_models), .f = ~ {
    .x %>%
      # calls flexsurv:::tidy.flexsurvreg()
      flexsurv::tidy() %>%
      dplyr::select(term, estimate, se = std.error) %>%
      dplyr::mutate(cause = .y, .before = term)
  })
}

simulate_single_dataset(seed = 1) %>%
  create_cause_specific_datasets() %>%
  fit_cause_specific_models() %>%
  extract_coefs()

list_models <- simulate_single_dataset(seed = 1) %>%
  create_cause_specific_datasets() %>%
  fit_cause_specific_models()

estimate_CIF <- function(seed, time_seq, covariate_level = data.frame(x1 = 0)) {
  list_models <- seed %>%
    simulate_single_dataset(seed = .) %>%
    create_cause_specific_datasets() %>%
    fit_cause_specific_models()

  transition_matrix <- matrix(data = c(c(NA, 1, 2, 3), rep(NA, 12)), nrow = 4, ncol = 4, byrow = TRUE)

  multi_state_object <- flexsurv::fmsm(list_models[[1]], list_models[[2]], list_models[[3]], trans = transition_matrix)

  flexsurv::pmatrix.fs(multi_state_object, trans = transition_matrix, t = time_seq, newdata = covariate_level, tidy = TRUE) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(start == 1) %>%
    dplyr::select(-start) %>%
    tidyr::pivot_longer(cols = -time, names_to = "Distribution", values_to = "CIF")
}

estimate_CIF(seed = 3, time_seq = 1:10)

safely_estimate_CIF <- safely(estimate_CIF)

run_pipeline <- function(seed) {
  seed %>%
    simulate_single_dataset(seed = .) %>%
    create_cause_specific_datasets() %>%
    fit_cause_specific_models() %>%
    extract_coefs() %>%
    dplyr::mutate(rep = seed, .before = cause)
}

run_pipeline(seed = 3)

# estimate approx. runtime --------------
# system.time({test_sequential_run <- map(.x = 1:10, .f = ~ run_pipeline(seed = .x))})
#   user  system elapsed
# 16.793   2.649  15.600
# test_sequential_run %>% head(3)

# 10 runs take about 15.6 seconds
# so to run 5000 sims it would take
((15.6 / 10) * 1000) / (60 * 60) # ~ 2.2 hours for 5000, 26 mins (0.4333 hrs) for 1000, 13 mins (0.21666 hrs) for 500

(2.2 * 60) / 20
# using 20 workers for parallelization would take about 7 mins (note it took about 23 mins using purrr::in_parallel)
# (2.2 (hours) * 60 (minutes / hour)) / 20 (cores)

# run pipeline assess unbiasedness of estimates --------
# fails for seed = 125
run_pipeline(seed = 125)

safely_run_pipeline <- safely(run_pipeline)

safely_run_pipeline(seed = 125)

n_reps <- 5000

# list_parameter_estimates <- 1:n_reps %>%
#   map(.f = ~ {
#     if (.x %% 10 == 0) {
#       cat(glue::glue("rep = {.x}"), sep = "\n")
#     }
#     safely_run_pipeline(seed = .x)
#   })
#
# failed_runs <- keep(
#   .x = 1:length(list_parameter_estimates),
#   .p = ~ !is.null(list_parameter_estimates[[.x]]$error)
# )
#
# write_rds(failed_runs, get_path("failed_runs"))
#
# parameter_estimates <- list_parameter_estimates %>%
#   map(.f = ~ pluck(.x, "result")) %>%
#   compact() %>%
#   list_rbind()
#
# write_rds(parameter_estimates, get_path("parameter_estimates"))

# read results back in
failed_runs <- read_rds(get_path("failed_runs"))
failed_runs
# p(failed runs)
length(failed_runs) / n_reps

failed_parallel_runs <- read_rds(get_path("failed_parallel_runs"))
failed_parallel_runs
length(failed_parallel_runs) / n_reps

parameter_estimates <- read_rds(get_path("parameter_estimates")) %>%
  print()

parameter_estimates

parameter_estimates_parallel <- read_rds(get_path("parameter_estimates_parallel")) %>%
  print()

# pull in parallel run estimates
# compare prop failed runs
# and add to the density plot

# density plot of coefficients
bind_rows(
  parameter_estimates %>%
    mutate(Run = glue::glue("Sequential ({round(100 * length(failed_runs) / n_reps, 2)}% failed runs)")),
  parameter_estimates_parallel %>%
    mutate(Run = glue::glue("Parallel ({round(100 * length(failed_parallel_runs) / n_reps, 2)}% failed runs)"))
) %>%
  ggplot(aes(x = estimate, group = interaction(Run, cause, term), color = Run)) +
  geom_density() +
  geom_vline(data = true_params, aes(xintercept = estimate)) +
  xlab("Sampling Distribution of Simulation Parameters (5,000 datasets with n = 1,000)") +
  #facet_grid(rows = vars(cause), cols = vars(term), scales = "free", axes = "all")
  ggh4x::facet_grid2(rows = vars(cause), cols = vars(term), scales = "free", independent = TRUE) +
  theme(legend.position = "bottom")

# plot 20 CIFs for cause 1 from parametric and non-parametric fits with the true CIF overlaid -------
time_seq %>% str()
time_seq %>% summary()

estimate_CIF(seed = 3, time_seq = time_seq[time_seq  <= 10.0]) %>%
  summary()

time_seq_subset <- time_seq %>% keep(.p = ~ .x <= 10.0)

# sampling_distribution_CIFs <- 1:100 %>%
#   map(
#     .f = ~ {
#       print(glue::glue("rep: {.x}"))
#       safely_estimate_CIF(seed = .x, time_seq = time_seq_subset, covariate_level = data.frame(x1 = 0))
#     }
#   ) %>%
#   imap(
#     .f = ~ {
#       res <- pluck(.x, "result")
#       if (!is.null(res)) {
#         res %>%
#           mutate(rep = .y, .before = time)
#       } else {
#         res
#       }
#     }
#   ) %>%
#   compact() %>%
#   list_rbind() %>%
#   mutate(
#     Event = case_when(
#       Distribution == "V1" ~ "Event-free",
#       Distribution == "V2" ~ "Event 1: Weibull(shape = 1.6, scale = 5)",
#       Distribution == "V3" ~ "Event 2: Log-logistic(shape = 5, scale = 5)",
#       Distribution == "V4" ~ "Event 3: Exponential(rate = 5)"
#     ),
#     Event = fct_relevel(factor(Event), "Event-free", after = 0)
#   )
#
# sampling_distribution_CIFs %>%
#   print(n = 5)
#
# write_rds(sampling_distribution_CIFs, get_path("sampling_distribution_CIFs"))

sampling_distribution_CIFs <- read_rds(get_path("sampling_distribution_CIFs"))

sampling_distribution_CIFs_mean <- sampling_distribution_CIFs %>%
  summarise(CIF = mean(CIF), .by = c(Event, time))

sampling_distribution_CIF_plot <- sampling_distribution_CIFs %>%
  ggplot(aes(x = time, y = CIF, group = rep)) +
  geom_line(color = "gray80") +
  geom_line(
    data = cumulative_incidence_curves_plot_data %>%
      filter(Group == "x1 = 0", time <= 10.0) %>%
      mutate(Event = fct_relevel(factor(Event), "Event-free", after = 0)),
    aes(x = time, y = CIF),
    color = "gray20", linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = sampling_distribution_CIFs_mean,
    aes(x = time, y = CIF),
    color = "gray20", linewidth = 1.5, linetype = "dotted",
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "gray20") +
  xlab(glue::glue("Time (t)\nSampling distribution of curves ",
                  "estimated from 100 datasets with ",
                  "n = 1,000 (light gray); \ntrue curve ",
                  "(gray, solid); averaged Curve (gray, dotted);",
                  "\nmaximum Observed Time (gray, dashed)")) +
  ylab("Probability of being in state") +
  facet_wrap(~ Event, scales = "free_y")

sampling_distribution_CIF_plot

post_plot <- sampling_distribution_CIF_plot + xlab("Time (t)")
post_plot
ggsave(plot = post_plot, filename = "cif_plot.png",
       device = "png", units = "px", width = 900, height = 800, dpi = 130,
       path = base_path)

# simulate some data from survsim to benchmark speed -------
survsim::crisk.ncens.sim

survsim_sim_fn <- function(n = 1000) {
  survsim::crisk.sim(
    n = n,
    foltime = 5,
    # TTE distr
    nsit = 3,
    dist.ev = c("weibull", "llogistic", "weibull"),
    # anc.ev is 1 / shape, exp(beta0) is scale in flexsurv weibull AFT
    # anc.ev is shape, beta0.ev is shape * log(scale) in flexsurv llogis
    anc.ev = c(1 / 1.6, 5, 1),
    beta0.ev = c(log(5), 5 * log(5), log(5)),
    # covariate process
    x = list(c("bern", 0.5)),
    beta = list(c(1.4), c(0), c(0)),
    # censoring process, assume constant here to only apply administrative censoring
    dist.cens = "unif", beta0.cens = 500, anc.cens = 500
  ) %>%
    as_tibble()
}

# survsim_sim_fn()
#
# survsim_bench <- bench::mark(
#   survsim_sim_fn(), iterations = 10,
#   check = FALSE, time_unit = "s", memory = FALSE
# )
# write_rds(survsim_bench, file = get_path("survsim_bench"))

survsim_bench <- read_rds(get_path("survsim_bench"))

survsim_bench
survsim_bench %>% pluck("time", 1)

# my_implementation_bench <- bench::mark(
#   run_pipeline(sample(1:100, size = 1, replace = TRUE)),
#   iterations = 10, check = FALSE, time_unit = "s", memory = FALSE
# )
# write_rds(my_implementation_bench, file = get_path("my_implementation_bench"))

my_implementation_bench <- read_rds(file = get_path("my_implementation_bench"))

my_implementation_bench
my_implementation_bench %>% pluck("time", 1)

# dot plot of the run times
tibble(
  time = c(survsim_bench %>% pluck("time", 1), my_implementation_bench %>% pluck("time", 1)),
  group = rep(c("survsim", "handrolled"), each = 10)
) %>%
  ggplot(aes(x = time, y = group, fill = group, color = group)) +
  geom_dotplot(binwidth = 0.05,  binaxis = "x", method = "histodot", binpositions = "bygroup", stackdir = "centerwhole") +
  #geom_point(size = 3, position = position_jitter(height = 0.1, seed = 4)) +
  #theme(legend.position = "bottom", legend.title = element_blank()) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(1, 4, 1), labels = seq(1, 4, 1), limits = c(1, 4)) +
  xlab("Time (in seconds)") +
  ylab("")

# eCDFs for a large draw overlaid to compare the draws - they should be overlapping
# set.seed(23)
# survsim_sample <- survsim_sim_fn(n = 10000)
# set.seed(23)
# handrolled_sample <- simulate_single_dataset(n_rows = 10000, seed = NULL)
#
# write_rds(survsim_sample, file = get_path("survsim_sample"))
# write_rds(handrolled_sample, file = get_path("handrolled_sample"))

survsim_sample <- read_rds(file = get_path("survsim_sample"))
handrolled_sample <- read_rds(file = get_path("handrolled_sample"))

survsim_sample %>% count(cause)
handrolled_sample %>% count(cause)

bind_rows(
  survsim_sample %>%
    select(time) %>%
    mutate(group = "survsim"),
  handrolled_sample %>%
    select(time) %>%
    mutate(group = "handrolled")
) %>%
  mutate(rep = sample(1:100, size = 20000, replace = TRUE)) %>%
  ggplot(aes(x = time, group = interaction(group, rep), color = group)) +
  stat_ecdf(pad = FALSE, alpha = 0.8) +
  # ggplot(aes(x = time, group = group, color = group)) +
  # stat_ecdf(pad = FALSE, alpha = 1) +
  xlab("Simulated Times") +
  ylab("Empirical Distribution Function") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = c("orange", "gray10")) #+
  #facet_wrap(~ rep)

# compare uniroot with optimize and with optim ----------------
obj_fun

tibble(
  time = time_seq[time_seq < 10.0],
  obj_val = map_dbl(
    .x = time,
    .f = ~ obj_fun(.x, u = 0.8, x1 = 0)
  ),
  obj_val_sq = obj_val^2
) %>%
  ggplot(aes(x = time, y = obj_val_sq)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 3.771, linetype = "dashed") +
  xlab("Time (t)") +
  ylab("Objective function")

sq_obj_fun <- function(t, u, x1 = 0) {
  obj_fun(t=t, u=u, x1=x1)^2
}

all_cause_hazards <- function(t, x1) {
  sum(
    flexsurv::hweibull(x = t, shape = 1.6, scale = exp(log(5) + log(1.4) * x1), log = FALSE),
    flexsurv::hllogis(x = t, shape = 5, scale = 5, log = FALSE),
    flexsurv::hweibull(x = t, shape = 1, scale = 5, log = FALSE)
  )
}

sq_obj_fun_grad <- function(t, u, x1 = 0) {
  2 * obj_fun(t=t, u=u, x1=x1) * all_cause_hazards(t=t, x1=x1)
}

uniroot(f = obj_fun, interval = c(1e-05, 100), u = 0.3, x1 = 0)
optimise(f = sq_obj_fun, interval = c(1e-05, 100), u = 0.3, x1 = 0, maximum = FALSE)
optim(par = c(5), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = 0.3, x1 = 0, lower = 1e-05, upper = 100, method = "L-BFGS-B")
optim(par = c(3), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = 0.3, x1 = 0, lower = 1e-05, upper = 100, method = "L-BFGS-B")
optim(par = c(2), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = 0.3, x1 = 0, lower = 1e-05, upper = 100, method = "L-BFGS-B")
optim(par = c(20), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = 0.3, x1 = 0, lower = 1e-05, upper = 100, method = "L-BFGS-B")
optim(par = c(5), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = 0.3, x1 = 0, lower = 1e-05, upper = 10, method = "L-BFGS-B")
optim(par = c(5), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = 0.3, x1 = 0, lower = 1e-05, upper = 10, method = "L-BFGS-B", control = list(pgtol = 1e-5))

u_vec <- seq(0.05, 0.95, 0.05)
# large interval
bench::mark(
  map_dbl(.x = u_vec, .f = ~ uniroot(f = obj_fun, interval = c(1e-05, 100), u = .x, x1 = 0)$root),
  map_dbl(.x = u_vec, .f = ~ optimise(f = sq_obj_fun, interval = c(1e-05, 100), u = .x, x1 = 0, maximum = FALSE)$minimum),
  map_dbl(.x = u_vec, .f = ~ optim(par = c(5), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = .x, x1 = 0, lower = 1e-05, upper = 100, method = "L-BFGS-B")$par),
  iterations = 30, check = FALSE, memory = FALSE, filter_gc = TRUE
) %>% select(min, median, total_time)
# # A tibble: 3 × 3
#        min   median total_time
#   <bch:tm> <bch:tm>   <bch:tm>
# 1   17.9ms   18.7ms      453ms
# 2   18.1ms   19.3ms      459ms
# 3   30.1ms   31.9ms      572ms

# narrower interval
bench::mark(
  map_dbl(.x = u_vec, .f = ~ uniroot(f = obj_fun, interval = c(1e-05, 12), u = .x, x1 = 0)$root),
  map_dbl(.x = u_vec, .f = ~ optimise(f = sq_obj_fun, interval = c(1e-05, 12), u = .x, x1 = 0, maximum = FALSE)$minimum),
  map_dbl(.x = u_vec, .f = ~ optim(par = c(5), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = .x, x1 = 0, lower = 1e-05, upper = 12, method = "L-BFGS-B")$par),
  iterations = 30, check = FALSE, memory = FALSE, filter_gc = TRUE
) %>% select(min, median, total_time)
# # A tibble: 3 × 3
#        min   median total_time
#   <bch:tm> <bch:tm>   <bch:tm>
# 1     12ms   12.9ms      357ms
# 2   14.4ms     15ms      395ms
# 3   31.3ms   33.3ms      623ms

# setting same tol
bench::mark(
  map_dbl(.x = u_vec, .f = ~ uniroot(f = obj_fun, interval = c(1e-05, 12), u = .x, x1 = 0)$root, tol = 1e-5),
  map_dbl(.x = u_vec, .f = ~ optimise(f = sq_obj_fun, interval = c(1e-05, 12), u = .x, x1 = 0, maximum = FALSE)$minimum, tol = 1e-5),
  map_dbl(.x = u_vec, .f = ~ optim(par = c(5), fn = sq_obj_fun, gr = sq_obj_fun_grad, u = .x, x1 = 0, lower = 1e-05, upper = 12, method = "L-BFGS-B")$par,
          control = list(pgtol = 1e-5)),
  iterations = 30, check = FALSE, memory = FALSE, filter_gc = TRUE
) %>% select(min, median, total_time)
# # A tibble: 3 × 3
#        min   median total_time
#   <bch:tm> <bch:tm>   <bch:tm>
# 1   11.9ms   12.4ms      346ms
# 2     14ms   14.8ms      359ms
# 3   30.6ms   31.7ms      537ms

# check the estimated roots
uniroot_est <- map_dbl(.x = u_vec, .f = ~ uniroot(f = obj_fun, interval = c(1e-05, 100), u = .x, x1 = 0)$root)
optimize_est <- map_dbl(.x = u_vec, .f = ~ uniroot(f = obj_fun, interval = c(1e-05, 100), u = .x, x1 = 0)$root)

tibble(u = u_vec, uniroot_est)
