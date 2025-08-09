
library(tidyverse)
library(mirai)

run_pipeline_in_parallel <- function(seed) {
  `%>%` <- magrittr::`%>%`

  all_cause_cumulative_hazards <- function(t, x1) {
    sum(
      flexsurv::Hweibull(x = t, shape = 1.6, scale = exp(log(5) + log(1.4) * x1), log = FALSE),
      flexsurv::Hllogis(x = t, shape = 5, scale = 5, log = FALSE),
      flexsurv::Hweibull(x = t, shape = 1, scale = 5, log = FALSE)
    )
  }

  cause_specific_probabilities <- function(t, x1) {
    hazards <- c(
      flexsurv::hweibull(x = t, shape = 1.6, scale = exp(log(5) + log(1.4) * x1), log = FALSE),
      flexsurv::hllogis(x = t, shape = 5, scale = 5, log = FALSE),
      flexsurv::hweibull(x = t, shape = 1, scale = 5, log = FALSE)
    )

    hazards / sum(hazards)
  }

  obj_fun <- function(t, u, x1) {
    log(1 - u) + all_cause_cumulative_hazards(t, x1)
  }

  simulate_single_outcome <- function(u, x1, seed = 12345, admin_censor_time = 5, lower_bound_time = 1e-10) {

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
      #set.seed(seed)
      cause <- which.max(stats::rmultinom(1, 1, cause_specific_probabilities(t = time, x1 = x1)))
      iter <- uniroot_obj$iter
      fn_value_at_root <- uniroot_obj$f.root
    }

    dplyr::tibble(u, x1, admin_censor_time, iter, fn_value_at_root, time, cause)
  }

  simulate_single_dataset <- function(n_rows = 1000, prop_x1 = 0.5, seed = 12345) {
    set.seed(seed)
    x1 <- stats::rbinom(n = n_rows, size = 1, prob = prop_x1)
    u <- stats::runif(n_rows)
    purrr::map(
      .x = 1:n_rows,
      .f = ~ simulate_single_outcome(u = u[.x], x1 = x1[.x], seed = seed)
    ) %>%
      purrr::list_rbind()
  }

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

  fit_cause_specific_models <- function(cause_specific_datasets,
                                        distributions = c("weibull", "llogis", "weibull")) {
    purrr::map2(
      .x = cause_specific_datasets,
      .y = distributions,
      .f = ~ flexsurv::flexsurvreg(survival::Surv(time = time, event = event) ~ x1, data = .x, dist = .y))
  }

  extract_coefs <- function(list_models) {
    purrr::map2_dfr(.x = list_models, .y = 1:length(list_models), .f = ~ {
      .x %>%
        # calls flexsurv:::tidy.flexsurvreg()
        flexsurv::tidy() %>%
        dplyr::select(term, estimate, se = std.error) %>%
        dplyr::mutate(cause = .y, .before = term)
    })
  }

  # run the functions
  seed %>%
    simulate_single_dataset(seed = .) %>%
    create_cause_specific_datasets() %>%
    fit_cause_specific_models() %>%
    extract_coefs() %>%
    dplyr::mutate(rep = seed, .before = cause)
}

safely_run_pipeline_in_parallel <- safely(run_pipeline_in_parallel)

daemons(parallelly::availableCores())

# test for small n_rep
system.time({test_parallel_run <- map(.x = 1:10, .f = in_parallel(\(x) safely_run_pipeline_in_parallel(seed = x), safely_run_pipeline_in_parallel = safely_run_pipeline_in_parallel))})

test_parallel_run

# test if two separate runs are identical
foo <- map(
  .x = 1:10,
  .f = in_parallel(\(x) safely_run_pipeline_in_parallel(seed = x),
                   safely_run_pipeline_in_parallel = safely_run_pipeline_in_parallel)) %>%
  pluck(10) %>%
  pluck("result")

bar <- map(
  .x = 1:10,
  .f = in_parallel(\(x) safely_run_pipeline_in_parallel(seed = x),
                   safely_run_pipeline_in_parallel = safely_run_pipeline_in_parallel)) %>%
  pluck(10) %>%
  pluck("result")

identical(foo, bar)

# run for all n_reps
n_reps <- 5000

system.time({
  list_results <- map(
    .x = 1:n_reps,
    .f = in_parallel(\(x) safely_run_pipeline_in_parallel(x),
                     safely_run_pipeline_in_parallel = safely_run_pipeline_in_parallel)
  )
})
#   user   system  elapsed
# 35.955    3.581 1396.126

1396.126 / 60 # 23 mins for the run

daemons(0)

failed_parallel_runs <- keep(
  .x = 1:length(list_results),
  .p = ~ !is.null(list_results[[.x]]$error)
)

failed_parallel_runs

length(failed_parallel_runs) / n_reps

write_rds(failed_parallel_runs, "posts/simulating-and-modelling-competing-risk-data/failed_parallel_runs.rds")

parameter_estimates_parallel <- list_results %>%
  map(.f = ~ pluck(.x, "result")) %>%
  compact() %>%
  list_rbind()

parameter_estimates_parallel

write_rds(parameter_estimates_parallel, "posts/simulating-and-modelling-competing-risk-data/parameter_estimates_parallel.rds")

# code for experimenting with carrier::crate
#
# sim_fn_env <- rlang::env(
#   `%>%` = magrittr::`%>%`,
#   all_cause_cumulative_hazards = all_cause_cumulative_hazards,
#   cause_specific_probabilities = cause_specific_probabilities,
#   obj_fun = obj_fun,
#   simulate_single_outcome = simulate_single_outcome,
#   simulate_single_dataset = simulate_single_dataset,
#   create_cause_specific_datasets = create_cause_specific_datasets,
#   fit_cause_specific_models = fit_cause_specific_models,
#   extract_coefs = extract_coefs,
#   run_pipeline = run_pipeline
# )
#
# sim_crate <- carrier::crate(
#   function(x) run_pipeline(seed = x),
#   .parent_env = sim_fn_env
# )
#
# sim_crate
#
# # check the size difference with using the global env
# carrier::crate(
#   function(seed) run_pipeline(seed),
#   .parent_env = rlang::global_env()
# )
#
# test_parallel_run <- sim_crate(1:10)
#
# test_parallel_run
#
# system.time({test_parallel_run <- sim_crate(1:10)}) # only runs for the first input
#
# test_parallel_run
