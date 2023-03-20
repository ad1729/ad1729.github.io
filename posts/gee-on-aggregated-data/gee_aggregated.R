
library(tidyverse)
library(furrr)
library(geepack)
# library(gee)
# library(geeM)

theme_set(theme_bw())

# compare estimates from the ohio dataset ----

data(ohio)
summary(ohio)
glimpse(ohio)

mod1 <- geeglm(resp ~ smoke, data = ohio,
               family = "binomial", id = id, corstr = "exchangeable")
summary(mod1)

# geem(resp ~ smoke, data = ohio,
#      family = "binomial", id = id, corstr = "exchangeable") %>% summary()

agg_df <- ohio %>%
  group_by(id, smoke) %>%
  summarize(n = n(), resp = sum(resp), .groups = "drop") %>%
  mutate(resp_p = resp / n) %>%
  glimpse()

mod2 <- geeglm(resp_p ~ smoke, weights = n, data = agg_df,
               family = "binomial", id = id, corstr = "exchangeable")
summary(mod2)

bind_rows(
  broom::tidy(mod1) %>% mutate(data = 'raw data'),
  broom::tidy(mod2) %>% mutate(data = 'aggregated data')
) %>%
  arrange(term)

# simulations ----

simulate_data_and_get_estimates <- function(n, n_clus, seed, return_data = FALSE) {
  set.seed(seed)
  # cluster random intercept
  rand_intercept <- rnorm(n_clus, mean = 0, sd = 1)

  set.seed(seed)
  sim_data <- map_dfr(.x = 1:n_clus, .f = ~ {
    tibble(
      id = .x,
      x = rbinom(round(n / n_clus), size = 1, prob = runif(1, 0.4, 0.6)),
      y = rbinom(round(n / n_clus), size = 1, prob = plogis(-1.4 + 0.3 * x + rand_intercept[[.x]]))
    )
  })

  agg_data <- sim_data %>%
    group_by(id, x) %>%
    summarize(n = n(), y = sum(y), .groups = "drop") %>%
    mutate(y_prob = y / n)

  if (return_data) {
    return(list(sim_data = sim_data, agg_data = agg_data))
  }

  # fit model to raw data
  t0 <- Sys.time()
  mod1_sim <- geeglm(
    y ~ x,
    data = sim_data,
    family = "binomial", id = id, corstr = "exchangeable"
  )
  t1 <- Sys.time()
  t1 <- as.numeric(difftime(t1, t0, units = "secs"))

  # fit GLM model to raw data
  mod1_sim_glm <- glm(y ~ factor(id) + x - 1,
                      data = sim_data, family = "binomial")

  # fit model to aggregated data
  t2 <- Sys.time()
  mod2_sim <- geeglm(
    y_prob ~ x,
    data = agg_data, weights = n,
    family = "binomial", id = id, corstr = "exchangeable"
  )
  t3 <- Sys.time()
  t3 <- as.numeric(difftime(t3, t2, units = "secs"))

  # fit GLM model to aggregated data
  mod2_sim_glm <- glm(y_prob ~ factor(id) + x - 1, data = agg_data,
                      weights = n, family = "binomial")

  results <- bind_rows(
    # GEE results
    broom::tidy(mod1_sim) %>%
      mutate(data = 'raw data', time_secs = t1, estimator = "GEE"),
    broom::tidy(mod2_sim) %>%
      mutate(data = 'aggregated data', time_secs = t3, estimator = "GEE"),
    # GLM results
    broom::tidy(mod1_sim_glm) %>%
      mutate(data = 'raw data', time_secs = NA_real_, estimator = "GLM"),
    broom::tidy(mod2_sim_glm) %>%
      mutate(data = 'aggregated data', time_secs = NA_real_, estimator = "GLM")
  ) %>%
    mutate(n = n, n_clus = n_clus, seed = seed) %>%
    filter(!stringr::str_detect(term, pattern = "id"))

  return(results)
}

simulate_data_and_get_estimates(n = 100, n_clus = 5, seed = 10)

simulation_parameters <- expand_grid(
  n = c(100, 500, 1000, 2500, 5000, 7500, 10000),
  n_clus = c(5, 50, 250),
  seed = 1:20
) %>%
  # remove runs where number of clusters is larger than total sample size
  filter(n_clus < n)

# set up multicore processing
plan(multisession, workers = 15)

# parallelize purrr::pmap_dfr by using the {furrr} package
simulation_results <- future_pmap_dfr(
  .l = simulation_parameters,
  # have to pass args as ..1 or ..2, else it fails
  .f = ~ simulate_data_and_get_estimates(n = ..1, n_clus = ..2, seed = ..3),
  .options = furrr_options(seed = NULL)
)

plan(sequential) # turn off multicore

simulation_results %>% glimpse()

View(simulation_results)

# simulation_results %>%
#   write_csv(file = fs::path("posts", "gee-on-aggregated-data", "gee_glm_sim.csv"))

# analyze results ----

simulation_results <- readr::read_csv(
  file = fs::path("posts", "gee-on-aggregated-data", "gee_glm_sim.csv")
) %>%
  mutate(
    n_clus_raw = n_clus,
    n_clus = factor(paste0("# clusters: ", n_clus),
                    levels = paste0("# clusters: ", c(5, 50, 250))),
    term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "x" ~ "Slope",
      TRUE ~ term
    ))

simulation_results %>% glimpse()

simulation_results %>% count(n_clus)

# plot runtimes
simulation_results %>% glimpse()

runtime_data <- simulation_results %>%
  filter(data == "raw data", estimator == "GEE") %>%
  distinct(seed, time_secs, n, n_clus) %>%
  glimpse()

# this doesn't produce slope estimate of 3 because
# n is small enough here that n and n^2 terms are large enough relative to n^3
runtime_data %>%
  group_by(n_clus) %>%
  nest() %>%
  mutate(n_sims = map_int(.x = data, .f = ~ nrow(.x)),
         mod = map(.x = data,
                         .f = ~ {
                           .x %>%
                             lm(log(time_secs) ~ log(n), data = .) %>%
                             broom::tidy(conf.int = TRUE) %>%
                             filter(term == "log(n)")
                         }),
         est = map_dbl(mod, ~ pull(.x, estimate)),
         #se = map_dbl(mod, ~ pull(.x, std.error)),
         conf.low = map_dbl(mod, ~ pull(.x, conf.low)),
         conf.high = map_dbl(mod, ~ pull(.x, conf.high))) %>%
  select(-data, -mod) %>%
  print()

runtime_data %>%
  filter(n_clus == "# clusters: 5") %>%
  lm(time_secs ~ poly(n, degree = 3, raw = TRUE), data = .) %>%
  broom::tidy()

runtime_data %>%
  filter(n_clus == "# clusters: 5") %>%
  lm(time_secs ~ poly(n, degree = 3, raw = TRUE), data = .) %>%
  broom::glance() %>%
  glimpse()

# not sure this makes sense as we're trying to
# fit a 4 degree polynomial through 5 (design) points
runtime_plots <- simulation_results %>%
  filter(data == "raw data", estimator == "GEE") %>%
  # distinct here since the data frame contains a row for
  # each of the intercept and slope terms but the runtimes are the same
  # for both these terms
  distinct(seed, time_secs, n, n_clus) %>%
  ggplot(aes(x = n, y = time_secs, group = n_clus)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, degree = 3)) +
  geom_smooth(method = "lm", formula = y ~ I(x ^ 3), color = "forestgreen") +
  stat_summary(aes(group = interaction(n_clus, n)), geom = "point",
               fun = mean, color = "darkorange", size = 2) +
  xlab("Total sample size (n)") +
  ylab("Runtime (in seconds)") +
  facet_wrap(~ n_clus, scales = "free")

plot(runtime_plots)

plotly::ggplotly(runtime_plots)

# some convergence issues with small # of clusters
simulation_results %>%
  group_by(estimator, n_clus) %>%
  summarize(max(std.error))

simulation_results %>%
  filter(estimator == "GEE") %>%
  arrange(desc(std.error)) %>%
  head(n = 10)

simulation_results <- simulation_results %>%
  mutate(high_se = as.numeric(std.error > 5))

simulation_results %>%
  filter(estimator == "GEE") %>%
  select(-estimate, -std.error, -statistic, -p.value, -time_secs) %>%
  # put intercept and slope estimates in one row
  pivot_wider(everything(), names_from = term, values_from = high_se) %>%
  # if at least one of the intercept or slope terms have a very high se
  mutate(high_se = pmin(Intercept + Slope, 1)) %>%
  group_by(data, n_clus, n) %>%
  summarize(n_total = n(), n_failed = sum(high_se), .groups = "drop_last") %>%
  #print(n = Inf) %>%
  summarize(n_total = sum(n_total), n_failed = sum(n_failed), .groups = "drop") %>%
  mutate(percent_failed = 100 * n_failed / n_total)

slope_plot_data <- simulation_results %>%
  filter(term == "Slope", high_se == 0) %>%
  select(seed, term, data, estimate, n_clus, n, estimator) %>%
  tidyr::pivot_wider(id_cols = c(seed, term, n_clus, n, estimator),
                     names_from = data, values_from = estimate) %>%
  mutate(`Total sample size` = factor(n),
         id = interaction(estimator, n_clus, sep = ", ", lex.order = TRUE))

slope_plot <- slope_plot_data %>%
  filter(estimator == "GEE") %>%
  ggplot(aes(x = `raw data`, y = `aggregated data`, color = `Total sample size`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  # geom_point(data = tibble(x = 0.3, y = 0.3),
  #            aes(x = x, y = y),
  #            color = "black", size = 3, inherit.aes = FALSE) +
  xlab("Slope coefficient from the full dataset") +
  ylab("Slope coefficient from the aggregated dataset") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1)) +
  facet_wrap(~ id, ncol = 3, scales = "fixed")

plotly::ggplotly(p = slope_plot)

# GLM slope slope plot
slope_plot_data %>%
  filter(estimator == "GLM") %>%
  ggplot(aes(x = `raw data`, y = `aggregated data`, color = `Total sample size`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  # geom_point(data = tibble(x = 0.3, y = 0.3),
  #            aes(x = x, y = y),
  #            color = "black", size = 3, inherit.aes = FALSE) +
  xlab("Slope coefficient from the full dataset") +
  ylab("Slope coefficient from the aggregated dataset") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1)) +
  facet_wrap(~ id, ncol = 3, scales = "fixed")

dist_plot_data <- simulation_results %>%
  mutate(
    ci_lower = estimate - (qnorm(0.975) * std.error),
    ci_upper = estimate + (qnorm(0.975) * std.error),
    ci_width = ci_upper - ci_lower,
    ci_upper_half_width = ci_upper - estimate
  ) %>%
  filter(term == "Slope", high_se == 0, estimator == "GEE") %>%
  select(
    Estimate = estimate,
    `Std. Error` = std.error,
    `95% CI LL` = ci_lower,
    `95% CI UL` = ci_upper,
    `95% CI Width` = ci_width,
    `95% CI Upper HW` = ci_upper_half_width,
    n, n_clus_raw, n_clus, data, seed) %>%
  tidyr::pivot_longer(cols = Estimate:`95% CI Upper HW`,
                      names_to = "statistic",
                      values_to = "values") %>%
  mutate(
    statistic = factor(statistic,
                       levels = c("Estimate", "Std. Error",
                                  "95% CI LL", "95% CI UL",
                                  "95% CI Width", "95% CI Upper HW")),
    n = factor(n, ordered = TRUE),
    data = stringr::str_to_title(data)
  )

glimpse(dist_plot_data)

dist_plot_nclus_5 <- dist_plot_data %>%
  filter(n_clus_raw == 5) %>%
  ggplot(aes(x = n, y = values, color = data)) +
  #geom_point(position = position_dodge(width = 0.2)) +
  geom_boxplot(position = position_dodge(width = 1)) +
  xlab("Total sample size (n)") +
  ylab("") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_wrap(~ statistic, ncol = 2, scales = "free")

plot(dist_plot_nclus_5)

plotly::ggplotly(dist_plot_nclus_5)

# test this with geem and gee packages as well ----
test_data <- simulate_data_and_get_estimates(n = 1000, n_clus = 5, seed = 10, return_data = TRUE)

glimpse(test_data)

geeglm(y ~ x, id = id, data = test_data[["sim_data"]], family = "binomial", corstr = "exch") %>%
  summary()

geeglm(y_prob ~ x, id = id, weights = n, data = test_data[["agg_data"]], family = "binomial", corstr = "exch") %>%
  summary()

# these estimates are very different from the geeglm ones, no clue why
geeM::geem(y ~ x, id = id, data = test_data[["sim_data"]], family = "binomial", corstr = "exch") %>%
  summary()

# with weights, throws data type errors, no clue why as all the data types are correct
test_data[["agg_data"]] %>%
  mutate(across(-id, ~ as.double(.x))) %>%
  geeM::geem(y_prob ~ x, weights = n, id = id,
             data = ., family = "binomial", corstr = "exchangeable") %>%
  summary()

gee_mod1 <- gee::gee(y ~ x, id = id, data = test_data[["sim_data"]], family = binomial, corstr = "exchangeable")

gee_mod1 # same as geeglm estimates

summary(gee_mod1)

# different from the raw data estimates from gee::gee as well as geepack::geeglm on aggregated data
gee::gee(cbind(y, n-y) ~ x, id = id, data = test_data[["agg_data"]], family = binomial, corstr = "exchangeable") %>% summary()

