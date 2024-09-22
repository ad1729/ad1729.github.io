
library(tidyverse)

theme_set(theme_bw())

post_path <- fs::path(getwd(), "posts", "simulating-survival-data-with-non-PH-and-cure") %>%
  print()

# plot hazard functions for different parameters ----

hazard_function_loglogistic <- function(alpha, beta,
                                        min_t = 0,
                                        max_t = 24,
                                        step_t = 0.1) {

  time <- seq(min_t, max_t, step_t)
  numerator <- (beta / alpha) * ((time / alpha) ^ (beta - 1))
  denominator <- 1 + ((time / alpha) ^ beta)

  tibble(
    alpha = alpha,
    beta = beta,
    parameters = paste0("LL(", alpha,
                        ", ", beta, ")"),
    time = time,
    hazard = numerator / denominator
  )
}

hazard_function_loglogistic(alpha = 1, beta = 4, step_t = 1)

plot_hazard <- function(data, unit = "weeks") {
  data %>%
    ggplot(aes(x = time, y = hazard, group = parameters, color = parameters)) +
    geom_line() +
    labs(x = glue::glue("Time (in {unit}"), y = "Hazard") +
    theme(legend.title = element_blank(), legend.position = "bottom")
}

hazard_function_loglogistic(1, 4) %>% plot_hazard()

hazard_functions <- expand_grid(
  alpha = c(1, 2, 16),
  beta = c(1, 2, 16),
  time = seq(0, 25, 0.1)
) %>%
  mutate(
    parameters = paste0("a = ", alpha, ", ", "b = ", beta),
    numerator = (beta / alpha) * ((time / alpha) ^ (beta - 1)),
    denominator = 1 + ((time / alpha) ^ beta),
    hazard = numerator / denominator
  )

hazard_functions %>%
  ggplot(aes(x = time, y = hazard, group = parameters, color = parameters)) +
  geom_line() +
  labs(x = "Time (in weeks)", y = "Hazard") +
  #theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(legend.position = "none") +
  #facet_grid(alpha ~ beta, scales = "free_y")
  facet_wrap(~ parameters, scales = "free_y")

hazard_functions %>%
  filter(alpha %in% c(1, 2), beta %in% c(2)) %>%
  ggplot(aes(x = time, y = hazard, group = parameters, color = parameters)) +
  geom_line() +
  labs(x = "Time (in weeks)", y = "Hazard") +
  theme(legend.title = element_blank(), legend.position = "bottom")

# beta = 7, alpha = 4 and 5
bind_rows(
  hazard_function_loglogistic(alpha = 4, beta = 7, max_t = 240),
  hazard_function_loglogistic(alpha = 5, beta = 7, max_t = 240)
) %>%
  plot_hazard() +
  geom_vline(xintercept = 6, linetype = "dotted", color = "gray60", linewidth = 1.3) +
  xlab("Time (in weeks)")

# set beta to 7, alpha to 4 for treatment and 5 for control arms ----
# code below in the post

# evaluate and visualize the arm specific hazards
loglogistic_attributes <- bind_rows(
  hazard_function_loglogistic(alpha = 4, beta = 7),
  hazard_function_loglogistic(alpha = 5, beta = 7)
) %>%
  mutate(
    parameters = case_when(
      alpha == 4 ~ "Treatment group (T ~ LL(4, 7))",
      alpha == 5 ~ "Control group (T ~ LL(5, 7))",
    )
  )

loglogistic_attributes %>%
  ggplot(aes(x = time, y = hazard, group = parameters, color = parameters)) +
  geom_line(linewidth = 1.1) +
  # campaign end date
  geom_vline(xintercept = 6, linetype = "dotted", color = "gray40") +
  annotate(geom = "text", x = 6.2, y = 0.6, color = "gray40",
           label = "Campaign\nends after\n6 weeks", hjust = "left") +
  # impact of treatment disappears by 2.5 months
  geom_vline(xintercept = 10, linetype = "dotted", color = "gray40") +
  annotate(geom = "text", x = 9.8, y = 0.1, color = "gray40",
           label = "Impact of\ntreatment\ndies down", hjust = "right") +
  # analysis date (right-censoring at 3 months)
  geom_vline(xintercept = 12, linetype = "dotted", color = "gray40") +
  annotate(geom = "text", x = 12.2, y = 0.1, color = "gray40",
           label = glue::glue("Cut-off date for analysis\n",
                              "(administrative right-censoring ",
                              "at 3 months)"),
           hjust = "left") +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  scale_y_continuous(breaks = seq(0, 1.2, 0.2)) +
  labs(x = "Time since product launch (in weeks)", y = "Hazard rate") +
  theme(legend.position = "bottom", legend.title = element_blank())

# survivor function and its inverse ----
survivor_loglogistic <- function(alpha, beta, time) {
  1 / (1 + ((time / alpha) ^ beta))
}

loglogistic_attributes <- loglogistic_attributes %>%
  mutate(
    survival = survivor_loglogistic(alpha, beta, time),
    cdf = 1 - survival
  )

# visualize the CICs
cumulative_incidence_plot <- loglogistic_attributes %>%
  ggplot(aes(x = time, y = cdf, group = parameters, color = parameters)) +
  geom_line(linewidth = 1.1) +
  # campaign end date
  geom_vline(xintercept = 6, linetype = "dotted", color = "gray40") +
  annotate(geom = "text", x = 6.2, y = 0.6, color = "gray40",
           label = "Campaign\nends after\n6 weeks", hjust = "left") +
  # impact of treatment disappears by 2.5 months
  geom_vline(xintercept = 10, linetype = "dotted", color = "gray40") +
  annotate(geom = "text", x = 9.8, y = 0.1, color = "gray40",
           label = "Impact of\ntreatment\ndies down", hjust = "right") +
  # analysis date (right-censoring at 3 months)
  geom_vline(xintercept = 12, linetype = "dotted", color = "gray40") +
  annotate(geom = "text", x = 12.2, y = 0.1, color = "gray40",
           label = glue::glue("Cut-off date for analysis\n",
                              "(administrative right-censoring ",
                              "at 3 months)"),
           hjust = "left") +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Time since product launch (in weeks)",
       y = "Cumulative incidence") +
  theme(legend.position = "bottom", legend.title = element_blank())

cumulative_incidence_plot +
  scale_x_continuous(breaks = seq(0, 24, 4))

# simulate event times ----
n_sim <- 30
n_sample <- 1000
n_total <- n_sim * n_sample

set.seed(43)
loglogistic_samples <- tibble(
  sim_id = rep(1:n_sim, each = n_sample),
  id = rep(1:n_sample, times = n_sim),
  treatment = sample(c(0, 1), size = n_total,
                     replace = TRUE, prob = c(0.5, 0.5)),
  parameters = case_when(
    treatment == 1 ~ "Treatment group (T ~ LL(4, 7))",
    treatment == 0 ~ "Control group (T ~ LL(5, 7))",
  ),
  # -log(1.25) corresponds to alpha = 5 in control
  # and alpha = 4 in treatment
  linear_predictor_alpha = log(5) - log(1.25) * treatment,
  linear_predictor_beta =log(7),
  u = runif(n_total, min = 0, max = 1),
  time = (exp(linear_predictor_alpha) *
            (((1 / u) - 1) ^ (-1 / exp(linear_predictor_beta))))
)

loglogistic_samples %>% glimpse()

# plot randomly sampled cumulative incidence curves for each arm ----
cumulative_incidence_plot_samples <- cumulative_incidence_plot

# for this hacky solution to overlay the true curves over
# the sampled curves, see
# https://stackoverflow.com/questions/20249653/insert-layer-underneath-existing-layers-in-ggplot2-object
cumulative_incidence_plot_samples$layers <- c(
  stat_ecdf(
    data = loglogistic_samples,
    aes(x = time, linetype = parameters,
        group = interaction(sim_id, parameters)),
    inherit.aes = FALSE,
    color = "gray70"
  ),
  cumulative_incidence_plot$layers
)

# in a single plot
cumulative_incidence_plot_samples +
  coord_cartesian(xlim = c(2, 10)) +
  scale_x_continuous(breaks = seq(2, 10, 1)) +
  guides(
    color = guide_legend(nrow = 2),
    linetype = guide_legend(nrow = 2)
  )

# save this plot as the main plot from the post
ggsave(plot = last_plot(), filename = "post-image.png",
       device = "png", units = "px", width = 900, height = 800, dpi = 150,
       path = post_path)

# facet by arm
cumulative_incidence_plot_samples + facet_wrap(~parameters)

# add censoring and cure ----

# administrative censoring
loglogistic_samples %>%
  select(time) %>%
  mutate(
    latent = time,
    time3 = pmin(time, 3.0),
    cens3 = as.numeric(time > 3.0),
    time6 = pmin(time, 6.0),
    cens6 = as.numeric(time > 6.0),
    time12 = pmin(time, 12.0),
    cens12 = as.numeric(time > 12.0)
  ) %>%
  select(-time) %>%
  summary(digits = 2)

# overlay the event time distributions and random censoring
# times samples from Gamma(6, 1)
set.seed(43)
loglogistic_samples %>%
  select(parameters, time) %>%
  bind_rows(
    tibble(
      parameters = "Censoring times ~ Gamma(6, 1)",
      time = rgamma(floor(n_total / 2), shape = 6, scale = 1)
    )
  ) %>%
  ggplot(aes(x = time, fill = parameters, color = parameters)) +
  geom_density(alpha = 0.6) +
  #stat_ecdf(aes(color = parameters), geom = "line", linewidth = 1.3) +
  xlab("Time since product launch (in weeks)") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = c("purple3", "gray20", "forestgreen")) +
  scale_fill_manual(values = c("purple3", "gray20", "forestgreen"))

# simulate censoring time randomly for each individual and summarize
set.seed(43)
loglogistic_samples %>%
  select(parameters, latent = time) %>%
  mutate(
    censoring_time = rgamma(n_total, shape = 6, scale = 1),
    observed_time = pmin(latent, censoring_time),
    cens = as.numeric(latent > censoring_time)
  ) %>%
  summarize(
    across(
      .cols = everything(),
      .fns = ~ round(mean(.x), 2),
      .names = "mean_{.col}"
    ),
    .by = parameters
  )

# partition the overall sample into uncensored, censored and uncured,
# and censored and cured with fixed overall cure fraction of 30%
set.seed(43)
loglogistic_samples %>%
  mutate(
    cured = rbinom(n_total, 1, 0.3),
    time = ifelse(cured == 1, 10000, time),
    time6 = pmin(time, 6.0),
    cens6 = as.numeric(time > 6.0)
  ) %>%
  count(cured, cens6) %>%
  mutate(p = round(100 * (n / sum(n)), 1))

# adding a covariate that influences cure fraction
plogis(4 - 0.07 * seq(20, 90, 10)) %>%
  purrr::set_names(nm = as.character(seq(20, 90, 10)))

rgamma(n = n_total, shape = 40, rate = 0.7) %>%
  pmin(., 90) %>% pmax(20, .) %>% round(., 1) %>%
  summary()

loglogistic_samples %>%
  mutate(
    age = {
      rgamma(n = n_total, shape = 40, rate = 0.7) %>%
        pmin(., 90) %>%
        pmax(20, .) %>%
        round(., 1)
    },
    cured = rbinom(n_total, 1, plogis(4 - 0.07 * age)),
    time = ifelse(cured == 1, 10000, time),
    time6 = pmin(time, 6.0),
    cens6 = as.numeric(time > 6.0),
    event6 = 1 - cens6,
  )

# simulate data for the CIC comparisons ----
set.seed(43)
loglogistic_samples %>%
  mutate(
    # add cure status and modify time distribution accordingly
    latent = time,
    cured = rbinom(n_total, 1, 0.3),
    time = ifelse(cured == 1, 10000, time),
    # add censoring indicator
    time3 = pmin(time, 3.0),
    cens3 = as.numeric(time > 3.0),
    event3 = 1 - cens3,
    time6 = pmin(time, 6.0),
    cens6 = as.numeric(time > 6.0),
    event6 = 1 - cens6,
    time12 = pmin(time, 12.0),
    cens12 = as.numeric(time > 12.0),
    event12 = 1 - cens12,
    across(starts_with("event"), as.integer)
  ) %>%
  #slice_sample(n = 1000) %>%
  summarize(
    across(.cols = c(latent, time, starts_with("cens"), cured),
           .fns = mean),
    .by = parameters
  )
# very weird behaviour when seed is set to 43 here
# cure fraction ends up being 0 in the treatment group
# and 60% in the control group. Overall it's correct with value of 0.3
# this goes away by removing the set.seed(43) call before generating
# or using any other seed

set.seed(49)
censored_and_cured_samples <- loglogistic_samples %>%
  mutate(
    # add cure status and modify time distribution accordingly
    latent = time,
    cured = rbinom(n_total, 1, 0.3),
    time = ifelse(cured == 1, 10000, time),
    # add censoring indicator
    time3 = pmin(time, 3.0),
    cens3 = as.numeric(time > 3.0),
    event3 = 1 - cens3,
    time6 = pmin(time, 6.0),
    cens6 = as.numeric(time > 6.0),
    event6 = 1 - cens6,
    time12 = pmin(time, 12.0),
    cens12 = as.numeric(time > 12.0),
    event12 = 1 - cens12,
    across(starts_with("event"), as.integer)
  )

censored_and_cured_samples %>%
  select(latent, time, starts_with("cens")) %>%
  summary()

censored_and_cured_samples %>%
  #slice_sample(n = 1000) %>%
  summarize(
    across(.cols = c(latent, time, starts_with("cens"), cured),
           .fns = mean),
    .by = parameters
  )

twelve_weeks_data <- censored_and_cured_samples %>%
  select(group = parameters, time = time12, event = event12) %>%
  print(n = 5)

twelve_weeks_data %>%
  write_rds(
    file = fs::path(post_path, "twelve-weeks-data.rds")
  )

three_weeks_data <- censored_and_cured_samples %>%
  select(group = parameters, time = time3, event = event3) %>%
  print(n = 5)

three_weeks_data %>%
  write_rds(
    file = fs::path(post_path, "three-weeks-data.rds")
  )

loglogistic_attributes %>%
  write_rds(
    file = fs::path(post_path, "loglogistic-attributes.rds")
  )
