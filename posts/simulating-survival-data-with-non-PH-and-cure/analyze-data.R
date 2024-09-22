
library(tidyverse)
# library(smcure)
# library(flexsurv)
library(casebase)
library(flexsurvcure)
library(ggsurvfit)
library(mgcv)

theme_set(theme_bw())

post_path <- fs::path(getwd(), "posts", "simulating-survival-data-with-non-PH-and-cure") %>%
  print()

# datasets are generated in simulate-data.R file
three_weeks_data <- read_rds(fs::path(post_path, "three-weeks-data.rds")) %>% print()
twelve_weeks_data <- read_rds(fs::path(post_path, "twelve-weeks-data.rds")) %>% print()
loglogistic_attributes <- read_rds(fs::path(post_path, "loglogistic-attributes.rds")) %>% print()

# log-log survival plot
log_log_survival_plot <- twelve_weeks_data %>%
  survfit2(Surv(time, event, type = "right") ~ group, data = .) %>%
  ggsurvfit(type = "cloglog") +
  scale_color_manual(values = c("gray20", "forestgreen"))

log_log_survival_plot +
  scale_x_continuous(breaks = seq(0, 12, 1)) +
  xlab("Time since product launch (in weeks)")

log_log_survival_plot +
  scale_ggsurvfit(
    x_scales = list(transform = "log"),
    y_scales = list(label = scales::label_percent(scale = 1, suffix = ""))
  ) +
  xlab("Time since product launch (in weeks, log scale)")

# proportional odds plot
# use ggsurvfit
twelve_weeks_data %>%
  survfit2(Surv(time, event, type = "right") ~ group, data = .) %>%
  ggsurvfit(type = "survival") +
  scale_ggsurvfit(
    x_scales = list(transform = "log"),
    y_scales = list(transform = "logit")
  ) +
  scale_color_manual(values = c("gray20", "forestgreen"))

# make it manually
proportional_odds_plot <- twelve_weeks_data %>%
  survfit2(Surv(time, event, type = "right") ~ group, data = .) %>%
  tidy_survfit(time = seq(0, 12, 0.05), type = "survival") %>%
  select(time, strata, estimate) %>%
  mutate(
    lntime = log(time),
    prop_odds = log(estimate / (1 - estimate))
  ) %>%
  ggplot(aes(x = lntime, y = prop_odds, color = strata, group = strata)) +
  geom_step() +
  labs(x = "Log time", y = "Log survival odds") +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  coord_cartesian(xlim = c(-0.5, 2.5)) +
  theme(legend.position = "bottom", legend.title = element_blank())

proportional_odds_plot

gridExtra::grid.arrange(log_log_survival_plot, proportional_odds_plot)

# can also make it as a facetted plot
twelve_weeks_data %>%
  survfit2(Surv(time, event, type = "right") ~ group, data = .) %>%
  tidy_survfit(time = seq(0, 12, 0.05), type = "survival") %>%
  select(time, strata, estimate) %>%
  mutate(
    lntime = log(time),
    `Log Minus Log Survival` = log(-log(estimate)),
    `Log Survival Odds` = log(estimate / (1 - estimate))
  ) %>%
  pivot_longer(cols = starts_with("Log ")) %>%
  ggplot(aes(x = lntime, y = value, color = strata, group = interaction(name, strata))) +
  geom_step() +
  labs(x = "Log time", y = "") +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  coord_cartesian(xlim = c(-0.5, 2.5)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  facet_wrap(~name, scales = "free_y")

# smcure code, gives errors for AFT model
# cure_model <- censored_and_cured_samples %>%
#   select(group = parameters, time12, event12) %>%
#   smcure(
#     formula = formula(Surv(time12, event12) ~ group),
#     cureform = ~ group,
#     data = .,
#     model = "aft", # gives nonsensical error with "aft", works with "ph"
#     link = "logit",
#     Var = FALSE,
#     emmax = 10
#   )

fit_loglogistic_cure_model <- function(data) {
  data %>%
    flexsurvcure(
      formula = Surv(time, event, type = "right") ~ group,
      data = .,
      dist = "llogis",
      link = "logistic",
      # treatment doesn't impact shape parameter beta, only scale parameter alpha,
      # but let's add that to the model too
      # covariates on the survival distribution parameters go here
      # coefficients are interpretable on log scale parameter alpha
      # and log shape parameter beta
      anc = list(
        scale = ~ group,
        shape = ~ group
      ),
      mixture = TRUE
    )
}

cure_model_twelve_weeks <- fit_loglogistic_cure_model(twelve_weeks_data)

cure_model_twelve_weeks

# get the non-exponentiated coefficients
cure_model_twelve_weeks$coefficients

# check the output data
cure_model_twelve_weeks %>%
  summary(type = "survival", t = seq(0, 12, 1 / 7), ci = TRUE, tidy = TRUE) %>%
  group_by(group) %>%
  slice(seq(10, 50, 10)) %>%
  arrange(time)

get_cumulative_incidence_curve_from_cure_model <- function(cure_model, max_t = 12) {
  cure_model %>%
    summary(type = "survival", t = seq(0, max_t, 1 / 7), ci = TRUE, tidy = TRUE) %>%
    mutate(
      CI_est = 1 - est,
      CI_lcl = 1 - ucl,
      CI_ucl = 1 - lcl,
      type = "Log-logistic mixture cure model"
    ) %>%
    select(type, group, time, starts_with("CI_")) %>%
    as_tibble()
}

cure_model_cumulative_incidence <- cure_model_twelve_weeks %>%
  get_cumulative_incidence_curve_from_cure_model() %>%
  print()

# this is the CDF for the uncured individuals, so it's defined on [0, 1]
# for some proportion of cured individuals, we would need to rescale these to [0, 1 - cure fraction]
# because for any time axis, the CDF will be cut off at the cure fraction
# so rescale [0, 1] to [0, 0.7]
# which corresponds to a cure fraction of 0.3 in each group
# for the simple scaling formula, see
# https://stats.stackexchange.com/questions/281162/scale-a-number-between-a-range
true_cumulative_incidence <- loglogistic_attributes %>%
  mutate(CI_est = cdf * 0.7, type = "Truth") %>%
  select(type, group = parameters, time, CI_est) %>%
  print()

# fit KM model, cox model, casebase package model (splines for hazard function)
# and get the survival curve to derive the CI curve (or get CI cruve directly)
fit_KM_and_get_cumulative_incidence <- function(data, max_t = 12) {
  data %>%
    survfit2(Surv(time, event, type = "right") ~ group, data = .) %>%
    tidy_survfit(times = seq(0, max_t, 1 / 7), type = "risk") %>%
    mutate(group = as.character(strata), type = "Kaplan-Meier fit") %>%
    select(type, group, time, CI_est = estimate,
           CI_lcl = conf.low, CI_ucl = conf.high)
}

km_fit_cumulative_incidence <- twelve_weeks_data %>%
  fit_KM_and_get_cumulative_incidence() %>%
  print(n = 10)

fit_stratified_cox_and_get_cumulative_incidence <- function(data, max_t = 12) {
  data %>%
    coxph(Surv(time, event, type = "right") ~ strata(group), data = .) %>%
    survfit2() %>%
    tidy_survfit(times = seq(0, max_t, 1 / 7), type = "risk") %>%
    mutate(group = as.character(strata), type = "Stratified Cox model") %>%
    select(type, group, time, CI_est = estimate,
           CI_lcl = conf.low, CI_ucl = conf.high)
}

cox_fit_cumulative_incidence <- twelve_weeks_data %>%
  fit_stratified_cox_and_get_cumulative_incidence() %>%
  glimpse()

# sanity check by getting group-wise five number summaries
cox_fit_cumulative_incidence %>%
  group_split(group) %>%
  map(.f = ~ summary(.x))

# fit separate functions of time for each arm
casebase_fit <- twelve_weeks_data %>%
  # need to convert to a factor to avoid errors
  mutate(group = factor(group)) %>%
  fitSmoothHazard(event ~ group + s(time, by = group),
                  data = ., family = "gam",
                  time = "time", ratio = 10)

summary(casebase_fit)

# for exploring the functions to get CI curves and uncertainty intervals
twelve_weeks_data %>%
  distinct(group) %>%
  mutate(group = factor(group)) %>%
  absoluteRisk(
    object = casebase_fit,
    time = seq(0, 12, 1 / 7),
    newdata = .
  ) %>%
  confint(parm = casebase_fit, nboot = 5) %>%
  as_tibble()

# when no confidence intervals are requested
# casebase_fit_cumulative_incidence <- twelve_weeks_data %>%
#   distinct(group) %>%
#   mutate(group = factor(group)) %>%
#   absoluteRisk(
#     object = casebase_fit,
#     time = seq(0, 12, 1 / 7),
#     newdata = .
#   ) %>%
#   as_tibble() %>%
#   pivot_longer(-time) %>%
#   mutate(
#     group = case_when(
#       name == "V2" ~ "Control group (T ~ LL(5, 7))",
#       name == "V3" ~ "Treatment group (T ~ LL(4, 7))"
#     ),
#     across(.cols = c(time, value), .fns = ~ as.numeric(.x)),
#     type = "Smooth-in-Time model"
#   ) %>%
#   select(type, group, time, CI = value) %>%
#   print(n = 10)

# try on a small subset of the data
twelve_weeks_data %>%
  mutate(group = factor(group)) %>%
  slice_sample(n = 5000) %>%
  fitSmoothHazard(
    event ~ group + s(time, by = group),
    data = ., family = "gam",
    time = "time", ratio = 10
  ) %>%
  summary()

# compare with running this on the full data
casebase_fit %>% summary()

# wrap it up in a function
fit_casebase_model_and_get_cumulative_incidence <- function(data, nboot = 20, max_t = 12) {
  data <- data %>%
    # need to have this as a factor for it to work
    # correctly with mgcv for factor smooth interactions
    mutate(group = factor(group))

  cat("\nFitting model...")

  # can take a while to fit (depends on ratio, higher takes longer)
  model <- data %>%
    fitSmoothHazard(event ~ group + s(time, by = group),
                    data = ., family = "gam",
                    time = "time", ratio = 10)

  arms <- data %>% distinct(group)

  cat("\nEstimating confidence intervals...")

  curves <- arms %>%
    absoluteRisk(object = model, time = seq(0, max_t, 1 / 7), newdata = ., type = "CI") %>%
    confint(object = ., parm = model, nboot = nboot) %>%
    as_tibble() %>%
    mutate(type = "Smooth-in-Time model",
           group = as.character(group)) %>%
    select(type, group, time, CI_est = estimate,
           CI_lcl = conf.low, CI_ucl = conf.high)

  return(curves)
}

# test function
twelve_weeks_data %>%
  slice_sample(n = 5000) %>%
  fit_casebase_model_and_get_cumulative_incidence(nboot = 5)

# fit one model with large nboot (takes a long time)
# the saved results can be read back in
casebase_fit_cumulative_incidence <- twelve_weeks_data %>%
  fit_casebase_model_and_get_cumulative_incidence(nboot = 100) %>%
  print(n = 10)
#
# casebase_fit_cumulative_incidence %>%
#   write_rds(file = fs::path(post_path, "casebase-12weeks.rds"))
#
# casebase_fit_cumulative_incidence <- read_rds(file = fs::path(post_path, "casebase-12weeks.rds")) %>%
#   glimpse()

plot_data <- bind_rows(
  true_cumulative_incidence,
  cure_model_cumulative_incidence,
  cox_fit_cumulative_incidence,
  km_fit_cumulative_incidence,
  casebase_fit_cumulative_incidence
) %>%
  mutate(
    type = factor(type, levels = c(
      "Truth", "Log-logistic mixture cure model", "Stratified Cox model",
      "Kaplan-Meier fit", "Smooth-in-Time model"
    ))
  ) %>%
  glimpse()

plot_data %>%
  write_rds(file = fs::path(post_path, "plot-data-12-weeks.rds"))

# plot the true cure CDF vs the predicted cure CDF
# for the time period observed
# (extrapolation is not the goal)
CI_curve_plots <- plot_data %>%
  ggplot(aes(x = time, y = CI_est, color = group)) +
  geom_line(aes(linetype = type, linewidth = type)) +
  # add cure fraction line and text
  geom_hline(yintercept = 0.7, color = "gray40", linetype = "solid") +
  annotate(geom = "text", x = 0.5, y = 0.72, color = "gray40",
           label = "Cure fraction of 70% in both groups", hjust = "left") +
  labs(x = "Time since product launch (in weeks)",
       y = "% of conversions") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed", "dotted", "longdash")) +
  scale_linewidth_manual(values = c(0.4, 1.1, 1.1, 1.1, 1.1)) +
  scale_x_continuous(breaks = seq(0, 12, 1)) +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 1)) +
  guides(
    color = guide_legend(nrow = 2),
    linetype = guide_legend(nrow = 2)
  )

CI_curve_plots
CI_curve_plots + facet_wrap(~ group)
CI_curve_plots + facet_wrap(~ type)

# plot the CI curves with the confidence bands
CI_curve_plots +
  # add pointwise confidence bands
  geom_ribbon(aes(ymin = CI_lcl, ymax = CI_ucl, fill = group,
                  linetype = type, linewidth = type), alpha = 0.2) +
  scale_fill_manual(values = c("gray20", "forestgreen"))

# same as previous but truth plotted in each facet separately but facetting plots by type
plot_data %>%
  filter(type != "Truth") %>%
  mutate(type = fct_drop(type)) %>%
  ggplot(aes(x = time, y = CI_est, color = group)) +
  geom_line(linetype = "dashed", linewidth = 1.1) +
  geom_line(data = select(filter(plot_data, type == "Truth"), -type),
            inherit.aes = FALSE,
            aes(x = time, y = CI_est, color = group),
            linetype = "solid", linewidth = 0.4) +
  # add pointwise confidence bands
  geom_ribbon(aes(ymin = CI_lcl, ymax = CI_ucl, fill = group), linetype = "dashed", alpha = 0.2) +
  # add cure fraction line and text
  geom_hline(yintercept = 0.7, color = "gray40", linetype = "solid") +
  annotate(geom = "text", x = 0.5, y = 0.72, color = "gray40",
           label = "Cure fraction of 70% in both groups", hjust = "left") +
  labs(x = "Time since product launch (in weeks)",
       y = "% of conversions") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  scale_fill_manual(values = c("gray20", "forestgreen")) +
  scale_x_continuous(breaks = seq(0, 12, 1)) +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 1)) +
  facet_wrap(~ type)

# plot by arm but compare estimates and intervals facet by group
# plot_data %>%
#   filter(type != "Truth") %>%
#   ggplot(aes(x = time, y = CI_est, group = type)) +
#   geom_line(aes(color = type)) +
#   geom_ribbon(aes(ymin = CI_lcl, ymax = CI_ucl, fill = type), alpha = 0.2) +
#   geom_line(data = true_cumulative_incidence, inherit.aes = TRUE,
#             color = "black", linewidth = 1.1) +
#   # add cure fraction line and text
#   geom_hline(yintercept = 0.7, color = "gray40", linetype = "solid") +
#   annotate(geom = "text", x = 0.5, y = 0.72, color = "gray40",
#            label = "Cure fraction of 70% in both groups", hjust = "left") +
#   labs(x = "Time since product launch (in weeks)",
#        y = "% of conversions") +
#   theme(legend.position = "bottom", legend.title = element_blank()) +
#   facet_wrap(~ group) +
#   scale_x_continuous(breaks = seq(0, 12, 1)) +
#   scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.2)) +
#   coord_cartesian(xlim = c(0, 12), ylim = c(0, 1))

# finally, the same AFT mixture cure model is fit after 3 weeks instead of 12
# kaplan-meier fit
three_weeks_data %>%
  survfit2(Surv(time, event, type = "right") ~ group, data = .) %>%
  ggsurvfit() +
  add_censor_mark() +
  coord_cartesian(ylim = c(0, 1))

cure_model_three_weeks <- fit_loglogistic_cure_model(three_weeks_data)

cure_model_three_weeks

cure_model_three_weeks %>%
  write_rds(file = fs::path(post_path, "cure-model-3-weeks.rds"))

interim_coefs <- cure_model_three_weeks$coefficients

interim_coefs

plogis(c(interim_coefs[1], interim_coefs[1] + interim_coefs[4]))

# very large SEs, also these are different from the ones from printing cure_model_three_weeks
# don't know why, I'm probably doing something wrong
cure_model_three_weeks %>%
  vcov() %>% diag() %>% sqrt() %>%
  set_names(
    nm = ~ str_replace(
      string = .x,
      pattern = "groupTreatment group \\(T ~ LL\\(4, 7\\)\\)",
      replacement = "trt"
    )
  ) # compare with last column in the next df
cure_model_three_weeks %>% pluck("res")

plot_data_three_weeks <- bind_rows(
  true_cumulative_incidence,
  get_cumulative_incidence_curve_from_cure_model(cure_model = cure_model_three_weeks, max_t = 3),
  fit_stratified_cox_and_get_cumulative_incidence(data = three_weeks_data, max_t = 3),
  fit_KM_and_get_cumulative_incidence(data = three_weeks_data, max_t = 3),
  fit_casebase_model_and_get_cumulative_incidence(data = three_weeks_data, nboot = 100, max_t = 3)
) %>%
  mutate(
    type = factor(type, levels = c(
      "Truth", "Log-logistic mixture cure model", "Stratified Cox model",
      "Kaplan-Meier fit", "Smooth-in-Time model"
    ))
  ) %>%
  glimpse()

plot_data_three_weeks

plot_data_three_weeks %>%
  write_rds(file = fs::path(post_path, "plot-data-3-weeks.rds"))

# separate method in each facet
three_weeks_plot_facetted <- plot_data_three_weeks %>%
  filter(type != "Truth") %>%
  mutate(type = fct_drop(type)) %>%
  ggplot(aes(x = time, y = CI_est, color = group)) +
  geom_line(linetype = "dashed", linewidth = 1.1) +
  geom_line(data = {
    plot_data_three_weeks %>%
      filter(type == "Truth", time <= 3) %>%
      select(-type)
  },
  inherit.aes = FALSE,
  aes(x = time, y = CI_est, color = group),
  linetype = "solid", linewidth = 0.4) +
  # add cure fraction line and text
  # geom_hline(yintercept = 0.7, color = "gray40", linetype = "solid") +
  # annotate(geom = "text", x = 0.5, y = 0.72, color = "gray40",
  #          label = "Cure fraction of 70% in both groups", hjust = "left") +
  labs(x = "Time since product launch (in weeks)",
       y = "% of conversions") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  scale_x_continuous(breaks = seq(0, 3, 1)) +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 0.1, 0.02)) +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 0.1)) +
  facet_wrap(~ type)

three_weeks_plot_facetted

three_weeks_plot_facetted + coord_cartesian(xlim = c(2, 3), ylim = c(0, 0.1))

# all in one facet
three_weeks_plot_single <- plot_data_three_weeks %>%
  filter(time <= 3) %>%
  ggplot(aes(x = time, y = CI_est, color = group)) +
  geom_line(aes(linetype = type, linewidth = type)) +
  labs(x = "Time since product launch (in weeks)",
       y = "% of conversions") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed", "dotted", "longdash")) +
  scale_linewidth_manual(values = c(0.4, 1.1, 1.1, 1.1, 1.1)) +
  scale_x_continuous(breaks = seq(0, 3, 1)) +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 0.1, 0.02)) +
  guides(
    color = guide_legend(nrow = 2),
    linetype = guide_legend(nrow = 2)
  )

three_weeks_plot_single

three_weeks_plot_single + coord_cartesian(xlim = c(0, 3), ylim = c(0, 0.1))

three_weeks_plot_single + coord_cartesian(xlim = c(2, 3), ylim = c(0, 0.1))

# map type to color and group to linetype
# plot_data_three_weeks %>%
#   filter(time <= 3) %>%
#   ggplot(aes(x = time, y = CI_est)) +
#   geom_line(aes(color = type, linetype = group), linewidth = 1.1) +
#   labs(x = "Time since product launch (in weeks)",
#        y = "% of conversions") +
#   theme(legend.position = "bottom", legend.title = element_blank()) +
#   scale_linetype_manual(values = c("dotdash", "solid")) +
#   scale_x_continuous(breaks = seq(0, 3, 1)) +
#   scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 0.1, 0.02)) +
#   guides(
#     color = guide_legend(nrow = 2),
#     linetype = guide_legend(nrow = 2)
#   ) +
#   coord_cartesian(xlim = c(2, 3), ylim = c(0, 0.1))

# same as before but facet by type and plot CIs
plot_data_three_weeks %>%
  filter(time <= 3, type != "Truth") %>%
  ggplot(aes(x = time, y = CI_est)) +
  # plot the true line on each panel
  geom_line(
    data = {
      plot_data_three_weeks %>%
        filter(type == "Truth", time <= 3) %>%
        select(-type)
    },
    inherit.aes = FALSE,
    aes(x = time, y = CI_est, group = group, color = group),
    linetype = "solid", linewidth = 1
  ) +
  # plot estimated curves and their confidence intervals
  geom_line(aes(color = group), linetype = "dashed", linewidth = 1.1) +
  geom_ribbon(aes(ymin = CI_lcl, ymax = CI_ucl, color = group, fill = group), alpha = 0.3) +
  labs(x = "Time since product launch (in weeks)",
       y = "% of conversions") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_color_manual(values = c("gray20", "forestgreen")) +
  scale_fill_manual(values = c("gray20", "forestgreen")) +
  scale_x_continuous(breaks = seq(2, 3, 0.2)) +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 0.1, 0.02)) +
  facet_wrap(~type) +
  coord_cartesian(xlim = c(2, 3), ylim = c(0, 0.1))
