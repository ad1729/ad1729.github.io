
# extra code for checking that the posterior predictive distribution
# is unbiased for estimates from different simulated datasets
# when the true parameters are used instead of estimated parameters from a
# single dataset
#
# this file can be run standalone, and has no dependencies on the main index.qmd file
#
# the parameters here 0.0882, 0.0082, 0.8, 0.001
# are the true probabilities for the main dataset used

library(tidyverse)

simulate_data <- function(n = 100000L,
                          pY0 = 0.01,
                          risk_diff = 0.1,
                          seed = 23,
                          p_compliers = 0.8,
                          p_treatment = 0.7,
                          pY_non_compliers_factor = 0.1) {

  pY1 <- pY0 + risk_diff
  pY_nc <- pY_non_compliers_factor * pY0

  set.seed(seed)
  data <- tibble(
    id = 1:n,
    # the underlying population can be stratified into
    # never-takers and compliers
    complier = rbinom(n, 1, prob = p_compliers),
    # generate individual potential outcomes
    # under control and treatment, i.e., Pr[Y^0 = 1]
    #Y0 = rbinom(n, 1, prob = pY0),
    Y0 = case_when(
      complier == 0 ~ rbinom(n, 1, pY_nc),
      complier == 1 ~ rbinom(n, 1, pY0),
    ),
    # assuming a constant effect of +10 percentage points
    # among the compliers, and no average effect under the never-takers
    Y1 = case_when(
      complier == 0 ~ rbinom(n, 1, pY_nc),
      complier == 1 ~ rbinom(n, 1, pY1)
    ),
    # treatment assigned at random
    # 70-30 split into treatment / control
    Z = rbinom(n, 1, prob = p_treatment),
    # treatment uptake depends on
    # being assigned to treatment (Z = 1)
    # AND being a complier (C = 1)
    A = Z * complier,
    # generate observed response using the
    # consistency equation
    Y = (1 - Z) * Y0 + Z * Y1
  )

  return(data)
}

profit_distributions_sensitivity <- expand_grid(
  id = 1:10000,
  compliance = c(0.5, 0.8, 1.0),
  split = c(0.7, 1.0, 0.0),
  # rescale the ITT probabilities by p(compliance)
  # shave off the effect of never takers on this probability
  # to recover the p(Y = 1 | Z = z, C = 1)
  p_treated = (0.0882 / 0.8) - (0.001 * (0.2 / 0.8)),
  p_control = (0.0082 / 0.8) - (0.001 * (0.2 / 0.8)),
  n = 100000
) %>%
  mutate(
    # rescale the IV estimates for conversion prob
    # under each treatment by the compliance
    p_treated = ((compliance * p_treated) +
                   ((1 - compliance) * (0.001 / 0.8))),
    p_control = ((compliance * p_control) +
                   ((1 - compliance) * (0.001 / 0.8))),
    n_treated = round(split * n),
    n_control = n - n_treated,
    s_treated = rbinom(n(), n_treated, p_treated),
    s_control = rbinom(n(), n_control, p_control),
    total_profit = 100 * (s_treated + s_control)
  ) %>%
  glimpse()

simulated_profit_sensitivity <- profit_distributions_sensitivity %>%
  distinct(compliance, split) %>%
  expand_grid(id = seq(from = 20, by = 1, length.out = 50)) %>%
  pmap_dfr(
    .f = function(compliance, split, id) {
      simulate_data(p_compliers = compliance,
                    p_treatment = split,
                    seed = id) %>%
        summarise(total_profit = 100 * sum(Y)) %>%
        tibble(compliance, split, .)
    }
  ) %>%
  glimpse()

profit_distributions_plot_data_sensitivity <- list(
  "extrapolated" = profit_distributions_sensitivity %>%
    select(compliance, split, total_profit),
  "simulated" = simulated_profit_sensitivity
) %>%
  map(
    .f = ~ {
      .x %>%
        mutate(
          total_profit = total_profit / 1e5,
          compliance = paste0("Compliers: ", compliance * 100, "%"),
          compliance = factor(compliance,
                              levels = c("Compliers: 50%",
                                         "Compliers: 80%",
                                         "Compliers: 100%")),
          split = paste0("Treatment: ", split * 100, "%"),
          split = factor(split,
                         levels = c("Treatment: 0%",
                                    "Treatment: 70%",
                                    "Treatment: 100%"))
        )
    }
  ) %>%
  glimpse()

profit_distributions_plot_data_sensitivity %>%
  pluck("extrapolated") %>%
  ggplot(aes(x = total_profit)) +
  geom_density(trim = TRUE) +
  geom_vline(
    data = profit_distributions_plot_data_sensitivity %>% pluck("simulated"),
    aes(xintercept = total_profit), color = "maroon"
  ) +
  theme_bw() +
  facet_wrap(vars(split, compliance), scales = "free") +
  #facet_grid(compliance ~ split, scales = "fixed") +
  xlab("Total profit (multiples of 100,000 EUR, per 100,000 individuals)")
