
library(tidyverse)

n <- 1e6

# tweedie distribution parameters
p <- 1.3 # variance power
mu <- 200 # mean
phi <- 350 # dispersion, so variance = phi * mu ^ p

set.seed(45)
sim_data <- tibble(
  w = runif(n, min = 0.5, max = 3),
  y = tweedie::rtweedie(n = n, mu = mu, phi = phi / w, power = p)
)

# summary(sim_data)

sim_data %>% slice_head(n = 3)

parameter_grid <- seq(50, 300, 5)

estimating_equation_curves <- tibble(
  parameter = parameter_grid,
   `Mean of ratios` = map_dbl(.x = parameter_grid, .f = ~ {
    sim_data %>%
      mutate(s = w * (y - .x)) %>%
      pull(s) %>%
      sum()
  }),
  `Ratio of means` = map_dbl(.x = parameter_grid, .f = ~ {
    sim_data %>%
      summarize(across(.cols = c(w, y), .fns = sum)) %>%
      mutate(s = y - (.x * w)) %>%
      pull(s)
  }),
  Incorrect = map_dbl(.x = parameter_grid, .f = ~ {
    sim_data %>%
      mutate(s = (w ^ (1 - 1.3)) * (y - w * .x)) %>%
      pull(s) %>%
      sum()
  })
) %>%
  glimpse()

plt <- estimating_equation_curves %>%
  pivot_longer(-parameter) %>%
  ggplot(aes(x = parameter, y = value, linetype = name)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotdash") +
  # geom_vline(xintercept = 114.5185, linetype = 2) +
  # geom_vline(xintercept = 121.7328, linetype = 1) +
  # geom_vline(xintercept = 200, linetype = 11) +
  theme_bw() +
  ylab("Score function") +
  xlab("Parameter") +
  theme(legend.position = "bottom", legend.title = element_blank())

plt

ggsave(plot = plt, filename = "tweedie-image.png",
       device = "png", units = "px", width = 900, height = 800,
       path = fs::path("posts", "tweedie-with-identity-link-and-offset"))
