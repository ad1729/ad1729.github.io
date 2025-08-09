
library(tidyverse)
# library(flexsurv)

theme_set(theme_classic())

# trying this out for the simplest distribution - the exponential
summary(-log(runif(1e7)) / 5)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.05752 0.13868 0.20007 0.27735 3.16216
summary(rexp(1e7, rate = 5))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00000 0.05753 0.13863 0.19991 0.27709 3.08009

qexp(0.5, rate = 5)
log(1 - 0.5) + flexsurv::Hexp(0.1386294, rate = 5)
# CDFs for checking that the same parametrization is used
flexsurv::Hexp
pexp(0.1386294, rate = 5)
1 - exp(-flexsurv::Hexp(0.1386294, rate = 5))

sim_exponential <- function(u, lambda = 5) {
  -log(u) / lambda
}

obj_fn <- function(t, u) {
  log(1 - u) + flexsurv::Hexp(t, rate = 5)
}

grad_fn <- function(t, u) {
  flexsurv::hexp(t, rate = 5)
}

# so the median event time is 0.1386294 for Y ~ exponential(5)
# for solving this 1d optimization problem we can plot the objective function
# and the derivative for u = 0.5 while varying t
tibble(
  t = seq(0.001, 1.0, 0.001),
  Objective = obj_fn(t = t, u = 0.5),
  Gradient = grad_fn(t = t, u = 0.5)
) %>%
  pivot_longer(cols = c(Objective, Gradient), names_to = "var", values_to = "val") %>%
  ggplot(aes(x = t, y = val, group = var)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  #scale_x_continuous(transform = "log") +
  facet_wrap(~var)

uniroot_solution <- function(u) {
  sim_time <- uniroot(obj_fn, interval = c(0, 10), u = u)
  tibble(iter = sim_time$iter, obj_val = sim_time$f.root, time = sim_time$root)
}

qexp(0.5, rate = 5)
uniroot(obj_fn, interval = c(0, 10), u = 0.5)
optim(par = c(1), fn = obj_fn, gr = grad_fn, method = "L-BFGS-B", lower = 0, u = 0.5)
