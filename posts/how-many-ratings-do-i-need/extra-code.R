
# make sure all the code in the index.qmd file is run before running the code below

simulate_paths <- function(n_ratings, probs) {
  rmultinom(1, n_ratings, probs) %>%
    t() %>%
    as_tibble(.name_repair = ~ paste("v", 1:length(probs), sep = ""))
}

map_dfr(.x = 1:3, .f = ~ simulate_paths(1000, c(1, 1, 1)))


simulated_paths <- possible_ratings %>%
  select(x1:x5) %>%
  # create a group id for each of the possible ratings vectors
  mutate(group = row_number(), .before = x1) %>%
  # number of additional ratings
  expand_grid(horizon = seq(100, 1000, 100),
              sim_id = 1:10000) %>%
  rowwise() %>%
  mutate(simulate_paths(horizon, c(x1, x2, x3, x4, x5))) %>%
  ungroup()

glimpse(simulated_paths)


mean_ratings <- simulated_paths %>%
  mutate(total1 = x1 + v1,
         total2 = x2 + v2,
         total3 = x3 + v3,
         total4 = x4 + v4,
         total5 = x5 + v5) %>%
  rowwise() %>%
  mutate(mean_rating = weighted.mean(1:5, w = c(total1, total2, total3, total4, total5))) %>%
  ungroup()

mean_ratings


mean_ratings %>%
  slice_sample(n = 30000) %>%  # sample subset of simulations to speed up plotting
  mutate(color = ifelse(mean_rating >= 4.6, "orange", "black"),
         alpha = ifelse(mean_rating >= 4.6, 0.7, 0.2)) %>%
  ggplot(aes(x = horizon, y = mean_rating, color = color, alpha = alpha)) +
  geom_point(position = position_jitter(width = 10, seed = 1234)) +
  geom_hline(yintercept = 4.6) +
  theme_classic() +
  facet_wrap(~ group, ncol= 3) +
  scale_color_identity() +
  scale_alpha_identity()


mean_ratings %>%
  group_by(group, horizon) %>%
  summarize(target = mean(mean_rating >= 4.6)) %>%
  mutate(group = factor(group)) %>%
  ggplot(aes(x = horizon, y = target, color = group, group = group)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  scale_y_continuous(labels = scales::label_percent())
