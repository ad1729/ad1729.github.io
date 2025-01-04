
# The claims datasets used below are downloaded from openml.
# claim frequency data: https://www.openml.org/d/41214
# claim severity data: https://www.openml.org/d/41215
#
# These are available in the CASdatasets R package
# https://dutangc.github.io/CASdatasets/reference/freMTPL.html
# The counts for these data mentioned in the R package doc page differ slightly
# from the openml ones.

library(tidyverse)

post_dir <- fs::path("posts", "fitting-tweedie-models-to-claims-data") %>% print()
get_path <- function(fname) fs::path(post_dir, fname)

freq <- read_csv(get_path("claim_freq_raw.csv")) %>% glimpse()

severity <- read_csv(get_path("claim_severity_raw.csv")) %>% glimpse()

claims_subset <- read_csv(get_path("claims_subset.csv")) %>% glimpse()

# single row per policy?
freq %>% distinct(IDpol) %>% nrow()
freq %>% nrow()

# single claim amount per policy? - no, as expected
severity %>% distinct(IDpol) %>% nrow()
severity %>% nrow()

# aggregate the claim amounts data to get total claim amounts per policy
severity <- severity %>%
  group_by(IDpol) %>%
  summarise(ClaimAmount = sum(ClaimAmount), .groups = 'drop') %>%
  glimpse()

claims_data <- freq %>%
  left_join(x = ., y = severity, by = "IDpol") %>%
  mutate(
    across(.cols = c("ClaimAmount"), .fns = ~ replace_na(.x, 0))
  ) %>%
  glimpse()

# frequency distribution of claims
claims_data %>% count(ClaimNb)

# final check after the join
claims_data %>% distinct(IDpol) %>% nrow()
claims_data %>% nrow()

# compare the entries for a handful of overlapping individuals present in both
# the subset and the full dataset
set.seed(4)
sample_ids <- claims_subset %>%
  select(IDpol, ClaimAmount) %>%
  mutate(is_zero = ClaimAmount == 0) %>%
  group_by(is_zero) %>%
  slice_sample(n = 7) %>%
  ungroup() %>%
  print(n = Inf) %>%
  pull(IDpol) %>%
  print()

# claim amounts seem to match from the two datasets
bind_rows(
  claims_subset %>%
    filter(IDpol %in% sample_ids) %>%
    mutate(type = "subset", .before = 1),
  claims_data %>%
    filter(IDpol %in% sample_ids) %>%
    mutate(type = "full", .before = 1)
) %>%
  arrange(IDpol) %>%
  View(title = "sanity check")

claims_data %>%
  #select(IDpol, Exposure, ClaimAmount) %>%
  write_rds(x = ., file = get_path("claims.rds"))
