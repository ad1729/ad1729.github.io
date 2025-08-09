
library(dplyr)
library(purrr)
library(rlang)
library(mirai)

foo <- function(x) {
  x + 2
}

bar <- function(y) {
  y %>%
    foo() %>%
    exp()
}

baz <- function(z) {
  z %>%
    bar() %>%
    log() %>%
    {. - 2}
}

baz(7)

baz_slow <- function(z) {
  Sys.sleep(1)
  baz(z)
}

baz_slow(7)

iris_tbl <- iris %>% as_tibble() %>% mutate(row_id = row_number())
iris_tbl

slice_df <- function(z) {
  row_id <- baz_slow(z)
  start <- row_id + 10
  end <- start + 3
  iris_tbl %>%
    dplyr::slice(start:end) %>%
    list()
}

slice_df(1)

test_env <- rlang::env(
  foo = foo,
  bar = bar,
  baz = baz,
  baz_slow = baz_slow,
  slice_df = slice_df
)

env_print(test_env)

daemons(5)

crate_obj <- carrier::crate(
  #function(z) baz_slow(z),
  function(z) slice_df(z),
  # baz_slow = baz_slow,
  # baz = baz,
  #.parent_env = global_env()
  .parent_env = test_env
)

crate_obj

crate_obj(1:5)

env_names(global_env())

system.time({res_seq <- map(.x = 1:5, .f = \(z) slice_df(z))})
system.time(res_par <- crate_obj(1:5))
system.time({res_par2 <- map(.x = 1:5, .f = ~ crate_obj(.x))})

res_seq
res_par
identical(res_seq, res_par)

class(res_par)

map(.x = 1:5, .f = crate_obj)

system.time({map(.x = 1:5, .f = crate_obj)})

daemons(0)

# reprex for crate package issue
reprex::reprex({
  library(dplyr)
  library(purrr)
  library(carrier)
  library(mirai)
  library(rlang)

  c("dplyr", "purrr", "carrier", "mirai", "rlang") %>%
    set_names(nm = ~ .x) %>%
    map(.f = ~ packageVersion(.x))

  iris_tbl <- iris %>%
    as_tibble() %>%
    mutate(row_id = row_number())

  helper_function <- function(x) {
   x %>% exp() %>% log()
  }

  slice_df_single_row <- function(z) {
    Sys.sleep(1)
    z <- helper_function(z)
    iris_tbl %>%
      dplyr::slice(z)
  }

  slice_df_multiple_rows <- function(z) {
    Sys.sleep(1)
    z <- helper_function(z)
    start <- z + 10
    end <- start + 3
    iris_tbl %>%
      dplyr::slice(start:end)
  }

  slice_df_single_row(5)
  slice_df_multiple_rows(5)

  system.time({test_seq_single <- map_dfr(.x = 1:5, .f = ~ slice_df_single_row(.x))})
  test_seq_single

  system.time({test_seq_multiple <- map_dfr(.x = 1:5, .f = ~ slice_df_multiple_rows(.x))})
  test_seq_multiple

  daemons(5)

  test_env <- env(
    `%>%` = magrittr::`%>%`,
    helper_function = helper_function,
    slice_df_single_row = slice_df_single_row,
    slice_df_multiple_rows = slice_df_multiple_rows
  )

  crate_obj_single <- crate(function(z) slice_df_single_row(z), .parent_env = test_env)
  crate_obj_single

  crate_obj_multiple <- crate(function(z) slice_df_multiple_rows(z), .parent_env = test_env)
  crate_obj_multiple

  system.time({test_par_single <- crate_obj_single(1:5)})
  test_par_single

  system.time({test_par_multiple <- crate_obj_multiple(1:5)})
  test_par_multiple

  daemons(0)
})
