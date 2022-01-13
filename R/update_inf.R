library(dplyr)
library(tidyr)
library(purrr)

source("https://gist.githubusercontent.com/sbfnk/d2900c745312219e3e48e08adde47cde/raw/c98fbdd738eafa71af12d05af0d3e068cf5b607b/get_covid19_nowcasts.r")

inf <- list()
inf[["cases"]] <- get_covid19_nowcasts("subnational/united-kingdom/cases", variable = "cases_by_infection")
inf[["admissions"]] <- get_covid19_nowcasts("subnational/united-kingdom/admissions", variable = "cases_by_infection")
inf[["deaths"]] <- get_covid19_nowcasts("subnational/united-kingdom/deaths", variable = "cases_by_infection")

gamma_approx <- function(x, n = 100) {
  return(tibble(sample = seq_len(n), 
                value = rgamma(n, shape = x$mean^2 / x$sd^2,
                              rate = x$mean / x$sd^2)))
}

inf_aggregate <- bind_rows(inf, .id = "data_source") %>%
  rename(geography = region) %>%
  filter(!is.na(sd), sd > 0) %>%
  nest(data = c(-data_source, -geography, -date)) %>%
  mutate(samples = map(data, gamma_approx)) %>%
  select(-data) %>%
  unnest(samples) %>%
  mutate(value = pmax(0, value))

dir.create(here::here("data-raw"), showWarnings = FALSE, recursive = TRUE)
saveRDS(inf_aggregate, here::here("data-raw", "inf.rds"))

