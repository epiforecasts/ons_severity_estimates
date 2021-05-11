library("here")
library("readxl")
library("rvest")
library("purrr")
library("tidyr")
library("dplyr")
library("janitor")
library("socialmixr")
library("lubridate")
library("ggplot2")
library("data.table")

# Read in infection estimates
infections <- data.table::as.data.table(readRDS("inf.rds"))


# Read in latest ONS estimates
url <- paste0("https://www.ons.gov.uk/peoplepopulationandcommunity/",
              "healthandsocialcare/conditionsanddiseases/datasets/",
              "coronaviruscovid19infectionsurveydata")
session <- session(url)

file_url <- session %>%
 html_nodes(xpath = paste0(
  "//a[contains(concat(' ', ",
  "normalize-space(@class),' '),' btn--thick ')]"
 )) %>%
 html_attr("href") %>%
 pluck(1)

file_name <- sub("^.*/([^/]+)$", "\\1", file_url)

cis_dir <- here::here()

if (!file.exists(file.path(cis_dir, file_name))) {
 download.file(paste0("https://www.ons.gov.uk", file_url),
               file.path(cis_dir, file_name))
}

age_region_data_raw <- read_excel(file.path(cis_dir, file_name),
                                  sheet = "1a",
                                  skip = 6)

ard <- age_region_data_raw %>%
 dplyr::select(1:4)

colnames(ard) <- c("period", "middle", "lower", "upper")

ard <- ard %>%
 dplyr::filter(!is.na(middle)) %>%
 dplyr::mutate(start_date = dmy(sub("^(.*) to .*$", "\\1", period)),
        end_date = start_date + 13) %>%
 dplyr::select(-period)

ons_eng <- as.data.table(ard)

ons_eng$geography <- "England"

# Variable for midpoint of 14 day intervals
ons_eng[, date := lubridate::ymd(start_date) + 6.5]

ons_pop <- fread(here("ons_data.csv"))
ons_eng$population_size <- sum(ons_pop$ons_population)

# A function to provide estimates of the detection probability for any continuous value of time since infection
# Have this data table in your environment
dt <- data.table::fread(here("fitted_params.csv"))

prob_cont <- function(time_since_inf, q = 0.5) {
  dt[, p := boot::inv.logit(beta1 + beta2 * (time_since_inf - cutpoint) + (time_since_inf - cutpoint) * beta3 * beta2 * fifelse(time_since_inf - cutpoint > 0, 1, 0))]
  out <- quantile(dt$p, probs = q)
  return(out)
}

# A nice speedy way of calculating number of detectable cases from
# a time series of daily infections
pvec <- vapply(X = 0:(365 * 2), FUN = prob_cont, FUN.VALUE = 1, q = 0.5)

detectable_cases <- function(x) {
  
  out <- rep(0, length(x))
  for(i in 1:length(x)) {
    for(j in i:1) {
      out[i] <- out[i] + x[i] * pvec[i - j + 1]
    }
  }
  
  return(out)
}

infection <- infections[, value_cum := detectable_cases(value), by = c("data_source", "geography", "sample")]

# Select relevant dates and geographies from infection samples and summarise over all samples
rel_inf <- infection[date %in% lubridate::ymd(ons_eng$date)][, .(inf_median = median(value_cum), inf_upper = quantile(value_cum, 0.975), 
                                                                 inf_lower = quantile(value_cum, 0.025)), by = c("data_source", "geography", "date")]


rel_inf[, date := lubridate::ymd(date)]
ons_eng[, date := lubridate::ymd(date)]
dt <- merge.data.table(rel_inf, ons_eng, by = c("geography", "date"))

# Calculate prevalence from infection samples
dt[, inf_median := 100 * inf_median / population_size
][,inf_lower := 100 * inf_lower / population_size
][,inf_upper := 100 * inf_upper / population_size]

dt[, data_source_f := factor(data_source, levels = c("cases", "admissions", "deaths"))]

# Plot IHR / ICR / IFR
p_out <- dt[date > "2020-10-01"] %>%
 ggplot(aes(x = date, 
            y =  100 * inf_median / middle,
            ymin = 100 * inf_lower / upper,
            ymax = 100 * inf_upper / lower)) +
 geom_point() +
 geom_errorbar() + 
 # geom_line() +
 facet_grid(data_source_f ~ geography, scales = "free") +
 labs(x = "ONS 14-day interval midpoint", 
      y = "Estimated IFR / IHR / ICR (%)") +
 cowplot::theme_minimal_grid()  +
 scale_x_date(breaks = "1 month", date_labels = "%b %Y") +
 theme(strip.text.x = element_blank())

ggsave(p_out, filename = "ons_severity.pdf")
 
