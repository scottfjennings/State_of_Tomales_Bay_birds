

library(tidyverse)
library(here)
library(birdnames)
library(MASS)

source("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/ACR_waterbird_data_management/code/utils.r")
source("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/waterbird_data_work/code/utility/waterbird_utility_functions.R")
source("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/waterbird_analyses/Tomales_waterbird_trends_2020/code/analysis_utilities.R")


source(here("code/utilities.R"))

custom_bird_list <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/my_R_general/birdnames_support/data/custom_bird_list")


wbird_keep_taxa <- c("AMCOGRSCLESCBUFF", "AMCO", "COGA", "Anseriformes", "Alcidae", "Gaviidae", "Pelecanidae", "Podicipediformes", "Sterninae", "Suliformes")

options(scipen = 999)
# data prep ----
# several dates had data quirks either in the field or on data sheets and should be excluded from analysis
exclude_dates <- as.Date(c("2017-12-16", # Beaufort 3-4
                           "2011-12-17", # Bivalve count incomplete; ponds only
                           "2006-02-11", # WCD not surveyed due to fog
                           "2004-12-18", # inverness, millerton not surveyed
                           "1999-12-18", # survey done north to south
                           "1999-01-30", # strong wind, not entered in raw tallies but including to have a complete coded reference to excluded dates
                           #"1996-12-21", # no WCD precount, never many birds there so maybe ok to include. Also multiple data sheet ambiguities, entered raw tally data likely ok but difficult to be sure
                           "1997-12-20", # apparent large swell and strong wind in north bay
                           #"1993-02-06", # no CGP precount, never many birds there so maybe ok to include
                           #"1993-01-23", # no WCD precount, never many birds there so maybe ok to include.
                           "1992-01-11", # survey apparently done north to south
                           "1991-12-14", # substantial ambiguity in precount data sheets. only millerton, wcd entered
                           "1990-12-15", # unclear labeling of section/transect on data sheets and can't reverse engineer. Data not fully entered; date should be excluded from analysis
                           "1989-12-16" # funky protocol. No precounts. Should probably exclude from analysis 
))

# add study year, combine groupd spp for analysis, and recalculate the total number of birds for each species or species group each day
spp_day_total <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/ACR_waterbird_data_management/data_files/working_rds/new_neg_machine_bay_total") %>% 
  wbird_add_study_day() %>% # from waterbird_utility_functions.R
  filter(!date %in% exclude_dates, study.year > 1991) %>% 
  bird_taxa_filter(wbird_keep_taxa) %>% 
  mutate(alpha.code = ifelse(alpha.code %in% c("GRSC", "LESC"), "SCAUP", alpha.code),
         alpha.code = ifelse(alpha.code %in% c("WEGR", "CLGR"), "WCGR", alpha.code)) %>% 
  group_by(study.year, date, alpha.code) %>% 
  summarise(bay.total = sum(bay.total)) %>% 
  ungroup()


# same as above but without SCAUP and WCGR lumped
spp_day_total_ungrouped <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/ACR_waterbird_data_management/data_files/working_rds/new_neg_machine_bay_total") %>% 
  wbird_add_study_day() %>% # from waterbird_utility_functions.R
  filter(!date %in% exclude_dates, study.year > 1991) %>% 
  bird_taxa_filter(wbird_keep_taxa) %>% 
  mutate(alpha.code = ifelse(alpha.code == "ACGO", "CACG", alpha.code)) %>% 
  group_by(study.year, date, alpha.code) %>% 
  summarise(bay.total = sum(bay.total)) %>% 
  ungroup()

# add up all species each year (includes non trend species), combine with by-species data, and calculate the 75th percentile of the individual day totals for each species/species group each year
spp_annual <- spp_day_total %>%
  group_by(study.year, date) %>% 
  summarise(bay.total = sum(bay.total)) %>% 
  mutate(alpha.code = "All") %>% 
  bind_rows(spp_day_total %>% dplyr::select(study.year, date, alpha.code, bay.total)) %>% 
  group_by(study.year, alpha.code) %>% 
  summarise(p75.abund = floor(quantile(bay.total, 0.75)))


# list of species detected >= 20 years. created by C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/waterbird_analyses/Tomales_waterbird_trends_2020/code/analysis1_prepare_data.R
wbird_trend_spp <- readRDS(here("data/wbird_trend_spp")) %>% 
  mutate(alpha.code = ifelse(alpha.code == "ALL", "All", alpha.code))


# fill in 0 for years when spp not detected
spp_annual_full <- spp_annual %>%
  filter(alpha.code %in% wbird_trend_spp$alpha.code) %>% 
  pivot_wider(id_cols = alpha.code, names_from = study.year, values_from = p75.abund) %>% 
  pivot_longer(-alpha.code, names_to = "study.year", values_to = "p75.abund") %>% 
  mutate(p75.abund = ifelse(is.na(p75.abund), 0, p75.abund),
         study.year = as.numeric(study.year))


# join bird data with predictors
spp_annual_full_preds <- spp_annual_full %>% 
  full_join(., readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/waterbird_analyses/Tomales_waterbird_trends_2020/data_files/predictors") %>% 
              rename(moci = mean.moci,
                     fresh = annual.freshwater) %>% 
              mutate(moci = scale(moci, scale = TRUE, center = TRUE)[,1],
                     fresh = scale(fresh, scale = TRUE, center = TRUE)[,1])) %>% 
  data.frame()


#saveRDS(spp_annual_full_preds, here("data_files/spp_annual_full_preds"))


fit_wbird_mod <- function(zspp) {
  zspp_annual <- spp_annual_full_preds %>% 
    filter(alpha.code == zspp)
  # trend and both environmental
  year2_fresh_moci <- glm.nb(p75.abund ~ poly(study.year, 2) + fresh + moci, data = zspp_annual)
  
  return(year2_fresh_moci)
  
}

all_spp_mod<- map(wbird_trend_spp$alpha.code, fit_wbird_mod)

names(all_spp_mod) <- wbird_trend_spp$alpha.code

# saveRDS(all_spp_mods, here("fitted_models/all_spp_mods"))


# model coeficients and 95% CI ----


coef_ci <- map2_df(all_spp_mod, names(all_spp_mod), get_coefs_cis)
saveRDS(coef_ci, "data/wbird_coef_ci")





# wbird_mod_preds adapted from mod_predictions_link in C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/waterbird_analyses/Tomales_waterbird_trends_2020/code/analysis_utilities.R


#' get_wbird_mod_preds
#' 
#' calculate model estimates (predictions) for zspp species and zmod model on the scale of the link function
#'
#' @param zspp 
#' @param zmod 
#'
#' @return data frame
#' @export
#'
#' @examples
get_wbird_mod_preds <- function(zmod, zmod.name) {
  
  ana_table <- spp_annual_full_preds %>% 
    filter(alpha.code == zmod.name)
  
  znewdat = data.frame(study.year = seq(min(ana_table$study.year), max(ana_table$study.year)),
                       fresh = mean(ana_table$fresh),
                       moci = mean(ana_table$moci))
  
  ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    bind_cols(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit))) %>% 
    full_join(ana_table %>% dplyr::select(study.year, alpha.code, p75.abund)) %>% 
    mutate(alpha.code = ifelse(is.na(alpha.code), zmod.name, alpha.code))
  
}


wbird_preds <- map2_df(all_spp_mod, names(all_spp_mod), get_wbird_mod_preds)


saveRDS(wbird_preds, here("data/wbird_preds"))

#' get_wbird_fresh_preds
#' 
#' calculate model estimates (predictions) for zspp species and zmod model on the scale of the link function
#'
#' @param zspp 
#' @param zmod 
#'
#' @return data frame
#' @export
#'
#' @examples
get_wbird_fresh_preds <- function(zmod, zmod.name) {
  
  ana_table <- spp_annual_full_preds %>% 
    filter(alpha.code == zmod.name)
  
#  znewdat = data.frame(study.year = mean(ana_table$study.year),
#                       fresh = seq(min(ana_table$fresh), max(ana_table$fresh), length.out = 10),
#                       moci = mean(ana_table$moci))
  
  
  #znewdat = data.frame(study.year = floor(mean(ana_table$study.year)),
  #                     fresh = c(mean(ana_table$fresh) - (0.5 * sd(ana_table$fresh)),
  #                               mean(ana_table$fresh) + (0.5 * sd(ana_table$fresh))),
  #                     moci = mean(ana_table$moci))
  
  #znewdat = data.frame(study.year = floor(mean(ana_table$study.year)),
  #                     fresh = c(quantile(ana_table$fresh, 0.25),
  #                               quantile(ana_table$fresh, 0.75)),
  #                     moci = mean(ana_table$moci))
  
  znewdat = data.frame(study.year = floor(mean(ana_table$study.year)),
                       fresh = c(-1, 0, 1),
                       moci = 0)
  
  
  ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    bind_cols(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit)),
           alpha.code = zmod.name)
  
}


wbird_fresh_preds <- map2_df(all_spp_mod, names(all_spp_mod), get_wbird_fresh_preds)

saveRDS(wbird_fresh_preds, here("data/wbird_fresh_preds"))


#' get_wbird_moci_preds
#' 
#' calculate model estimates (predictions) for zspp species and zmod model on the scale of the link function
#'
#' @param zspp 
#' @param zmod 
#'
#' @return data frame
#' @export
#'
#' @examples
get_wbird_moci_preds <- function(zmod, zmod.name) {
  
  ana_table <- spp_annual_full_preds %>% 
    filter(alpha.code == zmod.name)
  
# zmoci = seq(min(ana_table$moci), max(ana_table$moci), length.out = 10)
  
  #zmoci = seq(quantile(ana_table$moci, 0.25), quantile(ana_table$moci, 0.75))
  
  zmoci = c(-1, 0, 1)
  
  znewdat = data.frame(study.year = floor(mean(ana_table$study.year)),
                       fresh = 0,
                       moci = zmoci)
  

  
    ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    bind_cols(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit)),
           alpha.code = zmod.name) 
  
}


wbird_moci_preds <- map2_df(all_spp_mod, names(all_spp_mod), get_wbird_moci_preds)

saveRDS(wbird_moci_preds, here("data/wbird_moci_preds"))

