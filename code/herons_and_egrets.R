
library(tidyverse)
library(here)
library(birdnames)

options(scipen = 999)
# for functions to read in Access data and do standard cleaning
source("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/HEP/HEP_data_work/HEP_code/HEP_utility_functions.R")

source(here("code/utilities.R"))

hepdata_location = here("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/HEP/HEP_data_work/HEP_data/HEPDATA.accdb")

# tomales bay colonies
tbay_parent_codes <- c(50, 152, 119, 122, 143, 83, 114, 32, 160, 113)
# other west marin colonies
w_marin_parent_codes <- c(137, 94, 17, 1, 53, 186, 47, 71)


# prepare data ----
hep_sites <- hep_sites_from_access(hepdata_location) %>% 
  dplyr::select(code, parent.code, site.name, parent.site.name, subregion)


# calculate a cumulative rain index following Stenzel and Page 2018
rain_lag <- readRDS(here("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/HEP/HEP_data_work/HEP_data/subreg_rain")) %>% 
  dplyr::select(year = birdyear, subregion, subreg.name, subreg.rain) %>%
  data.frame() %>% 
  arrange(subregion, year) %>% 
  group_by(subregion) %>% 
  mutate(subreg.rain = subreg.rain + (lag(subreg.rain)/2) + (lag(subreg.rain, 2)/3),
         subreg.rain = scale(subreg.rain)[,1]) %>% 
  ungroup()

# total nests counted across the study area

hep_abund <- readRDS(here("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/HEP/HEP_data_work/HEP_data/hep_annual_nest_abundance")) %>% 
  mutate(common.name = translate_bird_names(species, "alpha.code", "common.name")) %>% 
  full_join(hep_sites)%>% 
  cut_never_nested() %>% 
  filter(peakactvnsts > 0, year >= 1990, year != 2020)

calculate_tot_nests <- function(df) {
  df <- df %>% 
    summarise(tot.nests = sum(peakactvnsts),
              num.colonies = n()) %>% 
    ungroup()
}


hep_abund_with_all <- hep_abund %>% 
  bind_rows(hep_abund %>% 
              filter(species %in% c("GBHE", "GREG")) %>% 
              mutate(species = "All",
                     common.name = "Large herons and egrets"))


# entire study area
hep_abund_study_area <- hep_abund_with_all %>% 
  group_by(year, species, common.name) %>% 
  calculate_tot_nests() %>% 
  mutate(subregion = "All")

# OUC subregion

hep_abund_ouc <- hep_abund_with_all %>% 
  filter(subregion == "OUC") %>% 
  group_by(year, subregion, species, common.name) %>% 
  calculate_tot_nests()


# OUC subregion without Tomales
hep_abund_ouc_no_tbay <- hep_abund_with_all %>% 
  filter(subregion == "OUC" & !parent.code %in% tbay_parent_codes) %>% 
  group_by(year, subregion, species, common.name) %>%
  calculate_tot_nests() %>% 
  mutate(area = "OUC no Tomales Bay")


# Tomales Bay only - need to keep subregion as OUC to join with that rain data
hep_abund_tbay <- hep_abund_with_all %>% 
  filter(parent.code %in% tbay_parent_codes) %>% 
  group_by(year, subregion, species, common.name) %>% 
  calculate_tot_nests() %>% 
  mutate(area = "Tomales Bay")



# combine into a single long format table with duplicates of data representing different scales of aggregation
trend_analysis_table <- bind_rows(hep_abund_study_area,
                                  hep_abund_ouc,
                                  hep_abund_tbay,
                                  hep_abund_ouc_no_tbay) %>%
  filter(species %in% c("GREG", "GBHE", "All", "DCCO")) %>% 
  dplyr::select(year, subregion, area, species, common.name, tot.nests) %>% 
  left_join(rain_lag) %>% 
  mutate(subreg.name = case_when(area == "Tomales Bay" ~ area,
                                 area == "OUC no Tomales Bay" ~ "Outer Coast, no Tomales Bay",
                                 TRUE ~ as.character(subreg.name)),
         subregion = case_when(area == "Tomales Bay" ~ "TBAY",
                               area == "OUC no Tomales Bay" ~ "OUCnoTBAY",
                               TRUE ~ as.character(subregion)))

#saveRDS(trend_analysis_table, here("data/hep_trend_analysis_table"))

spp_subreg <- distinct(trend_analysis_table, species, subregion) %>% 
  mutate(spp.subreg = paste(species, subregion, sep = "_"))


subreg_mean_rain_lag <- trend_analysis_table %>% 
  distinct(subregion, subreg.name, subreg.rain) %>% 
  group_by(subregion, subreg.name) %>% 
  summarise(mean.subreg.rain = mean(subreg.rain, na.rm = TRUE))

# model fitting ----
#' Fit negative binomial glm on untransformed nest abundance with year^2 as predictor
#'
#' @param zspp 
#' @param zsubreg 
#'
#' @return
#' @export
#' 
#' @details
#' modified from "C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/HEP/hep_analyses/how_are_the_egrets_doing/code/ms_analysis/hep_trend_utilities.R"
#' 
#'
#' @examples
#' spp_subreg_mods_glmnb <- map2(spp_subreg$species, spp_subreg$subregion, fit_mods_glmbn)
#' names(spp_subreg_mods_glmnb) <- spp_subreg$spp.subreg
fit_mods_glmbn_year2 <- function(zspp, zsubreg) {
  
  zdat <- trend_analysis_table %>% 
    filter(species == zspp, subregion == zsubreg)
  
  zmod <- MASS::glm.nb(data = zdat, formula = tot.nests ~ poly(year, 2) + subreg.rain)
  
  }


spp_subreg_mods <- map2(spp_subreg$species, spp_subreg$subregion, fit_mods_glmbn_year2)
names(spp_subreg_mods) <- spp_subreg$spp.subreg

# extract model output ----
# year estimates ----
# get_preds_glmnb copied from "C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/HEP/hep_analyses/how_are_the_egrets_doing/code/ms_analysis/hep_trend_utilities.R"


#' Extract and backtransform model predictions and their CI from the Negative Binomial GLMs
#'
#' @param zmod 
#' @param zmod.name 
#'
#' @return
#' @export
#'
#' @examples
#' preds_glmnb <- map2_df(spp_subreg_mods_glmnb, names(spp_subreg_mods_glmnb), get_preds_glmnb)
get_preds_glmnb <- function(zmod, zmod.name) { 
  
  spp_subreg <- zmod.name %>% 
    data.frame() %>% 
    rename(spp.subreg = 1) %>% 
    separate(spp.subreg, into = c("species", "subregion"), sep = "_")
  
  ana_table <- trend_analysis_table %>% 
    right_join(spp_subreg)
  
  sub_rain <- subreg_mean_rain_lag %>%  
    #readRDS(here("data/subreg_rain")) %>%
    right_join(spp_subreg) %>% 
    distinct(mean.subreg.rain)
  
  znewdat = data.frame(year = seq(min(ana_table$year), max(ana_table$year)),
                       subreg.rain = sub_rain$mean.subreg.rain)
  
  ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    cbind(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit))) %>% 
    bind_cols(spp_subreg)
}



preds <- map2_df(spp_subreg_mods, names(spp_subreg_mods), get_preds_glmnb) %>% 
  full_join(trend_analysis_table %>% 
              select(year, subregion, subreg.name, species, common.name, tot.nests)) %>% 
  mutate(#common.name = translate_bird_names(species, "alpha.code", "common.name"),
         subreg.name = factor(subreg.name, levels = c("Tomales Bay", "Outer Pacific Coast, North", "Outer Coast, no Tomales Bay", "Entire study area")))

saveRDS(preds, "data/HEP_predictions")

# model coeficients and 95% CI ----


coef_ci <- map2_df(spp_subreg_mods, names(spp_subreg_mods), get_coefs_cis)
saveRDS(coef_ci, "data/HEP_coef_ci")

# estimated abundance across the range of rain values ----


#' Extract and backtransform model predictions and their CI from the Negative Binomial GLMs
#'
#' @param zmod 
#' @param zmod.name 
#'
#' @return
#' @export
#'
#' @examples
#' preds_glmnb <- map2_df(spp_subreg_mods_glmnb, names(spp_subreg_mods_glmnb), get_preds_glmnb)
get_rain_preds_glmnb <- function(zmod, zmod.name) { 
  
  spp_subreg <- zmod.name %>% 
    data.frame() %>% 
    rename(spp.subreg = 1) %>% 
    separate(spp.subreg, into = c("species", "subregion"), sep = "_")
  
  ana_table <- trend_analysis_table %>% 
    right_join(spp_subreg)
  
  #zsubreg.rain = seq(min(ana_table$subreg.rain, na.rm = TRUE), max(ana_table$subreg.rain, na.rm = TRUE), length.out = 10)
  
  #zsubreg.rain = c(mean(ana_table$subreg.rain, na.rm = TRUE) - (0.5 * sd(ana_table$subreg.rain, na.rm = TRUE)), mean(ana_table$subreg.rain, na.rm = TRUE) + (0.5 * sd(ana_table$subreg.rain, na.rm = TRUE)))
  
  #zsubreg.rain = c(quantile(ana_table$subreg.rain, 0.25, na.rm = TRUE),
  #                                     quantile(ana_table$subreg.rain, 0.75, na.rm = TRUE))
  
  zsubreg.rain = c(-1, 0, 1)
  
  znewdat = data.frame(year = floor(mean(trend_analysis_table$year)),
                       subreg.rain = zsubreg.rain)
  
  
  ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    cbind(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit))) %>% 
    bind_cols(spp_subreg)
}



hep_rain_preds <- map2_df(spp_subreg_mods, names(spp_subreg_mods), get_rain_preds_glmnb)
saveRDS(hep_rain_preds, here("data/hep_rain_preds"))


