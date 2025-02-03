

library(tidyverse)
library(here)
library(birdnames)
library(MASS)
library(AICcmodavg)

#source("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/ACR_waterbird_data_management/code/utils.r")
#source("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/waterbird_data_work/code/utility/waterbird_utility_functions.R")
#source("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/waterbird_analyses/Tomales_waterbird_trends_2020/code/analysis_utilities.R")


source(here("code/utilities.R"))
source(here("code/prepare_predictors.R"))


custom_bird_list <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/my_R_general/birdnames_support/data/custom_bird_list")

wbird_keep_taxa <- c("AMCOGRSCLESCBUFF", "AMCO", "COGA", "Anseriformes", "Alcidae", "Gaviidae", "Pelecanidae", "Podicipediformes", "Sterninae", "Suliformes")

wbird_analysis_spp = c("ALL_WBIRD", "BRAN", "CANG", "GADW", "AMWI", "MALL", "NOPI", "GWTE", "SCAUP", "SUSC", "BLSC", "BUFF", "COGO", "COME", "RBME", "RUDU", "PBGR", "HOGR", "RNGR", "EAGR", "WCGR", "AMCO", "FOTE", "RTLO", "PALO", "COLO", "BRAC", "PECO", "DCCO", "BRPE")

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
wbird <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/water_birds/ACR_waterbird_data_management/data_files/working_rds/new_neg_machine_bay_total") %>% 
  add_study_year() %>%
  filter(!date %in% exclude_dates, study.year > 1991) %>% 
  bird_taxa_filter(wbird_keep_taxa) %>% 
  mutate(alpha.code = ifelse(alpha.code %in% c("GRSC", "LESC"), "SCAUP", alpha.code),
         alpha.code = ifelse(alpha.code %in% c("WEGR", "CLGR"), "WCGR", alpha.code)) 



wbird_out <- wbird %>% 
  mutate(alpha.code = "ALL_WBIRD") %>% 
  bind_rows(wbird) %>% 
  group_by(study.year, date, alpha.code) %>% 
  summarise(bay.total = sum(bay.total)) %>% 
  ungroup() %>% 
  filter(alpha.code %in% wbird_analysis_spp) %>% 
  group_by(study.year, alpha.code) %>% 
  summarise(p75.abund = floor(quantile(bay.total, 0.75))) %>% 
  ungroup()


# fill in 0 for years when spp not detected
wbird_out_full <- wbird_out %>%
  pivot_wider(id_cols = alpha.code, names_from = study.year, values_from = p75.abund) %>% 
  pivot_longer(-alpha.code, names_to = "study.year", values_to = "p75.abund") %>% 
  mutate(p75.abund = ifelse(is.na(p75.abund), 0, p75.abund),
         study.year = as.numeric(study.year)) %>% 
  left_join(mean_north_moci) %>% 
  left_join(lacr) %>% 
  filter(study.year >= 1995)

# fit models ----

zmaxit = 400

fit_wbird_mod_set <- function(zspp) {
  zspp_annual <- wbird_out_full %>% 
    filter(alpha.code == zspp)
  # trend and both environmental
  mods <- list(
  year2_flow_moci = glm.nb(p75.abund ~ poly(study.year, 2) + flow + moci, data = zspp_annual, maxit = zmaxit),
  year2_flow = glm.nb(p75.abund ~ poly(study.year, 2) + flow, data = zspp_annual, maxit = zmaxit),
  year2_moci = glm.nb(p75.abund ~ poly(study.year, 2) + moci, data = zspp_annual, maxit = zmaxit),
  
  year_flow_moci = glm.nb(p75.abund ~ study.year + flow + moci, data = zspp_annual, maxit = zmaxit),
  year_flow = glm.nb(p75.abund ~ study.year + flow, data = zspp_annual, maxit = zmaxit),
  year_moci = glm.nb(p75.abund ~ study.year + moci, data = zspp_annual, maxit = zmaxit),
  
  year = glm.nb(p75.abund ~ study.year, data = zspp_annual, maxit = zmaxit),
  flow = glm.nb(p75.abund ~ flow, data = zspp_annual, maxit = zmaxit),
  moci = glm.nb(p75.abund ~ moci, data = zspp_annual, maxit = zmaxit),
  
  intercept = glm.nb(p75.abund ~ 1, data = zspp_annual, maxit = zmaxit)
  )
  return(mods)
  
}


wbird_mods <- map(wbird_analysis_spp, fit_wbird_mod_set)

names(wbird_mods) <- wbird_analysis_spp


# model averaged estimates ----

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
wbird_mod_avg_predicter <- function(zspp) {
  
  spp_mods <- wbird_mods[[zspp]]
  
  znewdat = wbird_out_full %>% 
    filter(alpha.code == zspp) %>% 
    dplyr::select(study.year, alpha.code, p75.abund) %>% 
    mutate(flow = 0,
           moci = 0)
    
     
  
  all_best_pred = modavgPred(spp_mods, names(spp_mods), znewdat, se.fit=TRUE, type='response') %>% 
    data.frame() %>% 
    dplyr::select(-type, -contains("matrix")) %>% 
    cbind(znewdat)
  
  return(all_best_pred)
  
}


wbird_preds <- map_df(wbird_analysis_spp, wbird_mod_avg_predicter)
  
  
saveRDS(wbird_preds, here("model_objects/wbird_preds"))
