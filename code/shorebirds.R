

library(tidyverse)
library(here)
library(birdnames)
library(MASS)
library(AICcmodavg)

options(scipen = 999)


source(here("code/utilities.R"))
source(here("code/prepare_predictors.R"))

exclude_date <- as.Date(c("1990-01-04", "1990-02-12", "2010-01-18"))


sbird_analysis_spp <- c("ALL_SBIRD", "DUNL", "WESA", "MAGO", "LESA", "SAND", "WILL", "DOSP", "BBPL", "BLTU", "YELL", "SEPL", "KILL"
               #, "SPSA"
)


# data prep ----

sbirds <- readRDS(here("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/shorebirds/ACR_shorebird_data_management/data_files/rds/shorebirds_for_analysis")) %>% 
  filter(str_detect(season.year, "winter")) %>% 
  mutate(study.year = str_replace(season.year, "winter_", ""),
         study.year = as.numeric(study.year)) %>% 
  dplyr::select(-season.year)
  

# annual estimate each species and all species combined
sbirds_out <- sbirds %>% 
  mutate(alpha.code = "ALL_SBIRD") %>% 
  bind_rows(sbirds) %>% 
  group_by(study.year, date, alpha.code) %>% 
  summarize(bay.total = sum(count)) %>% 
  ungroup() %>% 
  group_by(study.year, alpha.code) %>% 
  summarise(p75.abund = quantile(bay.total, probs = c(0.75))) %>% 
  ungroup()


# fill 0s
sbirds_out_full <- sbirds_out %>% 
  pivot_wider(id_cols = alpha.code, names_from = study.year, values_from = p75.abund) %>% 
  pivot_longer(-alpha.code, names_to = "study.year", values_to = "p75.abund") %>% 
  mutate(p75.abund = ifelse(is.na(p75.abund), 0, p75.abund),
         study.year = as.numeric(study.year),
         giac.dummy = ifelse(study.year < 2009, 0, 1)) %>% 
  left_join(., rain)

zmaxit = 400


fit_sbird_mod_set <- function(zspp) {
  sbirds_winter <- sbirds_out_full %>% 
    filter(alpha.code == zspp)
  
  # big_mod <- glm.nb(floor(p75.abund) ~ poly(year, 2) * section + giac.dummy + seas.rain.mm + mean.bird.hunters, data = sbirds_winter, maxit = zmaxit)
  mods <- list(
  year2_rain_giac = glm.nb(floor(p75.abund) ~ poly(study.year, 2) + seas.rain.mm + giac.dummy, data = sbirds_winter, maxit = zmaxit),
  year2_rain = glm.nb(floor(p75.abund) ~ poly(study.year, 2) + seas.rain.mm, data = sbirds_winter, maxit = zmaxit),
  year2_giac = glm.nb(floor(p75.abund) ~ poly(study.year, 2) + giac.dummy, data = sbirds_winter, maxit = zmaxit),
  
  year_rain_giac = glm.nb(floor(p75.abund) ~ study.year + seas.rain.mm + giac.dummy, data = sbirds_winter, maxit = zmaxit),
  year_rain = glm.nb(floor(p75.abund) ~ study.year + seas.rain.mm, data = sbirds_winter, maxit = zmaxit),
  year_giac = glm.nb(floor(p75.abund) ~ study.year + giac.dummy, data = sbirds_winter, maxit = zmaxit),
  
  year2 = glm.nb(floor(p75.abund) ~ poly(study.year, 2), data = sbirds_winter, maxit = zmaxit),
  year = glm.nb(floor(p75.abund) ~ study.year, data = sbirds_winter, maxit = zmaxit),
  rain = glm.nb(floor(p75.abund) ~ seas.rain.mm, data = sbirds_winter, maxit = zmaxit),
  giac = glm.nb(floor(p75.abund) ~ giac.dummy, data = sbirds_winter, maxit = zmaxit),
  
  intercept = glm.nb(floor(p75.abund) ~ 1, data = sbirds_winter, maxit = zmaxit)
  )
}



sbird_mods <- map(sbird_analysis_spp, fit_sbird_mod_set)

names(sbird_mods) <- sbird_analysis_spp




# ---
sbird_mod_avg_predicter <- function(zspp){
  
  spp_mods <- sbird_mods[[zspp]]
  
  
  newdat <- sbirds_out_full %>% 
    filter(alpha.code == zspp) %>% 
    dplyr::select(study.year, alpha.code, p75.abund) %>% 
    mutate(giac.dummy = ifelse(study.year < 2009, 0, 1),
           seas.rain.mm = 0)  
  
  all_best_pred = modavgPred(spp_mods, names(spp_mods), newdat, se.fit=TRUE, type='response') %>% 
    data.frame() %>% 
    dplyr::select(-type, -contains("matrix")) %>% 
    cbind(newdat)
  
  return(all_best_pred)
}



sbird_preds <- map_df(sbird_analysis_spp, sbird_mod_avg_predicter)



saveRDS(sbird_preds, here("model_objects/sbird_preds"))


