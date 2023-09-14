

library(tidyverse)
library(here)
library(birdnames)
library(MASS)

options(scipen = 999)


exclude_date <- as.Date(c("1990-01-04", "1990-02-12", "2010-01-18"))

sbirds <- readRDS(here("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/shorebirds/ACR_shorebird_data_management/data_files/rds/shorebirds_for_analysis"))%>% 
  rename(section = North_South_Code)


# rain from C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/general_data_sources/Rainfall/get_PRISM.R
rain <- read.csv("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/general_data_sources/Rainfall/derived_data/tomales_mean_month_rain.csv") %>% 
  mutate(year = ifelse(month < 7, year - 1, year)) %>% 
  filter(month <= 2 | month > 7) %>% 
  group_by(year) %>% 
  summarise(seas.rain.mm = sum(mean.rain),
            seas.rain.mm = round(seas.rain.mm)) %>% 
  ungroup() %>% 
  filter(year > 1988) %>% 
  mutate(seas.rain.mm = scale(seas.rain.mm)[,1])


# annual estimate each species
analysis_table <- sbirds %>% 
  group_by(season.year, date, alpha.code) %>% 
  summarize(total.sbirds = sum(count)) %>% 
  ungroup() %>% 
  group_by(season.year, alpha.code) %>% 
  summarise(p75.total.sbirds = quantile(total.sbirds, probs = c(0.75)),
            num.surveys = n()) %>% 
  ungroup() %>%
  bind_rows(., sbirds %>% 
              group_by(season.year, date) %>% 
              summarize(total.sbirds = sum(count)) %>% 
              ungroup() %>% 
              group_by(season.year) %>% 
              summarise(p75.total.sbirds = quantile(total.sbirds, probs = c(0.75)),
                        num.surveys = n()) %>%
              ungroup() %>% 
              mutate(alpha.code = "All")) %>% 
  separate(season.year, c("season", "year"), remove = F) %>% 
  mutate(year = as.numeric(year),
         giac.dummy = ifelse(year < 2009, 0, 1)) %>% 
  full_join(., rain) %>% 
  filter(season == "winter")


zmaxit = 400

fit_big_mod <- function(zspp) {
  sbirds_winter <- analysis_table %>% 
    filter(alpha.code == zspp)
  
  # big_mod <- glm.nb(floor(p75.total.sbirds) ~ poly(year, 2) * section + giac.dummy + seas.rain.mm + mean.bird.hunters, data = sbirds_winter, maxit = zmaxit)
  
  big_mod <- glm.nb(floor(p75.total.sbirds) ~ poly(year, 2) + seas.rain.mm + giac.dummy, data = sbirds_winter, maxit = zmaxit)
  
}

# ---
sbird_predicter_glm <- function(zmod, zmod.name){
  
  
  newdat <- analysis_table %>% 
    distinct(year) %>% 
    mutate(giac.dummy = ifelse(year < 2009, 0, 1),
           seas.rain.mm = 0)  
  
  ilink <- family(zmod)$linkinv
  all_best_pred = predict(zmod, newdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    cbind(newdat) %>% 
    ungroup() %>% 
    mutate(predicted = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit)),
           alpha.code = zmod.name)
}



zspp_list <- c("All", "DUNL", "WESA", "MAGO", "LESA", "SAND", "WILL", "DOSP", "BBPL", "BLTU", "YELL", "SEPL", "KILL", "SPSA")




sbird_mods <- map(zspp_list, fit_big_mod)

names(sbird_mods) <- zspp_list


sbird_preds <- map2_df(sbird_mods, names(sbird_mods), sbird_predicter_glm) %>% 
  left_join(analysis_table %>% dplyr::select(year, alpha.code, p75.total.sbirds))



saveRDS(sbird_preds, here("data/sbird_preds"))

# estimated effect size of rainfall and Giacomini restoration

#' Extract and backtransform model predictions and their CI from the Negative Binomial GLMs
#'
#' @param zmod 
#' @param zmod.name 
#'
#' @return
#' @export
#'
#' @examples
get_sbird_rain_preds_glmnb <- function(zmod, zmod.name) { 
  
#  znewdat = data.frame(year = floor(mean(analysis_table$year)),
#                       seas.rain.mm = seq(min(analysis_table$seas.rain.mm), max(analysis_table$seas.rain.mm), length.out = 10),
#                       giac.dummy = 0)
  
  
  #znewdat = data.frame(year = floor(mean(analysis_table$year)),
  #                     seas.rain.mm = c(mean(analysis_table$seas.rain.mm) - (0.5 * sd(analysis_table$seas.rain.mm)),
  #                                      mean(analysis_table$seas.rain.mm) + (0.5 * sd(analysis_table$seas.rain.mm))),
  #                     giac.dummy = 0)
  
  
  #znewdat = data.frame(year = floor(mean(analysis_table$year)),
  #                     seas.rain.mm = c(quantile(analysis_table$seas.rain.mm, 0.25),
  #                                      quantile(analysis_table$seas.rain.mm, 0.75)),
  #                     giac.dummy = 0)
  
  
  znewdat = data.frame(year = floor(mean(analysis_table$year)),
                       seas.rain.mm = c(-1, 0, 1),
                       giac.dummy = 0)
  
  
  ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    cbind(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit)),
           alpha.code = zmod.name)
}



sbird_rain_preds <- map2_df(sbird_mods, names(sbird_mods), get_sbird_rain_preds_glmnb)
saveRDS(sbird_rain_preds, here("data/sbird_rain_preds"))



#' Extract and backtransform model predictions and their CI from the Negative Binomial GLMs
#'
#' @param zmod 
#' @param zmod.name 
#'
#' @return
#' @export
#'
#' @examples
get_sbird_giac_preds_glmnb <- function(zmod, zmod.name) { 
  
  znewdat = data.frame(year = 2009,
                       seas.rain.mm = 0,
                       giac.dummy = c(0, 1))
  
  ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    cbind(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit)),
           alpha.code = zmod.name)
}



sbird_giac_preds <- map2_df(sbird_mods, names(sbird_mods), get_sbird_giac_preds_glmnb)
saveRDS(sbird_giac_preds, here("data/sbird_giac_preds"))

