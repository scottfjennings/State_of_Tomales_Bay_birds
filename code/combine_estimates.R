
# this script combines model estimates from the 3 monitoring projects
# end is 2 objects:
# 1. combined estimates showing overall trends
# 2. combined estimates showing effects of predictors 

library(tidyverse)
library(here)
library(birdnames)


custom_bird_list <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/my_R_general/birdnames_support/data/custom_bird_list")

options(scipen = 999)


# overall trends ----
hep_preds <- readRDS(here("data/HEP_predictions"))
sbird_preds <- readRDS(here("data/sbird_preds"))
wbird_preds <- readRDS(here("data/wbird_preds"))


combined_trend_estimates <- bind_rows(hep_preds %>%
                                          filter(subregion == "TBAY") %>%
                                          rename("alpha.code" = species,
                                                 "abundance" = tot.nests) %>% 
                                          mutate(alpha.code = ifelse(alpha.code == "All", "All.hep", alpha.code),
                                                 common.name = ifelse(alpha.code == "All.hep", "Large herons and egrets", translate_bird_names(alpha.code, "alpha.code", "common.name")),
                                                 common.name = ifelse(alpha.code == "DCCO", paste("Breeding", common.name), common.name),
                                                 waterbird.group = ifelse(alpha.code == "DCCO", common.name, "Large herons and egrets")),
                                        sbird_preds %>% 
                                          rename("estimate" = predicted,
                                                 "abundance" = p75.total.sbirds) %>% 
                                          mutate(alpha.code = ifelse(alpha.code == "All", "All.sbird", alpha.code),
                                                 common.name = case_when(alpha.code == "All.sbird" ~ "All shorebirds",
                                                                         alpha.code == "DOSP" ~ "Dowitcher spp.",
                                                                         alpha.code == "YELL" ~ "Yellowlegs spp.",
                                                                         TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name")),
                                                 waterbird.group = "Shorebirds"),
                                        wbird_preds %>%
                                          rename("year" = study.year,
                                                 "abundance" = p75.abund) %>% 
                                          mutate(alpha.code = ifelse(alpha.code == "ALL", "All.wbird", alpha.code),
                                                 common.name = case_when(alpha.code == "All.wbird" ~ "All open water birds", 
                                                                         alpha.code == "SCAUP" ~ "Scaup spp.",
                                                                         TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name")),
                                                 waterbird.group = "Open water birds")
) %>% 
  select(year, waterbird.group, alpha.code, common.name, abundance, estimate, lci, uci) %>% 
  group_by(common.name) %>% 
  mutate(upper.y = 1.1 * max(uci, na.rm = TRUE)) %>% 
  ungroup()


saveRDS(combined_trend_estimates, here("data/combined_trend_estimates"))


# predictor effects ----
# rain effects
hep_rain_preds <- readRDS(here("data/hep_rain_preds")) %>% 
  filter(subregion == "TBAY", species %in% c("GBHE", "GREG", "All")) %>% 
  mutate(common.name = translate_bird_names(species, "alpha.code", "common.name")) %>% 
  dplyr::select(year, "predictor.value" = subreg.rain, estimate, lci, uci, "alpha.code" = species, common.name) %>% 
  mutate(waterbird.group = "Large herons and egrets")



sbird_rain_preds <- readRDS(here("data/sbird_rain_preds")) %>% 
  mutate(common.name = case_when(alpha.code == "All" ~ "All shorebirds",
                                 alpha.code == "DOSP" ~ "Dowitcher spp.",
                                 alpha.code == "YELL" ~ "Yellowlegs spp.",
                                 TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name"))) %>% 
  dplyr::select(year, "predictor.value" = seas.rain.mm, estimate, lci, uci, alpha.code, common.name) %>% 
  mutate(waterbird.group = "Shorebirds")


wbird_fresh_preds <- readRDS(here("data/wbird_fresh_preds")) %>% 
  mutate(common.name = case_when(alpha.code == "All" ~ "All open water birds",
                                 TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name"))) %>% 
  dplyr::select("year" = study.year, "predictor.value" = fresh, estimate, lci, uci, alpha.code, common.name) %>% 
  mutate(waterbird.group = "Open water birds")


rain_preds <- bind_rows(hep_rain_preds, sbird_rain_preds, wbird_fresh_preds) %>% 
  filter(!alpha.code %in% c("All", "ALL")) %>% 
  mutate(predictor.varb = "rain",
         predictor.label = case_when(predictor.value == -1 ~ "Dry year",
                                     predictor.value == 0 ~ "Average rain",
                                     predictor.value == 1 ~ "Wet year"),
         predictor.label = factor(predictor.label, levels = c("Dry year", "Average rain", "Wet year")))

# Gicaomini effects
sbird_giac_preds <- readRDS(here("data/sbird_giac_preds")) %>% 
  mutate(common.name = case_when(alpha.code == "All" ~ "All shorebirds",
                                 alpha.code == "DOSP" ~ "Dowitcher spp.",
                                 alpha.code == "YELL" ~ "Yellowlegs spp.",
                                 TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name"))) %>% 
dplyr::select(year, "predictor.value" = giac.dummy, estimate, lci, uci, alpha.code, common.name) %>% 
  mutate(waterbird.group = "Shorebirds",
         predictor.varb = "giac",
         predictor.label = case_when(predictor.value == 0 ~ "Pre-restoration",
                                     predictor.value == 1 ~ "Post-restoration"),
         predictor.label = factor(predictor.label, levels = c("Pre-restoration", "Post-restoration")))

# moci effects
wbird_moci_preds <- readRDS(here("data/wbird_moci_preds")) %>% 
  mutate(common.name = case_when(alpha.code == "ALL" ~ "All open water birds",
                                 alpha.code == "SCAUP" ~ "Scaup spp.",
                                 TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name")))  %>% 
  dplyr::select("year" = study.year, "predictor.value" = moci, estimate, lci, uci, alpha.code, common.name) %>% 
  mutate(waterbird.group = "Open water birds",
         predictor.varb = "moci",
         predictor.label = case_when(predictor.value == -1 ~ "Cold ocean/strong upwelling",
                                     predictor.value == 0 ~ "Average ocean conditions",
                                     predictor.value == 1 ~ "Warm ocean/weak upwelling"),
         predictor.label = factor(predictor.label, levels = c("Cold ocean/strong upwelling", "Average ocean conditions", "Warm ocean/weak upwelling")))



all_preds_base <- bind_rows(rain_preds, sbird_giac_preds, wbird_moci_preds)


combined_predictor_estimates <- all_preds_base %>% 
  filter(predictor.value == 0) %>% 
  select("base.estimate" = estimate, alpha.code, common.name, waterbird.group, predictor.varb) %>% 
  full_join(all_preds) %>% 
  mutate(across(c(estimate, lci, uci), ~. - base.estimate))


saveRDS(combined_predictor_estimates, here("data/combined_predictor_estimates"))


