
# this script combines model estimates from the 3 monitoring projects
# end is 2 objects:
# 1. combined estimates showing overall trends
# 2. combined estimates showing effects of predictors 

library(tidyverse)
library(here)
library(birdnames)


custom_bird_list <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/my_R_general/birdnames_support/data/custom_bird_list")

options(scipen = 999)

# model coefficients for non-year variables ----
hep_coef <- readRDS(here("data/HEP_coef_ci")) %>% 
  separate(model, into = c("species", "subregion"), sep = "_") %>%
  filter(subregion == "TBAY") %>%
  rename("alpha.code" = species) %>% 
  mutate(alpha.code = ifelse(alpha.code == "All", "All.hep", alpha.code),
         common.name = ifelse(alpha.code == "All.hep", 
                              "Large herons and egrets", 
                              translate_bird_names(alpha.code, "alpha.code", "common.name")),
         common.name = ifelse(alpha.code == "DCCO", 
                              paste("Breeding", common.name), 
                              common.name),
         waterbird.group = ifelse(alpha.code == "DCCO", common.name, "Large herons and egrets")) %>% 
  select(-subregion)

sbird_coef <- readRDS(here("data/sbird_coef_ci")) %>%
  rename("alpha.code" = model) %>% 
  mutate(alpha.code = ifelse(alpha.code == "All", "All.sbird", alpha.code),
         common.name = case_when(alpha.code == "All.sbird" ~ "All shorebirds",
                                 alpha.code == "DOSP" ~ "Dowitcher spp.",
                                 alpha.code == "YELL" ~ "Yellowlegs spp.",
                                 TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name")),
         waterbird.group = "Shorebirds")

wbird_coef <- readRDS(here("data/wbird_coef_ci")) %>%
  rename("alpha.code" = model) %>% 
  mutate(alpha.code = ifelse(alpha.code == "ALL", "All.wbird", alpha.code),
         common.name = case_when(alpha.code == "All.wbird" ~ "All Water birds", 
                                 alpha.code == "SCAUP" ~ "Scaup spp.",
                                 TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name")),
         waterbird.group = "Water birds")




combined_predictor_coefs <- bind_rows(hep_coef, sbird_coef, wbird_coef)  %>% 
  filter(varb %in% c("subreg.rain", "seas.rain.mm", "giac.dummy", "fresh", "moci")) %>% 
  mutate(varb = case_when(grepl("rain|fresh", varb) ~ "rain",
                          grepl("giac", varb) ~ "giac",
                          TRUE ~ as.character(varb))) %>% 
  select(-coef, -lci, -uci) %>% 
  rename("predictor.varb" = varb,
         "estimate" = per.change.estimate,
         "lci" = per.change.lci,
         "uci" = per.change.uci) %>% 
  mutate(predictor.varb.label = case_when(predictor.varb == "rain" ~ "Response to rain",
                                          predictor.varb == "giac" ~ "Response to Giacomini restoration",
                                          predictor.varb == "moci" ~ "Response to ocean conditions"))

saveRDS(combined_predictor_coefs, here("data/combined_predictor_coefs"))

combined_predictor_coefs_up_down <- bind_rows(combined_predictor_coefs %>%
                                   mutate(predictor.label = case_when(predictor.varb == "rain" ~ "Wet year",
                                                                      predictor.varb == "giac" ~ "Post-\nrestoration",
                                                                      predictor.varb == "moci" ~ "Warm ocean/\nweak upwelling"),
                                          predictor.value = 1),
                                 combined_predictor_coefs %>%
                                   mutate(predictor.label = case_when(predictor.varb == "rain" ~ "Average rain year",
                                                                      predictor.varb == "giac" ~ "Pre-\nrestoration",
                                                                      predictor.varb == "moci" ~ "Average ocean\nconditions"),
                                          predictor.value = 0,
                                          lci = lci - estimate,
                                          uci = uci - estimate,
                                          estimate = estimate - estimate),
                                 combined_predictor_coefs %>%
                                   filter(predictor.varb != "giac") %>% 
                                   mutate(predictor.label = case_when(predictor.varb == "rain" ~ "Dry year",
                                                                      predictor.varb == "moci" ~ "Cool ocean/\nstrong upwelling"),
                                          predictor.value = -1,
                                          lci = lci - (2 * estimate),
                                          uci = uci - (2 * estimate),
                                          estimate = estimate - (2 * estimate))) %>% 
  mutate(predictor.label = factor(predictor.label, levels = c("Dry year", "Average rain year", "Wet year", "Pre-\nrestoration", "Post-\nrestoration", "Warm ocean/\nweak upwelling", "Average ocean\nconditions", "Cool ocean/\nstrong upwelling"))) %>% 
  arrange(waterbird.group, alpha.code, predictor.value)

saveRDS(combined_predictor_coefs_up_down, here("data/combined_predictor_coefs_up_down"))


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
                                          mutate(alpha.code = ifelse(alpha.code == "All", "All.wbird", alpha.code),
                                                 common.name = case_when(alpha.code == "All.wbird" ~ "All waterbirds", 
                                                                         alpha.code == "SCAUP" ~ "Scaup spp.",
                                                                         TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name")),
                                                 waterbird.group = "Water birds")
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
  mutate(common.name = case_when(alpha.code == "All" ~ "All Water birds",
                                 TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name"))) %>% 
  dplyr::select("year" = study.year, "predictor.value" = fresh, estimate, lci, uci, alpha.code, common.name) %>% 
  mutate(waterbird.group = "Water birds")


rain_preds <- bind_rows(hep_rain_preds, sbird_rain_preds, wbird_fresh_preds) %>% 
  mutate(predictor.varb = "rain",
         predictor.varb.label = "Response to rain",
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
         predictor.varb.label = "Response to Giacomini restoration",
         predictor.label = case_when(predictor.value == 0 ~ "Pre-\nrestoration",
                                     predictor.value == 1 ~ "Post-\nrestoration"),
         predictor.label = factor(predictor.label, levels = c("Pre-\nrestoration", "Post-\nrestoration")))

# moci effects
wbird_moci_preds <- readRDS(here("data/wbird_moci_preds")) %>% 
  mutate(common.name = case_when(alpha.code == "ALL" ~ "All water birds",
                                 alpha.code == "SCAUP" ~ "Scaup spp.",
                                 TRUE ~ translate_bird_names(alpha.code, "alpha.code", "common.name")))  %>% 
  dplyr::select("year" = study.year, "predictor.value" = moci, estimate, lci, uci, alpha.code, common.name) %>% 
  mutate(waterbird.group = "Water birds",
         predictor.varb = "moci",
         predictor.varb.label = "Response to ocean conditions",
         predictor.label = case_when(predictor.value == -1 ~ "Cold ocean/\nstrong upwelling",
                                     predictor.value == 0 ~ "Average ocean\nconditions",
                                     predictor.value == 1 ~ "Warm ocean/\nweak upwelling"),
         predictor.label = factor(predictor.label, levels = c("Cold ocean/\nstrong upwelling", "Average ocean\nconditions", "Warm ocean/\nweak upwelling")))



combined_predictor_estimates <- bind_rows(rain_preds, sbird_giac_preds, wbird_moci_preds) %>% 
  mutate(predictor.label = factor(predictor.label, levels = c("Dry year", "Average rain", "Wet year", 
                                                              "Pre-\nrestoration", "Post-\nrestoration", 
                                                              "Cold ocean/\nstrong upwelling", "Average ocean\nconditions", "Warm ocean/\nweak upwelling")))
  

saveRDS(combined_predictor_estimates, here("data/combined_predictor_estimates"))





