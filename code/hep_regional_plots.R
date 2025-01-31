

library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(here)
library(birdnames)


custom_bird_list <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/my_R_general/birdnames_support/data/custom_bird_list")

source(here("code/utilities.R"))

options(scipen = 999)


hep_preds <- readRDS(here("data/HEP_predictions"))


combined_predictors <- readRDS(here("data/combined_predictor_coefs_up_down"))


# We also have heron and egret monitoring data for the rest of the north SF Bay Area. How does Tomales Bay compare to trends across larger areas? We can summarize our monitoring data to show all nesting colonies along the outer Pacific coast of Marin and Sonoma Counties, and all colonies in the entire study area. For Great Blue Heron, 
hep_preds %>% 
  filter(!subregion %in% c("OUCnoTBAY"), !is.na(subreg.name), species %in% c("All")) %>%
  mutate(subreg.name = ifelse(subregion == "All", "Entire SF Bay area", as.character(subreg.name))) %>% 
  ggplot() +
  geom_line(aes(x = year, y = estimate)) +
  geom_ribbon(aes(x = year, ymin = lci, ymax = uci), alpha = 0.25) +
  geom_point(aes(x = year, y = tot.nests))  + 
  facet_wrap(~subreg.name, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(8), limits = c(1990, NA)) +
  labs(y = "Abundance",
       x = "Year",
       title = "Large heron and egret regional trends",
       color = "",
       fill = "") +
  theme_bw()  +
  theme_bw(base_size = 14)

+
  theme(legend.position=c(.23,.75),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))

ggsave(here("figures/HEP_tbay_ouc_all_trends.png"), dpi = 300, width = 8)
