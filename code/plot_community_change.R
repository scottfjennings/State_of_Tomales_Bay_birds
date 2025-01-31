

# Load package
library(networkD3)

# Load energy projection data
URL <- "https://cdn.rawgit.com/christophergandrud/networkD3/master/JSONdata/energy.json"
Energy <- jsonlite::fromJSON(URL)


# Now we have 2 data frames: a 'links' data frame with 3 columns (from, to, value), and a 'nodes' data frame that gives the name of each node.
head( Energy$links )
head( Energy$nodes )


p <- sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
                   Target = "target", Value = "value", NodeID = "name",
                   units = "TWh", fontSize = 12, nodeWidth = 30)
p





library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(here)
library(birdnames)


custom_bird_list <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/my_R_general/birdnames_support/data/custom_bird_list")

source(here("code/utilities.R"))

options(scipen = 999)




start_end_trend_estimates <- readRDS(here("data/combined_trend_estimates")) %>% 
  filter(year >= 1992, !grepl("All", alpha.code)) %>% 
  group_by(waterbird.group) %>% 
  mutate(zyear = case_when(year == min(year) ~ "30 years ago", 
                           year == max(year) ~ "Now")) %>% 
  ungroup() %>% 
  filter(!is.na(zyear)) %>% 
  mutate(estimate = floor(estimate)) %>% 
  group_by(waterbird.group, zyear) %>% 
  mutate(tot.estimate.group = sum(estimate)) %>%
  ungroup() %>% 
  group_by(zyear) %>% 
  mutate(tot.estimate = sum(estimate),
         zyear = factor(zyear, levels = c("30 years ago", "Now"))) %>% 
  ungroup() %>% 
  arrange(zyear, alpha.code)

start_end_trend_estimates %>% 
  distinct(waterbird.group, zyear, tot.estimate.group) %>% 
  ggplot() +
  geom_col(aes(x = zyear, y = tot.estimate.group, fill = waterbird.group, group = waterbird.group))

start_end_trend_estimates %>% 
  distinct(waterbird.group, zyear, tot.estimate.group) %>% 
  arrange(zyear, desc(waterbird.group)) %>% 
  group_by(zyear) %>% 
  mutate(ztop = cumsum(tot.estimate.group),
         zbottom = ztop - tot.estimate.group) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(ztop, zbottom)) %>% 
  mutate(zorder = case_when(zyear == "30 years ago" & name == "ztop" ~ 1,
                            zyear == "30 years ago" & name == "zbottom" ~ 2,
                            zyear == "Now" & name == "zbottom" ~ 3,
                            zyear == "Now" & name == "ztop" ~ 4)) %>%
  arrange(waterbird.group, zorder) %>% 
  filter(waterbird.group != "Breeding Double-crested Cormorant") %>% 
  ggplot() +
  geom_polygon(aes(x = zyear, y = value, fill = waterbird.group, group = waterbird.group))+
  scale_y_continuous(breaks = seq(0, 40000, by = 5000), labels = seq(0, 40000, by = 5000)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme_bw(base_size = 18) +
  labs(x = "",
       y = "Total birds",
       fill = "")

ggsave(here("figures/total_bird_change.png"), width = 8)  

distinct(start_end_trend_estimates, zyear, tot.estimate)
36240-24832


  start_end_trend_estimates %>% 
    distinct(waterbird.group, zyear, tot.estimate.group) %>% 
  ggplot() +
    geom_point(aes(x = zyear, y = tot.estimate.group, color = waterbird.group)) +
    geom_line(aes(x = zyear, y = tot.estimate.group, color = waterbird.group, group = waterbird.group))
  