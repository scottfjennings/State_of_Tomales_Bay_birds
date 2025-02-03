# prepare predictor variables 


library(tidyverse)
library(lubridate)
library(here)




# moci ----
# We use the annual mean value of the North section MOCI values for the AMJ, JAS, OND periods because:
# 1. transport of nutrients is more likely to have a downcoast component than an upcoast component, thus ocean conditions fro Pt Reyes North are more likely to influece TB than from Pt Reyes south. 
# 2. we lack conclusive prior research suggesting MOCI effect would be strongest at a particular lag, so averaging the prior year's values takes a broad, exploratory view of this variable's importance

# 3. Conditions during JFM may act on current year and next year. so to avoid complicated interpretation just exclude that quarter

mean_north_moci <- read.csv("https://www.faralloninstitute.org/s/CaliforniaMOCI.csv") %>% 
  dplyr::select(study.year = Year, season = Season, moci = North.California..38.42N.) %>% 
  filter(season != "JFM") %>% 
  group_by(study.year) %>% 
  summarise(moci = mean(moci)) %>% 
  ungroup() %>% 
  mutate(moci = scale(moci, scale = TRUE, center = TRUE)[,1])
#

# freshwater inflow ----
# Although Walker Creek flow is meaningful, the Lagunitas and Walker Creek flows are highly correlated (Pearson R = 0.95) but the Lagunitas Flow has higher spikes and better represents the combined freshwater flow into the bay. Using just the Lagunitas flow in this analysis to simplify   
lacr <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/general_data_sources/rivers/data/flow_data/lagunitas_daily_discharge") %>% 
  dplyr::select(date = Date, mean.daily.cfs = Flow) %>% 
  mutate(study.year = ifelse(month(date) < 3, year(date) - 1, year(date))) %>% 
  filter(study.year > 1991) %>% 
  data.frame() %>% 
  mutate(flow.start = paste(study.year, "-10-01", sep = ""),
         flow.end = paste(study.year + 1, "-02-15", sep = "")) %>% 
  filter(date >= flow.start & date <= flow.end) %>% 
  group_by(study.year) %>% 
  summarise(flow = mean(mean.daily.cfs)) %>% 
  mutate(flow = scale(flow, scale = TRUE, center = TRUE)[,1])

# ggplot(annual_lacr) + geom_line(aes(x = study.year, y = flow))
 
# PRISM rainfall for shorebirds
# rain from C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/general_data_sources/Rainfall/
rain <- read.csv("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/general_data_sources/Rainfall/data/derived_data/tomales_mean_month_rain.csv") %>% 
  mutate(study.year = ifelse(month < 7, year - 1, year)) %>% 
  filter(month <= 2 | month > 7) %>% 
  group_by(study.year) %>% 
  summarise(seas.rain.mm = sum(tomales.rain.mm),
            seas.rain.mm = round(seas.rain.mm)) %>% 
  ungroup() %>% 
  filter(study.year > 1988) %>% 
  mutate(seas.rain.mm = scale(seas.rain.mm, scale = TRUE, center = TRUE)[,1])
