

add_study_year = function(df) {
  df <- df %>% 
    mutate(study.year = ifelse(month(date) > 6, year(date), year(date) - 1))
}



set_1month_winter_period <- function(df) {
  df <- mutate(year.day = lubridate::yday(date),
               winter.period = case_when(between(year.day, 304, 365) ~ "Nov-Dec", 
                                         between(year.day, 1, 31) ~ "Jan",
                                         between(year.day, 32, 60) ~ "Feb",
                                         TRUE ~ NA),
               winter.period = factor(winter.period, levels = c("Nov-Dec", "Jan", "Feb")))
}



set_2week_winter_period <- function(df) {
df <- mutate(year.day = yday(as.Date(date)),
             winter.period = case_when(between(year.day, 304, 320) ~ "Nov 1-15", 
                                       between(year.day, 321, 335) ~ "Nov 16-30", 
                                       between(year.day, 336, 350) ~ "Dec 1-15", 
                                       between(year.day, 351, 365) ~ "Dec 16-31", 
                                       between(year.day, 1, 15) ~ "Jan 1-15", 
                                       between(year.day, 16, 31) ~ "Jan 16-31",
                                       between(year.day, 32, 46) ~ "Feb 1-15",
                                       between(year.day, 47, 60) ~ "Feb 16-28(29)",
                                       TRUE ~ NA),
             winter.period = factor(winter.period, levels = c("Nov 1-15", "Nov 16-30", "Dec 1-15", "Dec 16-31", "Jan 1-15", "Jan 16-31", "Feb 1-15", "Feb 16-28(29)")))
}




#' trend_plotter
#' 
#' @description
#' plot estimated trend line for 1 or more species. Will facet_wrap if more than 1 species. Can optionally include raw data
#'
#' @param zspp character string for 4-letted species codes to plot
#' @param include.raw TRUE (default) to include points for raw data
#'
#' @return
#' @export
#'
#' @examples
trend_plotter <- function(zspp, include.raw = TRUE, include.spp.title = TRUE) {
  
  zpred <- filter(combined_trend_estimates, alpha.code %in% zspp)
  
  #ymax = ifelse(max(zpred$abundance) > max(zpred$uci), max(zpred$abundance), max(zpred$uci))
  
  zplot <- zpred %>% 
    ggplot(aes(x = study.year, y = mod.avg.pred)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL),alpha=0.3) +
    theme_bw(base_size = 18) +
    scale_y_continuous(breaks = scales::pretty_breaks(5), limits = c(0, NA)) +
    scale_x_continuous(breaks = scales::pretty_breaks(8), limits = c(1990, NA)) +
    xlab("Year") +
    ylab("Abundance")
  
  if(include.raw == TRUE) {
    zplot <- zplot   +
      geom_point(aes(x = study.year, y = p75.abund))
  }
  
  if(include.spp.title == TRUE) {
    zplot <- zplot   +
      labs(title = distinct(zpred, common.name))
  }
  
  
  return(zplot)
}


#' predictor_est_plotter
#' 
#' plot model estimates for non-year predictor variables
#'
#' @param pred_df data frame with model estimates
#'
#' @return a ggplot object
#'
#' @examples
predictor_est_plotter <- function(pred_df) {
  
  zhjust = -.25
  effects_plot <- pred_df  %>% 
    ggplot(aes(group = waterbird.group)) + 
    geom_line(aes(y = estimate, x = predictor.label)) +
    geom_ribbon(aes(x = predictor.label, ymin = lci, ymax = uci), alpha = 0.25) +
    theme_bw(base_size = 18) +
    guides(color = "none") +
    scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    labs(x = "",
         y = "% population change",
         fill = "",
         title = pred_df$predictor.varb.label) + 
    scale_x_discrete(expand=c(0.05, 0.15))
  return(effects_plot)
}


#' single_species_plotter
#' 
#' combine trend and predictor effect plots for a single species. is a wrapper for trend_plotter() and predictor_est_plotter()
#'
#' @param zalpha.code 4 letter code for the desired species
#'
#' @return a cowplot object, also saves the plot to figures/
#'
#' @examples
single_species_plotter <- function(zalpha.code) {
  
  common.name = filter(combined_trend_estimates, alpha.code == zalpha.code) %>% 
    distinct(common.name)
  
  zplot <- plot_grid(trend_plotter(zalpha.code) +
                       labs(title = paste(common.name, "\n\nPopulation change over time", sep = ""),
                            x = "",
                            y = "Abundance") + 
                       theme(strip.background = element_blank(),
                             strip.text.x = element_blank()),
                     combined_predictors %>%
                       filter(alpha.code == zalpha.code) %>% 
                       predictor_est_plotter() +
                       labs(title = "Effect of local variables") +
                       facet_wrap(~predictor.varb.label, nrow = 1, scales = "free"),
                     ncol = 1,
                     rel_heights = c(3, 3),
                     align = "v")
  
  
  zpath = as.character(here(paste("figures/", zalpha.code, "_combined.png", sep = "")))
  
  save_plot(zpath, zplot, dpi = 300, base_height = 8, base_width = 8)
  
  
  return(zplot)
  
}




#' extract model coefficients and their CI 
#'
#' @param zmod 
#' @param zmod.name 
#'
#' @return
#' @export
#'
#' @examples
#' coef_ci_log <- map2_df(spp_subreg_mods_log, names(spp_subreg_mods_log), get_coefs_cis)
#' coef_ci <- map2_df(spp_subreg_mods, names(spp_subreg_mods), get_coefs_cis)
#' coef_ci_glmnb <- map2_df(spp_subreg_mods_glmnb, names(spp_subreg_mods_glmnb), get_coefs_cis)
get_coefs_cis <- function(zmod, zmod.name) { 

  coefs <- coef(zmod) %>%
    data.frame() %>% 
    rename("coef" = 1) %>% 
    rownames_to_column("varb") 
  
  cis  <- confint(zmod) %>%
    data.frame() %>% 
    rename(lci = 1, uci = 2) %>% 
    rownames_to_column("varb")
  
  
  coef_ci <- full_join(coefs, cis) %>% 
    mutate(per.change.estimate = 100 * (exp(coef)-1),
           per.change.lci = 100 * (exp(lci)-1),
           per.change.uci = 100 * (exp(uci)-1),
           model = zmod.name) 
  
}







#' pred_plotter
#' 
#' generic plotting function for model estimates
#'
#' @param zestimates_df data frame with model estimates
#'
#' @return
#' @export
#'
#' @examples
pred_plotter <- function(zestimates_df) {
  min.x <- floor((min(zestimates_df$lci))/10) * 10
  max.x <- ceiling((max(zestimates_df$uci))/10) * 10
  
  zestimates_df <- zestimates_df%>% 
    mutate(zcolor = ifelse(estimate < 0, brewer.pal(11, "RdBu")[2], brewer.pal(11, "RdBu")[10]))   
  
  
  effects_plot <- zestimates_df  %>% 
    ggplot(group = c(common.name, predictor.label)) + 
    geom_segment(aes(x = -Inf, y = predictor.label, xend = mod.avg.pred, yend = predictor.label), color = 'gray85', size = 1) +
    geom_segment(aes(x = -Inf, y = predictor.label, xend = lower.CL, yend = predictor.label), size = 1, color = 'gray50') +
    geom_point(aes(x = estimate, y = predictor.label), color = zestimates_df$zcolor, size = 4) +
    geom_segment(aes(x = estimate, y = predictor.label, xend = 0, yend = predictor.label), color = zestimates_df$zcolor, size = 1) +
    geom_errorbar(aes(y = predictor.label, xmin = lci, xmax = uci), size = 2, color = zestimates_df$zcolor, width = 0.5) +
    geom_vline(xintercept = 0)  +
    theme_bw() +
    guides(color = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    labs(y = "",
         x = "",
         fill = "")
  return(effects_plot)
}





#' coef_plotter_horiz
#' 
#' plot model coefficients for local predictor variables oriented horizontally. not currently used
#'
#' @param coef_df 
#'
#' @return
#' @export
#'
#' @examples
coef_plotter_horiz <- function(coef_df) {
  effects_plot <- coef_df  %>% 
    ggplot() + 
    geom_vline(xintercept = 0, color = "gray50")  +
    geom_point(aes(x = per.change.estimate, y = predictor.label), color = coef_df$zcolor, size = 4) +
    geom_errorbar(aes(y = predictor.label, xmin = per.change.lci, xmax = per.change.uci), size = 1, color = coef_df$zcolor, width = 0.3) +
    theme_bw() +
    guides(color = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    labs(y = "",
         x = "% population change",
         fill = "",
         title = coef_df$predictor.varb.label) +
    scale_y_discrete(limits=rev) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  if("giac" %in% coef_df$varb){
    if(filter(coef_df, varb == "giac")$per.change.estimate < 0) {
      effects_plot <- effects_plot +
        geom_text(aes(x = per.change.estimate, y = predictor.label, label = lab.text,
                      vjust = ifelse(varb == "giac", 2 * predictor.label.vjust, 0), 
                      hjust = ifelse(varb == "giac", NA, 1.9)), color = coef_df$zcolor)
    } else {
      effects_plot <- effects_plot +
        geom_text(aes(x = per.change.estimate, y = predictor.label, label = lab.text, 
                      vjust = ifelse(varb == "giac", 2 * predictor.label.vjust, 0), 
                      hjust = ifelse(varb == "giac", NA, -0.5)), color = coef_df$zcolor)
    }
  } else if("moci" %in% coef_df$varb) { 
    xx <- filter(coef_df, varb == "moci", per.change.estimate > 0)
    
    if(xx$per.change.estimate * 2 > xx$uci)  {
      effects_plot <- effects_plot +
        geom_text(aes(x = ifelse(varb == "moci", per.change.estimate * 0.7, per.change.estimate), y = predictor.label, label = lab.text, 
                      vjust = 2 * predictor.label.vjust), color = coef_df$zcolor)
    } else {
      effects_plot <- effects_plot +
        geom_text(aes(x = per.change.estimate, y = predictor.label, label = lab.text,
                      vjust = 2 * predictor.label.vjust), color = coef_df$zcolor)
    }
  } else { 
    effects_plot <- effects_plot +
      geom_text(aes(x = per.change.estimate, y = predictor.label, label = lab.text, vjust = 2 * predictor.label.vjust), color = coef_df$zcolor)
  }
  
  return(effects_plot)
}



#' coef_plotter
#' 
#' plot model coefficients for local predictor variables. not currently used
#'
#' @param coef_df 
#'
#' @return
#' @export
#'
#' @examples
coef_plotter <- function(coef_df) {
  ann_df <- filter(coef_df, alpha.code == distinct(coef_df, alpha.code)$alpha.code[[1]])
  
  zhjust = -.25
  effects_plot <- coef_df  %>% 
    ggplot() + 
    geom_hline(yintercept = 0, color = "gray50")  +
    geom_point(aes(y = per.change.estimate, x = predictor.varb.label), color = coef_df$zcolor, size = 4, position = position_nudge(x = zhjust)) +
    geom_errorbar(aes(x = predictor.varb.label, ymin = per.change.lci, ymax = per.change.uci), size = 1, color = coef_df$zcolor, width = 0.1, position = position_nudge(x = zhjust)) +
    theme_bw() +
    guides(color = "none") +
    scale_y_continuous(breaks = scales::pretty_breaks(8)) +
    labs(x = "",
         y = "% population change",
         fill = "",
         title = "") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_text(aes(y = per.change.estimate, x = predictor.varb.label, label = trimws(lab.text), hjust = -0.1), color = coef_df$zcolor)
}



