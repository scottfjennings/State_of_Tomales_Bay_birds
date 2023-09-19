


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
    geom_segment(aes(x = -Inf, y = predictor.label, xend = estimate, yend = predictor.label), color = 'gray85', size = 1) +
    geom_segment(aes(x = -Inf, y = predictor.label, xend = lci, yend = predictor.label), size = 1, color = 'gray50') +
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


