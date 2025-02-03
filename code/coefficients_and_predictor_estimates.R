




########## shorebirds

# model coeficients and 95% CI ----
# DOSP model throwing error, need to get CI manually from SE
dosp.coef <- coef(sbird_mods$DOSP) %>%
  data.frame() %>% 
  rename("coef" = 1) %>% 
  rownames_to_column("varb") 

dosp_cis  <- summary(sbird_mods$DOSP) %>%
  coefficients() %>% 
  data.frame() %>% 
  rownames_to_column("varb") %>% 
  mutate(lci = Estimate - (1.96 * Std..Error),
         uci = Estimate + (1.96 * Std..Error)) %>% 
  dplyr::select(varb, "coef" = Estimate, lci, uci)

dosp_coef_cis <- full_join(dosp.coef, dosp_cis) %>% 
  mutate(per.change.estimate = 100 * (exp(coef)-1),
         per.change.lci = 100 * (exp(lci)-1),
         per.change.uci = 100 * (exp(uci)-1),
         model = "DOSP")

coef_ci <- map2_df(sbird_mods[-grep("DOSP", names(sbird_mods))], names(sbird_mods[-grep("DOSP", names(sbird_mods))]), get_coefs_cis) %>% 
  bind_rows(dosp_coef_cis)
saveRDS(coef_ci, "data/sbird_coef_ci")



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




############## waterbirds


#' get_wbird_fresh_preds
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
get_wbird_fresh_preds <- function(zmod, zmod.name) {
  
  ana_table <- spp_annual_full_preds %>% 
    filter(alpha.code == zmod.name)
  
  #  znewdat = data.frame(study.year = mean(ana_table$study.year),
  #                       fresh = seq(min(ana_table$fresh), max(ana_table$fresh), length.out = 10),
  #                       moci = mean(ana_table$moci))
  
  
  #znewdat = data.frame(study.year = floor(mean(ana_table$study.year)),
  #                     fresh = c(mean(ana_table$fresh) - (0.5 * sd(ana_table$fresh)),
  #                               mean(ana_table$fresh) + (0.5 * sd(ana_table$fresh))),
  #                     moci = mean(ana_table$moci))
  
  #znewdat = data.frame(study.year = floor(mean(ana_table$study.year)),
  #                     fresh = c(quantile(ana_table$fresh, 0.25),
  #                               quantile(ana_table$fresh, 0.75)),
  #                     moci = mean(ana_table$moci))
  
  znewdat = data.frame(study.year = floor(mean(ana_table$study.year)),
                       fresh = c(-1, 0, 1),
                       moci = 0)
  
  
  ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    bind_cols(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit)),
           alpha.code = zmod.name)
  
}


wbird_fresh_preds <- map2_df(all_spp_mod, names(all_spp_mod), get_wbird_fresh_preds)

saveRDS(wbird_fresh_preds, here("data/wbird_fresh_preds"))


#' get_wbird_moci_preds
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
get_wbird_moci_preds <- function(zmod, zmod.name) {
  
  ana_table <- spp_annual_full_preds %>% 
    filter(alpha.code == zmod.name)
  
  # zmoci = seq(min(ana_table$moci), max(ana_table$moci), length.out = 10)
  
  #zmoci = seq(quantile(ana_table$moci, 0.25), quantile(ana_table$moci, 0.75))
  
  zmoci = c(-1, 0, 1)
  
  znewdat = data.frame(study.year = floor(mean(ana_table$study.year)),
                       fresh = 0,
                       moci = zmoci)
  
  
  
  ilink <- family(zmod)$linkinv
  best_pred = predict(zmod, znewdat, se.fit=TRUE, type='link') %>% 
    data.frame() %>% 
    bind_cols(znewdat) %>% 
    ungroup() %>% 
    mutate(estimate = ilink(fit),
           lci = ilink(fit - (1.96 * se.fit)),
           uci = ilink(fit + (1.96 * se.fit)),
           alpha.code = zmod.name) 
  
}


wbird_moci_preds <- map2_df(all_spp_mod, names(all_spp_mod), get_wbird_moci_preds)

saveRDS(wbird_moci_preds, here("data/wbird_moci_preds"))

