


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
