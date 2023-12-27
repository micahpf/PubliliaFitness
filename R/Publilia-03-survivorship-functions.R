## Publilia-03-survivorship: regression functions

### fit survivorship model
fit_surv <- function(data_in, survNom, survDen, tenure, treatment, quad){
  data <- data_in %>% filter(!is.na(.data[[tenure]]))
  
  if (!quad) {
    formula.mai <- paste0(survNom,"/",survDen,"~",tenure,"*",treatment) %>% as.formula()
    formula.ma <- paste0(survNom,"/",survDen,"~",tenure,"+",treatment) %>% as.formula()
  } else if (quad) {
    formula.mai <- paste0(
      survNom,"/",survDen,"~",
      tenure,"*",treatment,"+", "(I(as.numeric(",tenure,")^2)*",treatment,")"
    ) %>% as.formula()
    formula.ma <- paste0(
      survNom,"/",survDen,"~",
      tenure,"+",treatment,"+", "(I(as.numeric(",tenure,")^2)*",treatment,")"
    ) %>% as.formula()
  }
  
  mai.b <- glm(formula = formula.mai, data = data, weights = data[[survDen]], family = binomial)
  mai.qb <- glm(formula = formula.mai, data = data, weights = data[[survDen]], family = quasibinomial)
  ma.qb <- glm(formula = formula.ma, data = data, weights = data[[survDen]], family = binomial)
  models <- list(mai.b = mai.b, mai.qb = mai.qb, ma.qb = ma.qb)
  return(models)
}

capture_glm_stats <- function(model) {
  capture.output(summary(model))[-(6:which(grepl("Deviance Residuals: ", capture.output(summary(model))))-1)]
}

print_glm_stats <- function(models) {
  cat("Interaction, Binomial: \n")
  print(capture_glm_stats(models$mai.b))
  cat("\nInteraction, Quasibinomial: \n")
  print(capture_glm_stats(models$mai.b))
  cat("\nNo Interaction, Quasibinomial: \n")
  print(capture_glm_stats(models$mai.b))
}


### Create confidence intervals in the link scale for plotting, and backtransform for plot
# Based on: https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
ci_for_plot <- function(data_in, survNom, survDen, 
                        tenure, treatment, model){
  
  data <- data_in %>%
    filter(!is.na(.data[[tenure]])) %>%
    select(location, all_of(tenure), all_of(treatment), all_of(survDen), all_of(survNom)) %>%
    add_column(fit = predict(model, newdata = ., type = 'response'))
  
  ## grad the inverse link function
  ilink <- family(model)$linkinv
  
  ## add fit and se.fit on the **link** scale
  data_fit_link <- data %>% 
    add_column(
      predict(model, data, se.fit = TRUE)[1:2] %>% 
        as_tibble() %>%
        setNames(c('fit_link','se_link'))) %>%
    ## create the interval and backtransform
    mutate(fit_resp  = ilink(.data$fit_link),
           conf_upr = ilink(.data$fit_link + (2 * .data$se_link)),
           conf_lwr = ilink(.data$fit_link - (2 * .data$se_link)))
  
  return(data_fit_link)
}

# Example
# test_fits <- fit_surv(data_in=surv_reg_data, tenure="des", treatment="ant_treatment",
#                  survNom="n_hatched", survDen="med_t5eggs", quad=TRUE)
# 
# test_ci_data <- ci_for_plot(data_in=surv_reg_data, tenure="des", treatment="ant_treatment",
#                  survNom="n_hatched", survDen="med_t5eggs",
#                  model = test_fits$mai.qb)

### Plot model results with confidence intervals
plot_ci <- function(ci_data, survNom, survDen, treatment, tenure){
  p <- ggplot(data = ci_data, 
                 aes(x = .data[[tenure]], y = fit, 
                     color = .data[[treatment]], fill = .data[[treatment]])) +
    theme_publilia() +
    geom_point(aes(y = .data[[survNom]]/.data[[survDen]]), alpha = 0.25, size = 4, shape=16, key_glyph=draw_key_rect) + 
    geom_line(aes(y = fit), size = 1, key_glyph=draw_key_rect) +
    geom_ribbon(aes(ymin = conf_upr, ymax = conf_lwr),
                alpha = 0.5, color = NA, key_glyph=draw_key_rect)
  
  return(p)
}

# Example
# test_fits <- fit_surv(data_in=surv_reg_data, tenure="des", treatment="ant_treatment",
#                  survNom="n_hatched", survDen="med_t5eggs", quad=TRUE)
# 
# test_ci_data <- ci_for_plot(data_in=surv_reg_data, tenure="des", treatment="ant_treatment",
#                  survNom="n_hatched", survDen="med_t5eggs",
#                  model = test_fits$mai.qb)
# 
# (p <- plot_ci(ci_data = test_ci_data, tenure="des", treatment="ant_treatment",
#                  survNom="n_hatched", survDen="med_t5eggs") +
#     xlab("Latest desertion date relative to hatch") +
#     ylab("Hatching success") +
#     scale_color_manual(values = c("brown1", "black")) +
#     scale_fill_manual(values = c("brown1", "black")))
# summary(test_fits$mai.qb)

### Wrapper for regression fitting and plotting
runTenureXOffspringModel <- function(surv_reg_data, survNom, survDen, tenure, treatment, quad = TRUE) {
  
  if (quad) {
    test_fits <- surv_reg_data %>%
      fit_surv(data_in=., tenure=tenure, treatment=treatment, 
               survNom=survNom, survDen=survDen, quad = TRUE)
  } else if (!quad) {
    test_fits <- surv_reg_data %>%
      fit_surv(data_in=., tenure=tenure, treatment=treatment, 
               survNom=survNom, survDen=survDen, quad = FALSE)
  }
  
  test_ci_data <- ci_for_plot(data_in=surv_reg_data, tenure=tenure, treatment=treatment,
                              survNom=survNom, survDen=survDen,
                              model = test_fits$mai.qb)
  
  p <- plot_ci(ci_data = test_ci_data, tenure = tenure, treatment=treatment,
               survNom = survNom, survDen = survDen) +
    scale_color_manual(values = c("black", "brown1")) +
    scale_fill_manual(values = c("black", "brown1"))
  
  stat_summary <- summary(test_fits$mai.qb)
  
  obj_out <- list(
    fits = test_fits,
    ci = test_ci_data,
    p = p,
    summary = stat_summary
  )
  return(obj_out)
}

# Example
# des_hatch.desL.ant_treatment.quad <- surv_reg_data %>%
#   runTenureXOffspringModel(survNom="n_hatched", survDen="med_t5eggs", 
#                            tenure="des", treatment="ant_treatment")

## Regression stats interpretation
# function to get fits and confint
ci_for_stats <- function(data_in, model){
  
  data <- eval(substitute(
    data_in %>%
      add_column(fit = predict(model, newdata = ., type = 'response'))
  ))
  
  ## grad the inverse link function
  ilink <- family(model)$linkinv
  
  ## add fit and se.fit on the **link** scale
  data_fit_link <- eval(substitute(
    data %>% add_column(
      setNames(as_tibble(predict(model, data, se.fit = TRUE)[1:2]),
               c('fit_link','se_link'))) %>%
      ## create the interval and backtransform
      mutate(fit_resp  = ilink(fit_link),
             conf_upr = ilink(fit_link + (2 * se_link)),
             conf_lwr = ilink(fit_link - (2 * se_link)))
  ))
  
  return(data_fit_link)
}

# function to get max fit, plus upper and lower ci
get_max_fit <- function(data_in, model){
  
  max_fit <- eval(substitute(ci_for_stats(data_in=data_in, 
                                          model) %>%
                               group_by(ant_treatment) %>%
                               mutate(max = max(fit)) %>% filter(fit==max) %>% 
                               select(ant_treatment, des, max)))
  
  max_upr <- eval(substitute(ci_for_stats(data_in=data_in, 
                                          model) %>%
                               group_by(ant_treatment) %>%
                               mutate(max_upr = max(conf_upr)) %>% filter(conf_upr==max_upr) %>% 
                               select(ant_treatment, des_upr=des, max_upr)))
  
  max_lwr <- eval(substitute(ci_for_stats(data_in=data_in, 
                                          model) %>%
                               group_by(ant_treatment) %>%
                               mutate(max_lwr = max(conf_lwr)) %>% filter(conf_lwr==max_lwr) %>% 
                               select(ant_treatment, des_lwr=des, max_lwr)))
  
  data_out <- max_fit %>% left_join(max_upr) %>% left_join(max_lwr)
  return(data_out)
}

## Custom function for calculating qaicc (requires bbmle)
get_qaicc <- function(model.fits) {
  model.qaicc <- qAICc(model.fits$mai.b,
                       dispersion=summary(model.fits$mai.qb)$dispersion,
                       nobs = nrow(model.fits$mai.qb$data))
  return(model.qaicc)
}
