## Publilia-05-fitness custom functions

# Confidence interval for estimating fitness
ci_for_fitness <- function(data_in, survNom, survDen, 
                           tenure, treatment, vis, model){
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

## Split violin plot
# Based on jan-glx answer on https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

## A function to run Kruskall-Wallis test and Dunn tests for each panel
fitness_KW_test <- function(data, y, group){
  # non-parametric: Kruskall-Wallis test followed by Dunn test
  # run KW test using rstatix::kruskal_test
  KWtest <- eval(substitute(data %>% kruskal_test(y ~ group)))
  # get effect size using effectsize::rank_epsilon_squared
  KWepsilon2 <- eval(substitute(rank_epsilon_squared(y ~ group, data=data)))
  # run Dunn test using rstatix::dunn_test
  dunn <- eval(substitute(data %>% dunn_test(y ~ group)))
  # get multiple comparison groups letters using multcompView::multcompLetters
  dunn.vec <- dunn$p.adj.signif
  names(dunn.vec) <- paste(dunn$group1,dunn$group2,sep="-")
  multcompLetters <- tibble(x=names(multcompLetters(dunn.vec)$Letters),
                            label=unname(multcompLetters(dunn.vec)$Letters))
  # return
  out <- list(KWtest=KWtest, KWepsilon2=KWepsilon2, dunn=dunn, multcompLetters=multcompLetters)
  return(out)
}

# Example
# fitness_life_strats_data <- est_fitness_life_strats %>% 
#   mutate(desert_cat_ants = ifelse(ants_removed, 
#                                   paste(desert_cat,"ants removed"),paste(desert_cat,"control")))
# fitness_KW_test(data=fitness_life_strats_data, y=est_eggs_life, group=desert_cat_ants)


## plot boxplots, including KWtest letters
fitness_plot_dunn <- function(data, y, group, itero=F, multcompLetters, letterSize, ymax=160){
  if(itero==F){
    (p <- eval(substitute(data %>%
                            filter(!females_removed, !is.na(group)) %>%
                            ggplot(data = ., aes(x = group, y=y, fill=ants_removed, color=ants_removed)) + 
                            ylim(-ymax * 0.05, ymax) + 
                            geom_split_violin(color = NA, alpha = 0.5, key_glyph = draw_key_rect) + 
                            geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.3), 
                                       pch = 16, size = 2, key_glyph = draw_key_rect) + 
                            xlab("Guarding strategy for first clutch") + theme_publilia() + 
                            scale_fill_manual(labels = c("Control", "Ant exclusion"), values = c("brown1", "black")) + 
                            scale_color_manual(labels = c("Control", "Ant exclusion"),  values = c("brown1", "black")) +
                            geom_text(data = multcompLetters, 
                                      aes(x = x, y = -ymax * 0.05, label = label), 
                                      fontface = "bold", size = letterSize, position=position_dodge(width=0.5)))))
    return(p)
  } else {
    # single-clutch females
    (p_semel <- eval(substitute(data %>%
                                  filter(!females_removed, !is.na(group),
                                         (n_clutches_with_eggs_iv==1)) %>%
                                  ggplot(data = ., aes(x = group, y=y, fill=ants_removed, color=ants_removed)) + 
                                  ylim(-ymax * 0.05, ymax) + 
                                  geom_split_violin(color = NA, alpha = 0.5, key_glyph = draw_key_rect) + 
                                  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.3), 
                                             pch = 16, size = 2, key_glyph = draw_key_rect) + 
                                  xlab("Guarding strategy for first clutch") + theme_publilia() + 
                                  scale_fill_manual(labels = c("Control", "Ant exclusion"), values = c("brown1", "black")) + 
                                  scale_color_manual(labels = c("Control", "Ant exclusion"),  values = c("brown1", "black")) +
                                  geom_text(data = multcompLetters, 
                                            aes(x = x, y = -ymax * 0.05, label = label), 
                                            fontface = "bold", size = letterSize, position=position_dodge(width=0.5)))))
    # multi-clutch females
    (p_itero <- eval(substitute(data %>%
                                  filter(!females_removed, !is.na(group),
                                         (n_clutches_with_eggs_iv>1)) %>%
                                  ggplot(data = ., aes(x = group, y=y, fill=ants_removed, color=ants_removed)) + 
                                  ylim(-ymax * 0.05, ymax) + 
                                  geom_split_violin(color = NA, alpha = 0.5, key_glyph = draw_key_rect) + 
                                  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.3), 
                                             pch = 16, size = 2, key_glyph = draw_key_rect) + 
                                  xlab("Guarding strategy for first clutch") + theme_publilia() + 
                                  scale_fill_manual(labels = c("Control", "Ant exclusion"), values = c("brown1", "black")) + 
                                  scale_color_manual(labels = c("Control", "Ant exclusion"),  values = c("brown1", "black")) +
                                  geom_text(data = multcompLetters, 
                                            aes(x = x, y = -ymax * 0.05, label = label), 
                                            fontface = "bold", size = letterSize, position=position_dodge(width=0.5)))))
    p_list<-list(semel=p_semel, itero=p_itero)
    return(p_list)
  }
}

# Example
# fitness_life_strats_data <- est_fitness_life_strats %>%
#   mutate(desert_cat_ants = ifelse(ants_removed,
#                                   paste(desert_cat,"ants removed"),paste(desert_cat,"control")))
# eggs_fitness_stats <- fitness_KW_test(data=fitness_life_strats_data, y=est_eggs_life, group=desert_cat_ants)
# 
# fitness_plot_dunn(data=fitness_life_strats_data, y=est_eggs_life, group=desert_cat_ants,
#                   itero=F, multcompLetters=eggs_fitness_stats$multcompLetters,
#                   letterSize=8)



