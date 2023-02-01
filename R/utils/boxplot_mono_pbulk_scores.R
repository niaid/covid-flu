
# plotting functions for the boxplot scores for Extended Fig5d, 5e, Supplementary Figure 2e
# use for D0 only
plotcelltype_group <- function(df){
  my_comparisons1 <- list( c("HC", "COVR"))
  
  # covid vs healthy
  p <- ggplot(df, aes(x = group, y = module.score.mean)) +
    geom_boxplot(aes(fill=sex,color=sex), alpha = 0.4, width = 0.3, outlier.shape = NA)+
    geom_point(aes(fill=sex), shape = 21, color = "white", size = 4, position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1)) +
    scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
    scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
    stat_compare_means(comparisons = my_comparisons1, method = "wilcox.test", size = 8, label.y = 0.65)+
    facet_wrap(~geneset, scale = "free", ncol = 5) +
    ylim(-0.7, 0.85)+
    theme_bw(base_size = 26)

  return(p)  
}


# all timepoints, but covid only
plotcelltype_group_covid <- function(df){
  my_comparisons1 <- list(c("Day 0", "Day 1"), c("Day 0", "Day 28"))
  tmp <- filter(df, group == "COVR")
  tmp.HC <- filter(df, group == "HC", visit %in% c("Day 0"))
  
  medians <- plyr::ddply(tmp.HC, .variables = "geneset", summarise, median = median(module.score.mean), 
                         quantile.25 = quantile(module.score.mean, .25), quantile.75 = quantile(module.score.mean, .75))
  
  # covid vs healthy
  sex.color = ifelse(unique(tmp$sex) == "Female", "#BD3B29", "#0472B6")
  p <- ggplot(tmp, aes(x = visit, y = module.score.mean)) +
    geom_boxplot(aes(fill = sex), color="grey10", alpha = 0.4, width = 0.5, size = 1, outlier.shape = NA)+
    geom_point(aes(color=sex), size = 1.5)+
    scale_color_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
    scale_fill_manual(values = c("Female"="#BD3B29", "Male"="#0472B6"))+
    geom_line(aes(group = alt.subject.id), color = sex.color, alpha = 0.75)+
    geom_hline(data = medians, aes(yintercept = median), alpha = 0.5, linetype = "dashed", color = sex.color, size = 1.5)+
    stat_pvalue_manual(pval, label = "P.Value", tip.length = 0.01,size = 7)+
    facet_wrap(~geneset, ncol = 5) +
    ylim(min(tmp$module.score.mean)-0.1, max(tmp$module.score.mean)+0.3)+
    theme_bw(base_size = 25)

  return(p)  
}
