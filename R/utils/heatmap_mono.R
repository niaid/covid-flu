
### heatmap function for plotting reset genes in monocytes and labeling interested genes
subject.hm.fm <- function(exprs_mtx=exprs_mtx, meta=meta, celltype="Mono_Classical", module="M4.0.M11.0.f.m.union", 
                              labelset=labelset, rowsplitLE){
  labelLE <- rownames(exprs_mtx) %in% labelset
  names(labelLE) <- rownames(exprs_mtx)
  
  # set colors
  col_fun = colorRamp2(c(-1.5, 0, 1.6), c("#1470F7", "#000000", "#FFFF00")) # blue to yellow
  group.color <- c("HC" = "#e5862d", "COVR" = "#707071")
  sex.color <- c("Female"="#BD3B29", "Male"="#0472B6")
  
  # set how many genes to label, due to space limit
  label <- c(which(rownames(exprs_mtx) %in% labelset))
  label_genes <- rownames(exprs_mtx)[label]
  # set heatmap annotation
  ha = HeatmapAnnotation(
    group = meta$group,
    # sex = meta$sex,
    col = list(group = group.color)
  )
  
  # plot
  p <- Heatmap(pheatmap:::scale_rows(exprs_mtx), name = paste(celltype, module, sep = "_"),
               show_column_names = FALSE, 
               top_annotation = ha, cluster_columns = FALSE,
               column_split = meta$Timepoint.group,
               heatmap_legend_param = list(title = 'LE.reset.M4.M11'),
               # row_split = coLE,
               col = col_fun,
               # row_names_gp = gpar(fontsize = 8),
               column_gap = unit(2, "mm")
  )+
    rowAnnotation(foo = anno_mark(at = label,
                                  labels = label_genes,
                                  labels_gp = gpar(fontsize = 11)))
  return(p)
}



### heatmap function for showing all genes and also highlight specific sets
### heatmap function for plotting acute COVID monocyte signature in  baseline samples
subject.hm <- function(exprs_mtx, meta, labelset, celltype, module, rowsplitLE){
  coLE <- rownames(exprs_mtx) %in% rowsplitLE
  names(coLE) <- rownames(exprs_mtx)
  
  # set colors
  col_fun = colorRamp2(c(-1.5, 0, 1.6), c("#1470F7", "#000000", "#FFFF00")) # blue to yellow
  group.color <- c("HC" = "#e5862d", "COVR" = "#707071")
  sex.color <- c("Female"="#BD3B29", "Male"="#0472B6")
  
  label <- c(which(rownames(exprs_mtx) %in% labelset))
  label_genes <- rownames(exprs_mtx)[label]
  # set heatmap annotation
  ha = HeatmapAnnotation(
    group = meta$group,
    sex = meta$sex,
    col = list(group = group.color,
               sex = sex.color)
  )
  
  # plot
  p <- Heatmap(pheatmap:::scale_rows(exprs_mtx), name = paste(celltype, module, sep = "_"),
               show_column_names = FALSE,
               show_row_names = TRUE,
               top_annotation = ha, 
               cluster_columns = FALSE,
               # cluster_rows = FALSE,
               column_split = meta$Timepoint.group,
               row_split = coLE,
               # row_km = 2,
               col = col_fun,
               row_names_gp = gpar(fontsize = 8),
               column_gap = unit(2, "mm")
  )
  rowAnnotation(foo = anno_mark(at = label,
                                labels = label_genes,
                                labels_gp = gpar(fontsize = 8)))
  
  return(p)
}

