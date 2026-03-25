########################################################################
# Self-contained analysis script for Serafini RGS4 mouse DRG study
# Starts from srt.rds produced by the scrna_pipeline
# Does NOT run SCENIC, cNMF or WGCNA
########################################################################

### ===== 1. LOAD LIBRARIES ================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(gplots)
  library(dplyr)
  library(Matrix)
  library(ggpubr)
  library(patchwork)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(fgsea)
  library(tidyr)
  library(harmony)
  library(viridis)
  library(BiocParallel)
  library(ggridges)
  library(clusterProfiler)
  library(tidyverse)
  library(rstatix)
  library(paletteer)
  library(circlize)
  library(GO.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
})

set.seed(1234)


### ===== 2. HELPER FUNCTIONS =============================================

## ggplot theme (from ggplot_aestetics.R) --------------------------------
gtheme <- theme(
  axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1),
  axis.line    = element_line(colour = 'black', size = .1),
  axis.ticks   = element_line(colour = "black", size = .1),
  panel.background = element_blank()
)

gtheme_no_rot <- theme(
  axis.line    = element_line(colour = 'black', size = .1),
  axis.ticks   = element_line(colour = "black", size = .1),
  panel.background = element_blank()
)

## ModScoreCor (from useful_functions.R) ----------------------------------
ModScoreCor <- function(seurat_obj, geneset_list, listName,
                        cor_threshold = NULL, pos_threshold = .1, outdir) {
  require(Seurat)
  message('Run AddModuleScore')
  seurat_obj <- AddModuleScore(seurat_obj, geneset_list)
  seurat_obj@meta.data <- seurat_obj@meta.data[,
    !colnames(seurat_obj@meta.data) %in% names(geneset_list)]
  colnames(seurat_obj@meta.data)[
    colnames(seurat_obj@meta.data) %in%
      paste0('Cluster', seq_along(geneset_list))] <- names(geneset_list)
  message(paste('Annotate cells based on highest module score and store in column:',
                paste0(listName, '_r', cor_threshold, '_max')))
  if (length(geneset_list) == 1) {
    seurat_obj@meta.data[, paste0(listName, '_r', cor_threshold, '_max')] <-
      ifelse(seurat_obj@meta.data[, names(geneset_list)] > pos_threshold,
             'pos', 'neg')
  } else {
    seurat_obj@meta.data[, paste0(listName, '_r', cor_threshold, '_max')] <-
      sapply(seq_along(colnames(seurat_obj)), function(x)
        colnames(seurat_obj@meta.data[, names(geneset_list)])[
          which.max(seurat_obj@meta.data[x, names(geneset_list)])])
  }
  if (!is.null(cor_threshold)) {
    message('cor_threshold provided! Filtering gene sets based on initial correlation to module score')
    filtered_geneset_list <- list()
    for (i in names(geneset_list)) {
      geneset_cor <- cor(seurat_obj@meta.data[, i],
        as.matrix(t(seurat_obj@assays$RNA@data[
          rownames(seurat_obj@assays$RNA@data) %in% geneset_list[[i]], ])))
      filtered_geneset_list[[i]] <-
        colnames(geneset_cor)[!is.na(geneset_cor) & geneset_cor > cor_threshold]
    }
    seurat_obj <- AddModuleScore(seurat_obj, filtered_geneset_list)
    seurat_obj@meta.data <- seurat_obj@meta.data[,
      !colnames(seurat_obj@meta.data) %in% paste0(names(geneset_list), '_r', cor_threshold)]
    colnames(seurat_obj@meta.data)[
      colnames(seurat_obj@meta.data) %in%
        paste0('Cluster', seq_along(geneset_list))] <-
      paste0(names(geneset_list), '_r', cor_threshold)
    if (length(geneset_list) == 1) {
      seurat_obj@meta.data[, paste0(listName, '_r', cor_threshold, '_max')] <-
        ifelse(seurat_obj@meta.data[,
          paste0(names(geneset_list), '_r', cor_threshold)] > pos_threshold,
          'pos', 'neg')
    } else {
      seurat_obj@meta.data[, paste0(listName, '_r', cor_threshold, '_max')] <-
        sapply(seq_along(colnames(seurat_obj)), function(x)
          colnames(seurat_obj@meta.data[,
            paste0(names(geneset_list), '_r', cor_threshold)])[
            which.max(seurat_obj@meta.data[x,
              paste0(names(geneset_list), '_r', cor_threshold)])])
    }
  }
  return(seurat_obj)
}

## Percent_Expressing replacement (scCustomize not available) -------------
Percent_Expressing <- function(seurat_object, features, group_by,
                               split_by = NULL, threshold = 0,
                               assay = 'RNA', layer = 'counts') {
  if ('layers' %in% slotNames(seurat_object@assays[[assay]])) {
    count_mat <- seurat_object@assays[[assay]]@layers$counts
    rownames(count_mat) <- rownames(seurat_object)
  } else {
    count_mat <- seurat_object@assays[[assay]]@counts
  }
  features_present <- features[features %in% rownames(count_mat)]
  count_mat <- count_mat[features_present, , drop = FALSE]
  meta <- seurat_object@meta.data
  if (is.null(split_by)) {
    groups <- as.character(meta[, group_by])
  } else {
    groups <- paste0(as.character(meta[, group_by]), '_',
                     as.character(meta[, split_by]))
  }
  unique_groups <- unique(groups)
  result <- sapply(unique_groups, function(g) {
    cells <- which(groups == g)
    apply(count_mat[, cells, drop = FALSE], 1,
          function(x) sum(x > threshold) / length(x))
  })
  if (is.vector(result))
    result <- matrix(result, nrow = 1,
                     dimnames = list(features_present, unique_groups))
  return(as.data.frame(result))
}

## geneDot (from useful_functions.R) --------------------------------------
geneDot <- function(seurat_obj = srt, gene = NULL, x = NULL, y = NULL,
                    x_name = 'genes', assay = 'RNA', y_name = 'clusters',
                    min_expression = 0, facet_ncol = 5, lim_expression = NULL,
                    scale.data = TRUE, plotcol = viridis::viridis(100),
                    include_NA = TRUE, swap_axes = FALSE, returnDF = FALSE) {
  require(tidyr)
  if (exists('levels_x')) rm('levels_x')
  if (exists('levels_y')) rm('levels_y')
  if (!is.null(y) && all(grepl('^\\d', seurat_obj@meta.data[, y])))
    seurat_obj@meta.data[, y] <- paste0('C', seurat_obj@meta.data[, y])
  if (is.factor(seurat_obj@meta.data[, x]))
    levels_x <- levels(seurat_obj@meta.data[, x]) else
    levels_x <- unique(seurat_obj@meta.data[, x])
  if (!is.null(y)) {
    if (is.factor(seurat_obj@meta.data[, y]))
      levels_y <- levels(seurat_obj@meta.data[, y]) else
      levels_y <- unique(seurat_obj@meta.data[, y])
  }
  seurat_obj@meta.data[, x] <- gsub('[-_ \\+]', '', seurat_obj@meta.data[, x])
  if (!is.null(y))
    seurat_obj@meta.data[, y] <- gsub('[-_ \\+]', '', seurat_obj@meta.data[, y])
  if (exists('levels_x')) levels_x <- gsub('[-_ \\+]', '', levels_x)
  if (exists('levels_y')) levels_y <- gsub('[-_ \\+]', '', levels_y)

  percent <- Percent_Expressing(seurat_object = seurat_obj, assay = assay,
    layer = 'counts', threshold = min_expression, features = gene,
    group_by = x, split_by = y)
  rownames(percent) <- gsub('_', '-', rownames(percent))
  if (length(unique(seurat_obj@meta.data[, x])) == 1)
    colnames(percent) <- gsub(paste0(unique(seurat_obj@meta.data[, x]), '_'), '',
                               colnames(percent))
  if (!is.null(y) && length(unique(seurat_obj@meta.data[, y])) == 1)
    colnames(percent) <- gsub(paste0('_', unique(seurat_obj@meta.data[, y])), '',
                               colnames(percent))
  percent$gene <- rownames(percent)
  percent <- gather(percent, groups, expression, 1:(ncol(percent) - 1))
  percent$key <- paste0(percent$gene, '_', percent$groups)
  colnames(percent)[colnames(percent) == 'expression'] <- 'percent'

  seurat_average <- log2(as.data.frame(
    AverageExpression(seurat_obj, assay = assay, features = gene,
                      group.by = c(x, if (!is.null(y)) y))[[1]]) + 1)
  if (length(gene) == 1) rownames(seurat_average) <- gene
  seurat_average$gene <- rownames(seurat_average)
  seurat_average <- seurat_average[gene, ]
  seurat_average <- gather(seurat_average, groups, expression,
                            1:(ncol(seurat_average) - 1))
  seurat_average$key <- paste0(seurat_average$gene, '_', seurat_average$groups)

  dot.df <- cbind(seurat_average,
                  percent[match(seurat_average$key, percent$key), ])
  dot.df <- dot.df[, !duplicated(colnames(dot.df))]
  if (scale.data)
    dot.df <- transform(dot.df, expression = ave(expression, gene,
      FUN = function(x) scale(x, scale = TRUE, center = TRUE)))
  dot.df$percent <- dot.df$percent * 100
  dot.df$expression[dot.df$percent == 0] <- NA
  dot.df$percent[dot.df$percent == 0]    <- NA
  dot.df$x_axis <- factor(sapply(dot.df$groups,
    function(x) unlist(strsplit(x, '_'))[1]), levels = unique(levels_x))
  dot.df$y_axis <- factor(dot.df$gene, levels = unique(gene))
  if (!is.null(y))
    dot.df$y_axis <- factor(sapply(dot.df$groups,
      function(x) unlist(strsplit(x, '_'))[2]), levels = levels_y)
  if (!is.null(y) && length(gene) > 1) {
    dot.df$y_axis <- factor(dot.df$gene, levels = unique(gene))
    dot.df$z_axis <- factor(sapply(dot.df$groups,
      function(x) unlist(strsplit(x, '_'))[2]), levels = unique(levels_y))
  }
  if (swap_axes)
    colnames(dot.df)[match(c('x_axis', 'y_axis'), colnames(dot.df))] <-
      c('y_axis', 'x_axis')
  p <- ggplot(data = dot.df, aes(x = x_axis, y = y_axis)) +
    geom_point(shape = 21, aes(fill = expression, size = percent),
               alpha = 0.7, colour = 'black', stroke = 0.3) +
    labs(x = x_name, y = y_name,
         title   = ifelse(length(gene) > 1, '', gene),
         subtitle = paste('Min expression >', min_expression)) +
    scale_shape(solid = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  if (!is.null(y) && length(gene) > 1) p <- p + facet_wrap(~z_axis)
  if (length(plotcol) > 3)
    p <- p + scale_fill_gradientn(colors = plotcol) else
    p <- p + scale_fill_gradient2(low = plotcol[1], mid = plotcol[2],
                                   high = plotcol[3])
  if (returnDF) return(dot.df) else return(p)
}

## cellComp (from useful_functions.R) -------------------------------------
cellComp <- function(seurat_obj = NULL, metaGroups = NULL,
                     plot_as = 'box', pal = NULL, prop = TRUE,
                     ptable_factor = 1, facet_ncol = 20,
                     facet_scales = 'free', subset_prop = NULL,
                     removeNA = TRUE, returnDF = FALSE,
                     metagroup_fill = NULL) {
  require(ggplot2); require(ggpubr)
  if (is.data.frame(seurat_obj))
    meta_groups_df <- seurat_obj[, metaGroups] else
    meta_groups_df <- seurat_obj@meta.data[, metaGroups]
  if (is.null(pal)) pal <- rainbow(length(unique(meta_groups_df[, 2])))
  if (prop)
    ccomp_df <- na.omit(as.data.frame(
      prop.table(table(meta_groups_df), ptable_factor))) else
    ccomp_df <- as.data.frame(table(meta_groups_df))
  if (removeNA) ccomp_df <- ccomp_df[ccomp_df$Freq != 0, ]
  if (!is.null(subset_prop)) {
    subset_col <- unlist(sapply(seq(ncol(ccomp_df)), function(x)
      if (any(ccomp_df[, x] %in% subset_prop)) colnames(ccomp_df)[x]))
    ccomp_df <- ccomp_df[ccomp_df[, subset_col] %in% subset_prop, ]
  }
  if (plot_as == 'box') {
    p <- ggplot(ccomp_df, aes_string(x = metaGroups[2], y = 'Freq')) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(values = pal) +
      xlab(metaGroups[2]) +
      ylab(ifelse(prop, 'proportion', 'counts'))
    if (length(metaGroups > 2))
      p <- p + geom_boxplot(aes_string(fill = metaGroups[3]),
                            color = 'black', outlier.size = .2,
                            alpha = 0.7, lwd = .2) else
      p <- p + geom_boxplot(aes_string(fill = metaGroups[2]),
                            color = 'black', outlier.size = .2,
                            alpha = 0.7, lwd = .2)
    if (length(metaGroups) > 3)
      p <- p + facet_wrap(as.formula(paste("~", metaGroups[4])),
                          scales = facet_scales, ncol = facet_ncol)
  }
  if (plot_as == 'bar') {
    metagroup_fill <- if (is.null(metagroup_fill)) metaGroups[2] else metagroup_fill
    p <- ggplot(ccomp_df, aes_string(x = metaGroups[1], y = 'Freq')) +
      geom_bar(position = "stack", stat = "identity",
               aes_string(fill = metagroup_fill), color = 'black') +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(values = pal) +
      xlab(metaGroups[2]) +
      ylab(ifelse(prop, 'proportion', 'counts'))
    if (length(metaGroups) == 3)
      p <- p + facet_wrap(as.formula(paste("~", metaGroups[3])),
                          scales = facet_scales, ncol = facet_ncol)
  }
  if (returnDF) return(ccomp_df) else return(p)
}

## NI.genes (from useful_functions.R) -------------------------------------
NI.genes <- function(all.genes = NULL, ni.goterms = NULL,
                     ni.genes = NULL, org = 'mouse') {
  require(GO.db)
  if (org == 'mouse') {
    require(org.Mm.eg.db)
    if (!is.null(ni.goterms)) {
      goterms <- toTable(GOTERM)
      goid    <- goterms[goterms$Term %in% ni.goterms, 'go_id']
      goids   <- toTable(revmap(org.Mm.egGO))
      ni.goterms.genes <- unique(goids[goids$go_id %in% goid, 'gene_id'])
      ni.goterms <- toTable(org.Mm.egSYMBOL[ni.goterms.genes])$symbol
    }
    if (!is.null(ni.genes)) {
      if (is.null(all.genes)) all.genes <- toTable(org.Mm.egSYMBOL)$symbol
      ni.genes <- unique(unlist(lapply(ni.genes,
        function(x) all.genes[grep(x, all.genes)])))
    }
    return(unique(c(ni.goterms, ni.genes)))
  }
  if (org == 'human') {
    require(org.Hs.eg.db)
    if (!is.null(ni.goterms)) {
      goterms <- toTable(GOTERM)
      goid    <- goterms[goterms$Term %in% ni.goterms, 'go_id']
      ni.goterms.genes <- unique(toTable(revmap(org.Hs.egGO)[goid])[, 'gene_id'])
      ni.goterms <- toTable(org.Hs.egSYMBOL[ni.goterms.genes])$symbol
    }
    if (!is.null(ni.genes)) {
      all.genes <- toTable(org.Hs.egSYMBOL)$symbol
      ni.genes  <- unique(unlist(lapply(ni.genes,
        function(x) all.genes[grep(x, all.genes)])))
    }
    return(unique(c(ni.goterms, ni.genes)))
  }
}

## dotGSEA (from useful_functions.R) --------------------------------------
dotGSEA <- function(enrichmentsTest_list, type = c('fgsea', 'enrich'),
                    top_pathways = NULL, padj_threshold = 0.05,
                    cluster_rows = TRUE, cluster_cols = TRUE,
                    remove_ns_modules = TRUE) {
  require(tidyr)
  if (type == 'fgsea') {
    sig_terms <- unique(unlist(sapply(enrichmentsTest_list,
      function(x) x[x$padj < padj_threshold, 'pathway'])))
    if (length(sig_terms) == 0) return(NULL)
    if (!is.null(top_pathways))
      sig_terms <- unique(unlist(sapply(enrichmentsTest_list, function(x) {
        x <- na.omit(x); x <- x[x$padj < padj_threshold, ]
        x_top <- head(x[order(x$padj), 'pathway'][order(-x[order(x$padj), 'NES'])],
                      top_pathways)
        x_bot <- head(x[order(x$padj), 'pathway'][order(x[order(x$padj),  'NES'])],
                      top_pathways)
        c(x_top, x_bot)
      })))
    if (length(sig_terms) == 1) {
      mat_pvalue  <- t(as.matrix(sapply(enrichmentsTest_list,
        function(x) x$pval[match(sig_terms, x$pathway)])))
      mat_sizeLog <- t(as.matrix(sapply(enrichmentsTest_list,
        function(x) x$NES[match(sig_terms, x$pathway)])))
    } else {
      mat_pvalue  <- sapply(enrichmentsTest_list,
        function(x) x$pval[match(sig_terms, x$pathway)])
      mat_sizeLog <- sapply(enrichmentsTest_list,
        function(x) x$NES[match(sig_terms, x$pathway)])
    }
  }
  if (type == 'enrich') {
    sig_terms <- na.omit(unique(unlist(sapply(enrichmentsTest_list,
      function(x) x[x$p.adjust < padj_threshold, 'ID']))))
    if (length(sig_terms) == 0) return(NULL)
    if (!is.null(top_pathways))
      sig_terms <- unique(as.vector(unlist(sapply(enrichmentsTest_list,
        function(x) na.omit(head(
          x[x$p.adjust < padj_threshold, 'ID'][
            order(x[x$p.adjust < padj_threshold, 'ID'])],
          top_pathways))))))
    if (length(sig_terms) == 1) {
      mat_pvalue  <- t(as.matrix(sapply(enrichmentsTest_list,
        function(x) x$p.adjust[match(sig_terms, x$ID)])))
      mat_sizeLog <- t(as.matrix(sapply(enrichmentsTest_list,
        function(x) x$Count[match(sig_terms, x$ID)])))
    } else {
      mat_pvalue  <- sapply(enrichmentsTest_list,
        function(x) x$p.adjust[match(sig_terms, x$ID)])
      mat_sizeLog <- sapply(enrichmentsTest_list,
        function(x) x$Count[match(sig_terms, x$ID)])
    }
  }
  mat_pvalue[is.na(mat_pvalue)]   <- 1
  mat_sizeLog[is.na(mat_sizeLog)] <- 0
  rownames(mat_pvalue)  <- sig_terms; colnames(mat_pvalue)  <- names(enrichmentsTest_list)
  rownames(mat_sizeLog) <- sig_terms; colnames(mat_sizeLog) <- names(enrichmentsTest_list)
  mat_sizeLog <- as.data.frame(mat_sizeLog)
  mat_pvalue  <- as.data.frame(mat_pvalue)
  if (cluster_rows && min(dim(mat_sizeLog)[1]) >= 2) {
    d_row <- dist(mat_sizeLog); d_row[is.na(d_row)] <- max(d_row)
    hc1_row     <- hclust(d_row, method = "ward.D")
    mat_sizeLog <- mat_sizeLog[hc1_row$order, ]
    mat_pvalue  <- mat_pvalue[hc1_row$order, ]
  }
  if (cluster_cols && min(dim(mat_sizeLog)[2]) >= 2) {
    d_col <- dist(t(mat_sizeLog)); d_col[is.na(d_col)] <- max(d_col)
    hc1_col     <- hclust(d_col, method = "ward.D")
    mat_sizeLog <- mat_sizeLog[, hc1_col$order]
    mat_pvalue  <- mat_pvalue[, hc1_col$order]
  }
  mat_sizeLog[mat_pvalue > padj_threshold] <- NA
  mat_pvalue[mat_pvalue  > padj_threshold] <- NA
  mat_pvalue <- -log10(as.data.frame(mat_pvalue))
  mat_pvalue$pathway <- rownames(mat_pvalue)
  if (remove_ns_modules) {
    keep_cols   <- apply(mat_pvalue, 2, function(x) !all(is.na(x)))
    mat_pvalue  <- mat_pvalue[, keep_cols, drop = FALSE]
    mat_sizeLog <- mat_sizeLog[, colnames(mat_pvalue)[colnames(mat_pvalue) != 'pathway'],
                               drop = FALSE]
  }
  if (type == 'fgsea') {
    gsea_df_pvalue <- gather(mat_pvalue, cluster, nlogPvalue,
      colnames(mat_pvalue)[1]:colnames(mat_pvalue)[ncol(mat_pvalue) - 1],
      factor_key = TRUE)
    gsea_df_NES <- gather(mat_sizeLog, cluster2, NES,
      colnames(mat_sizeLog)[1]:colnames(mat_sizeLog)[ncol(mat_sizeLog)],
      factor_key = TRUE)
    gsea_df <- cbind(gsea_df_NES, gsea_df_pvalue)
    gsea_df$pathway <- factor(gsea_df$pathway, levels = unique(gsea_df$pathway))
    p <- ggplot(data = gsea_df, aes(x = cluster, y = pathway,
                                    fill = NES, size = nlogPvalue)) +
      geom_point(shape = 21, color = 'black') +
      scale_fill_gradient2(low = "#67A9CF", mid = "#F7F7F7", high = "#EF8A62",
                           midpoint = 0) +
      labs(x = 'cluster', y = '-log10(pvalue)') +
      theme_minimal() +
      theme(text = element_text(size = 11),
            axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    return(p)
  }
  if (type == 'enrich') {
    enrich_df_pvalue <- gather(mat_pvalue, cluster, nlogPvalue,
      colnames(mat_pvalue)[1]:colnames(mat_pvalue)[ncol(mat_pvalue) - 1],
      factor_key = TRUE)
    enrich_df_Count <- gather(mat_sizeLog, cluster2, Count,
      colnames(mat_sizeLog)[1]:colnames(mat_sizeLog)[ncol(mat_sizeLog)],
      factor_key = TRUE)
    enrich_df <- cbind(enrich_df_Count, enrich_df_pvalue)
    enrich_df$pathway <- factor(enrich_df$pathway, levels = unique(enrich_df$pathway))
    p <- ggplot(data = enrich_df, aes(x = cluster, y = pathway,
                                      fill = nlogPvalue, size = log2(Count + 1))) +
      geom_point(shape = 21, color = 'black') +
      scale_fill_distiller(palette = "Spectral") +
      labs(x = 'cluster', y = '-log10(pvalue)') +
      theme_minimal() +
      theme(text = element_text(size = 11),
            axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    return(p)
  }
}

## diffClustHeat (from useful_functions.R) --------------------------------
diffClustHeat <- function(deg_genes_df = NULL, topGenes = 20,
                          pvalAdjTrheshold = 0.05,
                          pvalAdjTrheshold2 = 1e-2,
                          pvalAdjTrheshold3 = 1e-3,
                          sort_by = 'p_val_adj',
                          only_pos = FALSE,
                          addGene = NULL,
                          plotcol = NULL,
                          col_limit = NULL,
                          return_mat = FALSE, ...) {
  require(ComplexHeatmap); require(circlize)
  deg_genes_df1 <- deg_genes_df
  deg_genes_df1$p_val_adj[is.na(deg_genes_df1$p_val_adj)] <- 1
  sig_genes    <- unique(deg_genes_df1$gene[deg_genes_df1$p_val_adj < pvalAdjTrheshold])
  deg_genes_df <- deg_genes_df1[deg_genes_df1$gene %in% sig_genes, ]
  if (sort_by == 'p_val_adj') {
    geneNames <- unique(unlist(lapply(split(deg_genes_df, deg_genes_df$cluster),
      function(x) { x <- x[order(x[, sort_by]), ]; head(x, topGenes)$gene })))
  }
  if (sort_by == 'avg_log2FC') {
    geneNames <- unique(unlist(lapply(split(deg_genes_df, deg_genes_df$cluster),
      function(x) {
        x <- x[order(-x[, sort_by]), ]
        if (only_pos) head(x, topGenes)$gene else
          c(head(x, topGenes)$gene, tail(x, 2)$gene)
      })))
  }
  if (!is.null(addGene)) geneNames <- unique(append(geneNames, addGene))
  pval_mat <- split(deg_genes_df1, deg_genes_df1$cluster)
  pval_mat <- as.data.frame(sapply(pval_mat,
    function(x) x$p_val_adj[match(geneNames, x$gene)]))
  if (ncol(pval_mat) == 1) pval_mat <- as.data.frame(t(pval_mat[, , drop = FALSE]))
  rownames(pval_mat) <- geneNames
  pval_mat[is.na(pval_mat)] <- 1
  lfc_mat <- split(deg_genes_df1, deg_genes_df1$cluster)
  lfc_mat <- as.data.frame(sapply(lfc_mat,
    function(x) x$avg_log2FC[match(geneNames, x$gene)]))
  if (ncol(lfc_mat) == 1) lfc_mat <- as.data.frame(t(lfc_mat))
  rownames(lfc_mat) <- geneNames
  lfc_mat[is.na(lfc_mat)] <- 0
  if (is.null(col_limit)) col_limit <- max(abs(lfc_mat))
  if (is.null(plotcol)) {
    plotcol <- colorRamp2(c(-max(abs(lfc_mat)), 0, max(abs(lfc_mat))),
                          c("green", "white", "red"))
  } else {
    plotcol <- colorRamp2(c(-col_limit, 0, col_limit), plotcol)
  }
  labels_annotation <- HeatmapAnnotation(
    text = anno_text(rownames(lfc_mat), rot = 45,
                     location = unit(1, "npc"), just = "right",
                     gp = gpar(fontsize = 4)),
    annotation_height = max_text_width(rownames(lfc_mat))
  )
  ht <- Heatmap(t(lfc_mat),
    clustering_distance_columns = 'euclidean',
    clustering_distance_rows    = 'euclidean',
    col          = plotcol,
    top_annotation = labels_annotation,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 0),
    ...)
  return(ht)
}


### ===== 3. PATHS AND PARAMETERS =========================================
# -------------------------------------------------------------------------
# CONFIGURE THESE PATHS BEFORE RUNNING
# -------------------------------------------------------------------------

# Directory containing srt.rds downloaded from GEO (input only – read-only)
# Place the srt.rds file inside this directory.
# The reference object for label transfer (GSE154659_srt.rds) should also
# be placed here if you want to run label transfer from scratch;
# alternatively, a pre-computed prediction_scores_celltype.csv in the same
# directory will be used automatically.
srt_source_dir <- '/path/to/srt_directory'

# Directory where all outputs (plots, CSV files, RDS objects) will be written.
# Will be created if it does not exist.
projdir <- '/path/to/output_directory'

# Directory containing MSigDB GMT files (mouse).
# Expected sub-structure: <gsea_db_dir>/mouse/<annotation>.gmt
# Download from https://www.gsea-msigdb.org/gsea/downloads.jsp
gsea_db_dir <- '/path/to/GSEA_gmt_directory'

# -------------------------------------------------------------------------

setwd(projdir)
dir.create(file.path(projdir, 'Plots'), showWarnings = FALSE)

# Variables matching the scrna_pipeline processing parameters
reductionName     <- 'umap'
reductionSave     <- 'pca'
reductionGraphSnn <- 'RNA_snn'
reductionGraphKnn <- 'RNA_knn'
sigPCs            <- 15
org               <- 'mouse'


### ===== 4. LOAD SRT.RDS ==================================================

message('Loading srt.rds ...')
srt <- readRDS(file.path(srt_source_dir, 'srt.rds'))


### ===== 5. PALETTE SETUP =================================================

srt$Condition <- ifelse(grepl('SHAM', srt$condition), 'Sham', 'SNI')
if (!'genotype' %in% colnames(srt@meta.data)) {
  srt$genotype <- ifelse(grepl('KO', srt$condition), 'KO', 'WT')
}
srt$condgen  <- paste0(srt$Condition, '_', srt$genotype)
srt$condgen  <- gsub('_', ' ', srt$condgen)
srt$condgen  <- factor(srt$condgen,
                       levels = c('Sham WT', 'Sham KO', 'SNI WT', 'SNI KO'))
srt$condgen2 <- sub(' ', '_', as.character(srt$condgen))
srt$condgen2 <- factor(srt$condgen2,
                       levels = c('Sham_WT', 'Sham_KO', 'SNI_WT', 'SNI_KO'))

palette_condgen <- c('Sham WT' = 'grey', 'Sham KO' = 'firebrick1',
                     'SNI WT'  = 'darkslategrey', 'SNI KO' = 'violetred4')
condition_pal   <- setNames(c('darkseagreen4', 'orange'), c('Sham', 'SNI'))
genotype_pal    <- setNames(c('grey55', 'red'), c('WT', 'KO'))

srt$sampleID2 <- gsub('DRG-', '', srt$sampleID)
srt$sampleID2 <- gsub('Rgs4', '', srt$sampleID2)
srt$sampleID2 <- factor(srt$sampleID2,
  levels = c('Sham-WT-1','Sham-WT-2','Sham-WT-3','Sham-WT-4',
             'Sham-KO-1','Sham-KO-2','Sham-KO-3','Sham-KO-4',
             'SNI-WT-1','SNI-WT-2','SNI-WT-3','SNI-WT-4',
             'SNI-KO-1','SNI-KO-2','SNI-KO-3','SNI-KO-4'))

palette_condgen2 <- srt$condgen[!duplicated(srt$sampleID2)]
names(palette_condgen2) <- srt$sampleID2[!duplicated(srt$sampleID2)]
palette_condgen2 <- palette_condgen[palette_condgen2]
names(palette_condgen2) <- srt$sampleID2[!duplicated(srt$sampleID2)]


### ===== 6. LABEL TRANSFER ================================================

src_csv <- file.path(srt_source_dir, 'prediction_scores_celltype.csv')
out_csv <- file.path(projdir, 'prediction_scores_celltype.csv')

# Prefer already-computed predictions from source dir, then from projdir
if (file.exists(src_csv)) {
  message('Loading pre-computed label transfer predictions from source dir ...')
  predictionsMaj <- read.csv(src_csv)
  srt$celltype_predicted <- predictionsMaj$predicted.id
  if (!file.exists(out_csv)) file.copy(src_csv, out_csv)
} else if (file.exists(out_csv)) {
  message('Loading pre-computed label transfer predictions from projdir ...')
  predictionsMaj <- read.csv(out_csv)
  srt$celltype_predicted <- predictionsMaj$predicted.id
} else {
  ref_path <- file.path(srt_source_dir, 'GSE154659_srt.rds')
  if (file.exists(ref_path)) {
    ref <- readRDS(ref_path)
    message('Running label transfer ...')
    anchors       <- FindTransferAnchors(reference = ref, query = srt,
                                         dims = 1:15, npcs = 15)
    predictionsMaj <- TransferData(anchorset = anchors,
                                   refdata   = ref$celltype, dims = 1:15)
    write.csv(predictionsMaj, out_csv)
    srt$celltype_predicted <- predictionsMaj$predicted.id
  } else {
    message('No label transfer source available – skipping')
  }
}

# Manual cell-type annotation refinements
if ('celltype_predicted' %in% colnames(srt@meta.data)) {
  srt$celltype_predicted2 <- srt$celltype_predicted
  srt$celltype_predicted2[srt$RNA_snn_res.0.2 == 7] <- 'Macrophage'
  srt$celltype_predicted2[srt$RNA_snn_res.0.2 == 5] <- 'Endothelial'
  srt$celltype_predicted2[srt$RNA_snn_res.0.2 == 2] <- 'Fibroblast'
  srt$celltype_predicted2[srt$RNA_snn_res.0.2 == 8] <- 'Pericyte'
  srt$celltype_predicted2[srt$celltype_predicted2 %in% c('NF2','NF3')] <- 'NF2/NF3'
}


### ===== 7. CELL COMPOSITION BARPLOT =====================================

if ('celltype_predicted' %in% colnames(srt@meta.data)) {
  metadata <- srt@meta.data
  metadata <- metadata[
    metadata$celltype_predicted %in%
      names(table(metadata$celltype_predicted)[
        table(metadata$celltype_predicted) > 1]), ]
  bp <- cellComp(
    seurat_obj   = metadata,
    ptable_factor = 1,
    metaGroups   = c('celltype_predicted', 'condgen2'),
    plot_as      = 'bar',
    pal          = palette_condgen
  ) + gtheme

  pdf(file.path('Plots', 'cell_proportions_celltypes_conditions.pdf'),
      width = 6, height = 4)
  print(bp)
  dev.off()
}


### ===== 8. UMAP / tSNE DIMPLOTS =========================================

srt <- RunTSNE(object = srt, reduction = reductionSave, dims = 1:sigPCs)

p1 <- DimPlot(srt, group.by = 'celltype_predicted2', label = FALSE,
              label.size = 3) + theme_void()
p2 <- DimPlot(srt, group.by = 'Condition', label = FALSE,
              cols = condition_pal) + NoLegend() + theme_void()
p3 <- DimPlot(srt, group.by = 'genotype', label = FALSE,
              cols = genotype_pal) + NoLegend() + theme_void()
p4 <- DimPlot(srt, group.by = 'condgen', label = FALSE,
              cols = palette_condgen) + theme_void() +
  scale_fill_manual(values = palette_condgen)

png(file.path('Plots', 'label_transfer_celltype_umap.png'),
    width = 3200, height = 1000, res = 300)
wrap_plots(p1, p2, p4)
dev.off()

pdf(file.path('Plots', 'label_transfer_celltype_umap.pdf'), width = 12, height = 4)
wrap_plots(p1, p2, p4)
dev.off()

# Wt-only UMAP
srtWT <- RunUMAP(srt[, srt$condgen %in% c('Sham WT', 'SNI WT')],
                 reduction = reductionSave, dims = 1:sigPCs,
                 reduction.name = reductionName)
p_wt <- DimPlot(srtWT, group.by = 'celltype_predicted2',
                label = TRUE, label.size = 3) + theme_void()
pdf(file.path('Plots', 'label_transfer_celltype_WT_umap.pdf'), 4, 4)
print(p_wt)
dev.off()

png(file.path('Plots', 'label_transfer_celltype_WT_umap.png'),
    width = 1200, height = 1200, res = 300)
print(p_wt)
dev.off()

rm(srtWT)

# Clustering DimPlot
pdf(file.path('Plots', 'clustering_dimplot.pdf'), 4, 4)
DimPlot(srt, group.by = 'celltype_predicted2')
dev.off()


### ===== 9. QC VIOLIN PLOTS ==============================================

srt$nFeature_RNAL <- log10(srt$nFeature_RNA)
srt$nCount_RNAL   <- log10(srt$nCount_RNA)

vln_p <- VlnPlot(srt,
  features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"),
  combine = FALSE, group.by = 'sampleID2', pt.size = 0,
  ncol = 3, col = palette_condgen2)
vln_p <- lapply(vln_p, function(x)
  x + theme(axis.text = element_text(size = 9)) + NoLegend())

ccomp_df <- as.data.frame(table(srt$sampleID2))
cc_p2 <- ggplot(ccomp_df, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(position = "stack", stat = "identity") +
  ggtitle(paste('Tot cells', ncol(srt))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = palette_condgen2) +
  theme(axis.text = element_text(size = 9)) +
  NoLegend()

png("Plots/QC_nFeat_nCount_m.percent_vlnPlot.png",
    4200, 1000, res = 300)
print((cc_p2 | vln_p[[1]] | vln_p[[2]] | vln_p[[3]]) +
      plot_layout(widths = c(2, 2, 2, 2)))
dev.off()


### ===== 10. RGS4 EXPRESSION =============================================

pdf('Plots/RGS4_expression.pdf')
VlnPlot(srt, 'Rgs4', group.by = 'condgen', col = palette_condgen) +
  scale_fill_manual(values = palette_condgen)
dev.off()


### ===== 11. DEG STANDARD (per cluster markers) ==========================

enricher_universe <- 'all'
logfcThreshold    <- 0.25
pvalAdjTrheshold  <- 0.01
metaGroupName     <- paste0(reductionGraphSnn, '_res.0.8')
top_pathways      <- 5
force             <- FALSE
top_genes         <- 5

projdir_DEG <- paste0('DEG_', metaGroupName, '_logFC_', logfcThreshold,
                       '_padj_', pvalAdjTrheshold)
dir.create(file.path(projdir_DEG, 'Plots'), recursive = TRUE, showWarnings = FALSE)

if (is.numeric(srt@meta.data[, metaGroupName]))
  Idents(srt) <- factor(srt@meta.data[, metaGroupName],
    levels = unique(srt@meta.data[, metaGroupName])[
      order(as.numeric(as.character(unique(srt@meta.data[, metaGroupName]))))]) else
  Idents(srt) <- as.factor(srt@meta.data[, metaGroupName])

if (!file.exists(file.path(projdir_DEG, 'DEG.rds')) || force) {
  DefaultAssay(srt) <- 'RNA'
  degClusters <- FindAllMarkers(srt, max.cells.per.ident = 1000,
    min.pct = .1, logfc.threshold = logfcThreshold, verbose = TRUE)
  saveRDS(degClusters, file.path(projdir_DEG, 'DEG.rds'))
  write.csv(degClusters[degClusters$p_val_adj < pvalAdjTrheshold, ],
            file.path(projdir_DEG, 'DEG.csv'))
} else {
  message('DEG file found loading ...')
  degClusters <- readRDS(file.path(projdir_DEG, 'DEG.rds'))
}

degClusters_sig <- degClusters[degClusters$p_val_adj < pvalAdjTrheshold, ]

if (length(unique(degClusters_sig$cluster)) > 2) {
  top_deg <- degClusters_sig %>%
    filter(avg_log2FC > 0) %>%
    group_by(cluster) %>%
    top_n(top_genes, avg_log2FC)
} else {
  top_deg <- degClusters_sig[order(degClusters_sig$avg_log2FC), ]
  top_deg <- split(top_deg, top_deg$avg_log2FC > 0)
  if (length(top_deg) > 1)
    top_deg <- rbind(head(top_deg[[1]], top_genes), tail(top_deg[[2]], top_genes)) else
    top_deg <- head(top_deg[[1]], top_genes)
  top_deg <- as.data.frame(top_deg)
}

top_deg$cluster <- factor(top_deg$cluster,
  levels = names(table(srt@meta.data[, metaGroupName])[
    order(-table(srt@meta.data[, metaGroupName]))]))
top_deg <- top_deg[order(top_deg$cluster), ]

srt2    <- ScaleData(srt, features = unique(top_deg$gene))
Idents(srt2) <- factor(srt@meta.data[, metaGroupName],
  levels = names(table(srt@meta.data[, metaGroupName])[
    order(-table(srt@meta.data[, metaGroupName]))]))
heat_p  <- DoHeatmap(srt2, features = unique(top_deg$gene))

png(file.path(projdir_DEG, 'Plots', paste0('top_', top_genes, '_heatmap.png')),
    height = 6000, width = 16300, res = 400)
print(heat_p)
dev.off()
pdf(file.path(projdir_DEG, 'Plots', paste0('top_', top_genes, '_heatmap.pdf')),
    height = 6, width = 5)
print(heat_p)
dev.off()
rm(srt2)

# Feature plots
top_deg_fp <- degClusters_sig %>%
  arrange(-avg_log2FC) %>%
  group_by(cluster) %>%
  slice_head(n = top_genes)

feat_p <- lapply(seq_along(top_deg_fp$gene), function(x)
  FeaturePlot(srt, features = top_deg_fp$gene[x], keep.scale = 'all',
              order = FALSE, combine = FALSE, pt.size = .01,
              reduction = reductionName))
feat_p <- unlist(feat_p, recursive = FALSE)
for (i in seq_along(feat_p)) {
  feat_p[[i]] <- feat_p[[i]] + theme_void() + NoLegend() + NoAxes() +
    ggtitle(paste0('CL', top_deg_fp$cluster[i], ':', top_deg_fp$gene[i])) +
    scale_colour_gradientn(colours = viridis::turbo(100))
}
png(file.path(projdir_DEG, 'Plots', paste0('top_', top_genes, '_fplots.png')),
    (length(feat_p) * 50) + 1500, (length(feat_p) * 50) + 1500, res = 300)
print(wrap_plots(feat_p))
dev.off()

# Dotplot
top_deg_dp <- degClusters_sig %>%
  arrange(-avg_log2FC) %>%
  group_by(cluster) %>%
  slice_head(n = top_genes)
dp <- DotPlot(object = srt, features = rev(unique(top_deg_dp$gene)),
              scale = TRUE, group.by = metaGroupName) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major = element_line(colour = "gainsboro")) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.5)
dp$data$id <- factor(dp$data$id, levels = levels(top_deg_dp$cluster))

png(file.path(projdir_DEG, 'Plots', paste0('top_', top_genes, '_dotplot.png')),
    (length(unique(top_deg_dp$gene)) * 90) + 500,
    (length(unique(srt@meta.data[, metaGroupName])) * 102) + 250, res = 300)
print(dp)
dev.off()
pdf(file.path(projdir_DEG, 'Plots', paste0('top_', top_genes, '_dotplot.pdf')),
    width = 38, height = 14)
print(dp)
dev.off()

# Pathway enrichment (enricher)
gmt_annotations_enrich <- c('c2.cp.reactome.v7.1.symbol.gmt', 'c5.bp.v7.1.symbol.gmt')
if (!file.exists(file.path(projdir_DEG, 'Pathway_Enrichment_clusters.rds')) || force) {
  eu <- if (enricher_universe == 'all') rownames(srt) else
        if (enricher_universe == 'variable_genes') VariableFeatures(srt) else
        unique(degClusters_sig$gene)
  EnrichRResAll <- list()
  for (ann in gmt_annotations_enrich) {
    gmt.file <- file.path(gsea_db_dir, org, ann)
    pathways  <- read.gmt(gmt.file)
    message(paste('Compute enrichment per cluster using annotation:', ann))
    EnrichRResCluster <- list()
    for (i in unique(degClusters_sig$cluster)) {
      message(paste('EnrichR running cluster', i))
      degCluster <- degClusters_sig[degClusters_sig$cluster == i, ]
      sig_genes  <- degCluster[degCluster$avg_log2FC > 0, 'gene']
      egmt       <- enricher(sig_genes, TERM2GENE = pathways, universe = eu)
      EnrichRResCluster[[i]] <- egmt@result
    }
    EnrichRResAll[[ann]] <- EnrichRResCluster
  }
  saveRDS(EnrichRResAll, file.path(projdir_DEG, 'Pathway_Enrichment_clusters.rds'))
} else {
  EnrichRResAll <- readRDS(file.path(projdir_DEG, 'Pathway_Enrichment_clusters.rds'))
}

EnrichRes_dp <- lapply(EnrichRResAll, function(x)
  dotGSEA(enrichmentsTest_list = x, type = 'enrich',
          padj_threshold = pvalAdjTrheshold, top_pathways = top_pathways,
          cluster_rows = TRUE, cluster_cols = TRUE))

for (i in seq_along(EnrichRes_dp)) {
  ann_name <- names(EnrichRes_dp)[i]
  if (is.null(ann_name)) ann_name <- paste0('annotation_', i)
  if (!is.null(EnrichRes_dp[[i]])) {
    pdf(file.path(projdir_DEG, 'Plots',
        paste0('Pathway_Enrichment_', ann_name, '_dotplot.pdf')),
        width  = 3 + length(unique(EnrichRes_dp[[i]]$data$cluster)),
        height = 3 + length(unique(EnrichRes_dp[[i]]$data$pathway)) / 7)
    print(EnrichRes_dp[[i]])
    dev.off()
  }
}


### ===== 12. DEG2 – Sham WT vs Sham KO ===================================

force            <- FALSE
do.fgsea         <- TRUE
rankby           <- 'LFC'
logfcThreshold   <- 0
pvalAdjTrheshold <- 0.05
topGenes         <- 50
addGene          <- NULL
pval_column      <- 'pval'
metaGroupName1   <- 'celltype_predicted2'
metaGroupName2   <- 'condgen'
deg2Ident        <- c('Sham WT', 'Sham KO')

projdir_deg2 <- paste0('deg2_', metaGroupName2, '_', metaGroupName1, '_',
                        paste(deg2Ident, collapse = '_vs_'), '_lfc_', logfcThreshold)
dir.create(file.path(projdir_deg2, 'Plots'), recursive = TRUE, showWarnings = FALSE)

min.pct <- 0.1
if (!file.exists(file.path(projdir_deg2, 'DEG2_genes.rds')) || force) {
  deg2L <- list()
  Idents(srt) <- srt@meta.data[, metaGroupName1]
  for (i in unique(srt@meta.data[, metaGroupName1])) {
    message('Removing batchy genes')
    nigenes <- NI.genes(rownames(srt),
      ni.goterms = c('ribosome biogenesis','ncRNA processing',
                     'RNA processing','oxidation.reduction process','NADH oxidatio'),
      ni.genes   = c('^Mt', '^Rpl', '^mt-', '^Rps'))
    srt2 <- srt[!rownames(srt) %in% nigenes, ]
    if ('layers' %in% slotNames(srt@assays$RNA)) {
      count_mat <- srt@assays$RNA@layers$counts
      rownames(count_mat) <- rownames(srt)
    } else {
      count_mat <- srt2@assays$RNA@counts
    }
    is.expressed <- apply(count_mat, 1, function(x) sum(x > 0) / length(x))
    srt2 <- srt2[is.expressed >= min.pct, ]
    tryCatch({
      deg2Cluster <- FindMarkers(srt2, ident.1 = deg2Ident[1], ident.2 = deg2Ident[2],
        group.by = metaGroupName2, subset.ident = i, only.pos = FALSE,
        min.pct = min.pct, logfc.threshold = logfcThreshold, test.use = 'wilcox')
      if (nrow(deg2Cluster) < 1) stop('no genes')
      deg2Cluster$cluster <- i
      deg2Cluster$gene    <- rownames(deg2Cluster)
      deg2L[[i]] <- deg2Cluster
    }, error = function(e) message(paste('Skipping cluster', i, ':', e$message)))
  }
  deg2 <- do.call(rbind, deg2L)
  if (pval_column == 'pval') deg2$p_val_adj <- deg2$p_val
  write.csv(deg2[deg2$p_val_adj < pvalAdjTrheshold, ],
            file.path(projdir_deg2, 'DEG2_genes.csv'))
  saveRDS(deg2, file.path(projdir_deg2, 'DEG2_genes.rds'))
} else {
  message('DEG2-cluster object found!')
  deg2 <- readRDS(file.path(projdir_deg2, 'DEG2_genes.rds'))
}

deg2_sig <- deg2[deg2$p_val_adj < pvalAdjTrheshold, ]
clusters_idx <- unique(srt@meta.data[, metaGroupName1])

if (nrow(deg2_sig) > 0 && any(deg2_sig$avg_log2FC != 0)) {
  deg_sum_pos <- as.data.frame(table(deg2_sig$cluster[deg2_sig$avg_log2FC > 0]))
  rownames(deg_sum_pos) <- deg_sum_pos$Var1
  deg_sum_pos <- deg_sum_pos[clusters_idx, ]
  rownames(deg_sum_pos) <- clusters_idx
  deg_sum_neg <- as.data.frame(table(deg2_sig$cluster[deg2_sig$avg_log2FC < 0]))
  deg_sum_neg$Freq <- -1 * deg_sum_neg$Freq
  rownames(deg_sum_neg) <- deg_sum_neg$Var1
  deg_sum_neg <- deg_sum_neg[clusters_idx, ]
  rownames(deg_sum_neg) <- clusters_idx
  deg_sum <- cbind(deg_sum_pos, deg_sum_neg)
  colnames(deg_sum)[3:4] <- c('Var1_neg', 'Freq_neg')
  deg_bar <- ggplot(deg_sum, aes(x = rownames(deg_sum))) +
    geom_bar(aes(y = Freq), fill = 'red', color = 'black', stat = "identity") +
    geom_bar(aes(y = Freq_neg), fill = 'blue', color = 'black', stat = "identity",
             position = "identity") +
    ylim(c(-max(abs(na.omit(c(deg_sum$Freq_neg, deg_sum$Freq)))),
            max(na.omit(abs(c(deg_sum$Freq_neg, deg_sum$Freq)))))) +
    labs(title = "DEG per cluster", y = "Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  pdf(file.path(projdir_deg2, 'Plots', 'number_deg_cluster.pdf'), 4, 4)
  print(deg_bar)
  dev.off()
}

if (length(unique(deg2$cluster)) > 1 && nrow(deg2_sig) > 0) {
  dhm <- diffClustHeat(deg_genes_df = deg2, sort_by = 'p_val_adj',
    topGenes = topGenes, pvalAdjTrheshold = pvalAdjTrheshold, col_limit = NULL,
    plotcol = rev(RColorBrewer::brewer.pal(3, 'RdBu')),
    name   = paste(deg2Ident, collapse = '-'),
    cluster_columns = TRUE, cluster_rows = TRUE, border = TRUE)
  pdf(file.path(projdir_deg2, 'Plots',
      paste0('Top_', topGenes, '_genes_heatmap.pdf')),
      width = 3 + ncol(dhm@matrix) / 15, height = 3 + nrow(dhm@matrix) / 20)
  print(dhm)
  dev.off()

  deg2_top <- split(deg2, deg2$cluster)
  deg2_top <- lapply(deg2_top, function(x) {
    x <- x[order(-x$avg_log2FC), ]
    c(head(x$gene, topGenes), tail(x$gene, topGenes))
  })
  heat_p_per_cl <- lapply(seq_along(deg2_top), function(x) {
    srt_tmp <- srt[, srt@meta.data[, metaGroupName1] == names(deg2_top)[x]]
    srt_tmp <- ScaleData(srt_tmp, features = unique(unlist(deg2_top)))
    Idents(srt_tmp) <- srt_tmp@meta.data[, metaGroupName2]
    DoHeatmap(srt_tmp, features = deg2_top[[x]], raster = FALSE) +
      ggtitle(names(deg2_top)[x])
  })
  pdf(file.path(projdir_deg2, 'Plots',
      paste0('per_cluster_top_', topGenes, '_heatmap.pdf')), height = 12)
  print(heat_p_per_cl)
  dev.off()
}

# Volcano plots
logfcThreshold2 <- 0.5
vp <- list()
for (cl in unique(deg2$cluster)) {
  deg2_cl       <- deg2[deg2$cluster == cl, ]
  deg2_cl$label <- ''
  deg2_cl$label[abs(deg2_cl$avg_log2FC) > logfcThreshold2 &
                  deg2_cl$p_val_adj < pvalAdjTrheshold] <-
    deg2_cl$gene[abs(deg2_cl$avg_log2FC) > logfcThreshold2 &
                   deg2_cl$p_val_adj < pvalAdjTrheshold]
  deg2_cl$Direction <- ''
  deg2_cl$Direction[deg2_cl$avg_log2FC >  logfcThreshold2 &
                      deg2_cl$p_val_adj < pvalAdjTrheshold] <- deg2Ident[1]
  deg2_cl$Direction[deg2_cl$avg_log2FC < -logfcThreshold2 &
                      deg2_cl$p_val_adj < pvalAdjTrheshold] <- deg2Ident[2]
  vp[[cl]] <- ggplot(deg2_cl, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(size = 1, shape = 19, aes(color = Direction)) +
    geom_vline(xintercept =  logfcThreshold2, linetype = "dotted",
               color = "blue", size = 1) +
    geom_vline(xintercept = -logfcThreshold2, linetype = "dotted",
               color = "blue", size = 1) +
    geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype = "dotted",
               color = "blue", size = 1) +
    geom_text_repel(size = 2, data = deg2_cl, aes(label = label)) +
    ggtitle(cl) +
    scale_color_manual(values = c('grey77', 'red', 'green')) +
    theme_light()
}
png(file.path(projdir_deg2, 'Plots', 'volcano_plots.png'),
    width = 4900, height = 2000, res = 300)
print(wrap_plots(vp), ncol = 3)
dev.off()

# fGSEA
gmt_annotations_fgsea <- 'h.all.v7.1.symbol.gmt'
fgsea_file <- file.path(projdir_deg2,
  paste0('fGSEA_annotation_', gmt_annotations_fgsea, '_rankby_', rankby, '.rds'))
if (do.fgsea) {
  if (!file.exists(fgsea_file) || force) {
    fgseaResAll <- list()
    fgseaRanks  <- list()
    gmt.file  <- file.path(gsea_db_dir, org, gmt_annotations_fgsea)
    pathways  <- gmtPathways(gmt.file)
    pathways  <- lapply(pathways, function(x) x[!is.na(x) & x != 'NA' & x != ''])
    message(paste('fGSEA using annotation:', gmt_annotations_fgsea))
    fgseaResCluster <- list()
    for (i in unique(deg2$cluster)) {
      message(paste('fGSEA running cluster', i))
      deg2Cluster  <- deg2[deg2$cluster == i, ]
      fgsea_ranks  <- if (rankby == 'LFC') deg2Cluster$avg_log2FC else
        -log10(deg2Cluster$p_val + 1e-300) * sign(deg2Cluster$avg_log2FC)
      fgsea_ranks  <- setNames(fgsea_ranks, deg2Cluster$gene)
      fgsea_ranks  <- fgsea_ranks[fgsea_ranks != 0]
      tryCatch({
        fgseaRes   <- fgseaMultilevel(pathways, fgsea_ranks,
                                      minSize = 15, maxSize = 500, BPPARAM = NULL)
        fgseaResCol <- collapsePathways(fgseaRes, stats = fgsea_ranks,
                                        pathway = pathways)
        fgseaResCluster[[i]] <- fgseaRes[fgseaRes$pathway %in%
                                           fgseaResCol$mainPathways]
        fgseaRanks[[gmt_annotations_fgsea]][[i]] <- fgsea_ranks
      }, error = function(e)
        message(paste('fGSEA failed for cluster', i, ':', e$message)))
    }
    fgseaResAll[[gmt_annotations_fgsea]] <- fgseaResCluster
    saveRDS(fgseaResAll, fgsea_file)
    saveRDS(fgseaRanks,
      file.path(projdir_deg2,
        paste0('fgsea_ranks_', gmt_annotations_fgsea, '.rds')))
  } else {
    message('fGSEA object found!')
    fgseaResAll <- readRDS(fgsea_file)
    fgseaRanks  <- readRDS(file.path(projdir_deg2,
      paste0('fgsea_ranks_', gmt_annotations_fgsea, '.rds')))
  }

  if (length(fgseaResAll[[1]]) > 1) {
    fgseaResAll_dp <- lapply(fgseaResAll, function(y)
      dotGSEA(y, padj_threshold = pvalAdjTrheshold, type = 'fgsea',
              top_pathways = top_pathways, cluster_rows = TRUE, cluster_cols = TRUE))
    lapply(seq_along(fgseaResAll_dp), function(x) {
      if (!is.null(fgseaResAll_dp[[x]])) {
        pdf(file.path(projdir_deg2, 'Plots',
            paste0('fGSEA_annotation_', names(fgseaResAll_dp)[x], '_dotplots.pdf')),
            width = 8,
            height = 3 + length(unique(fgseaResAll_dp[[x]]$data$pathway)) / 7)
        print(fgseaResAll_dp[[x]])
        dev.off()
      }
    })
  }

  pdf(file.path(projdir_deg2, 'Plots', 'enrichment_plots.pdf'), 5, 3)
  for (x in seq_along(fgseaResAll)) {
    for (y in seq_along(fgseaResAll[[x]])) {
      sig_pathways <- fgseaResAll[[x]][[y]]$pathway[
        fgseaResAll[[x]][[y]]$padj < pvalAdjTrheshold]
      gmt.file  <- file.path(gsea_db_dir, org, names(fgseaResAll)[x])
      pathways  <- gmtPathways(gmt.file)
      for (z in sig_pathways) {
        print(plotEnrichment(pathways[[z]], fgseaRanks[[x]][[y]]) +
          labs(title = paste0(names(fgseaRanks[[x]])[y], '_', z)))
      }
    }
  }
  dev.off()
}


### ===== 13. DEG2 – SNI_RGS4_WT vs SNI_RGS4_KO ==========================

force            <- FALSE
do.fgsea         <- TRUE
rankby           <- 'LFC'
logfcThreshold   <- 0
pvalAdjTrheshold <- 0.05
topGenes         <- 50
addGene          <- NULL
pval_column      <- 'pval'
metaGroupName1   <- 'celltype_predicted2'
metaGroupName2   <- 'condition'
deg2Ident        <- c('SNI_RGS4_WT', 'SNI_RGS4_KO')

projdir_deg2 <- paste0('deg2_', metaGroupName2, '_', metaGroupName1, '_',
                        paste(deg2Ident, collapse = '_vs_'), '_lfc_', logfcThreshold)
dir.create(file.path(projdir_deg2, 'Plots'), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(projdir_deg2, 'DEG2_genes.rds')) || force) {
  deg2L <- list()
  Idents(srt) <- srt@meta.data[, metaGroupName1]
  for (i in unique(srt@meta.data[, metaGroupName1])) {
    nigenes <- NI.genes(rownames(srt),
      ni.goterms = c('ribosome biogenesis','ncRNA processing',
                     'RNA processing','oxidation.reduction process','NADH oxidatio'),
      ni.genes   = c('^Mt', '^Rpl', '^mt-', '^Rps'))
    srt2 <- srt[!rownames(srt) %in% nigenes, ]
    if ('layers' %in% slotNames(srt@assays$RNA)) {
      count_mat <- srt@assays$RNA@layers$counts
      rownames(count_mat) <- rownames(srt)
    } else {
      count_mat <- srt2@assays$RNA@counts
    }
    is.expressed <- apply(count_mat, 1, function(x) sum(x > 0) / length(x))
    srt2 <- srt2[is.expressed >= min.pct, ]
    tryCatch({
      deg2Cluster <- FindMarkers(srt2, ident.1 = deg2Ident[1], ident.2 = deg2Ident[2],
        group.by = metaGroupName2, subset.ident = i, only.pos = FALSE,
        min.pct = min.pct, logfc.threshold = logfcThreshold, test.use = 'wilcox')
      if (nrow(deg2Cluster) < 1) stop('no genes')
      deg2Cluster$cluster <- i
      deg2Cluster$gene    <- rownames(deg2Cluster)
      deg2L[[i]] <- deg2Cluster
    }, error = function(e) message(paste('Skipping cluster', i, ':', e$message)))
  }
  deg2 <- do.call(rbind, deg2L)
  if (pval_column == 'pval') deg2$p_val_adj <- deg2$p_val
  write.csv(deg2[deg2$p_val_adj < pvalAdjTrheshold, ],
            file.path(projdir_deg2, 'DEG2_genes.csv'))
  saveRDS(deg2, file.path(projdir_deg2, 'DEG2_genes.rds'))
} else {
  message('DEG2-cluster object found!')
  deg2 <- readRDS(file.path(projdir_deg2, 'DEG2_genes.rds'))
}

deg2_sig <- deg2[deg2$p_val_adj < pvalAdjTrheshold, ]

if (length(unique(deg2$cluster)) > 1 && nrow(deg2_sig) > 0) {
  dhm <- diffClustHeat(deg_genes_df = deg2, sort_by = 'p_val_adj',
    topGenes = topGenes, pvalAdjTrheshold = pvalAdjTrheshold, col_limit = NULL,
    plotcol = rev(RColorBrewer::brewer.pal(3, 'RdBu')),
    name   = paste(deg2Ident, collapse = '-'),
    cluster_columns = TRUE, cluster_rows = TRUE, border = TRUE)
  pdf(file.path(projdir_deg2, 'Plots',
      paste0('Top_', topGenes, '_genes_heatmap.pdf')),
      width = 3 + ncol(dhm@matrix) / 15, height = 3 + nrow(dhm@matrix) / 20)
  print(dhm)
  dev.off()
}

# Volcano plots
vp <- list()
for (cl in unique(deg2$cluster)) {
  deg2_cl       <- deg2[deg2$cluster == cl, ]
  deg2_cl$label <- ''
  deg2_cl$label[abs(deg2_cl$avg_log2FC) > 0.5 &
                  deg2_cl$p_val_adj < pvalAdjTrheshold] <-
    deg2_cl$gene[abs(deg2_cl$avg_log2FC) > 0.5 &
                   deg2_cl$p_val_adj < pvalAdjTrheshold]
  deg2_cl$Direction <- ''
  deg2_cl$Direction[deg2_cl$avg_log2FC >  0.5 &
                      deg2_cl$p_val_adj < pvalAdjTrheshold] <- deg2Ident[1]
  deg2_cl$Direction[deg2_cl$avg_log2FC < -0.5 &
                      deg2_cl$p_val_adj < pvalAdjTrheshold] <- deg2Ident[2]
  vp[[cl]] <- ggplot(deg2_cl, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(size = 1, shape = 19, aes(color = Direction)) +
    geom_vline(xintercept =  0.5, linetype = "dotted", color = "blue", size = 1) +
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "blue", size = 1) +
    geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype = "dotted",
               color = "blue", size = 1) +
    geom_text_repel(size = 2, data = deg2_cl, aes(label = label)) +
    ggtitle(cl) +
    scale_color_manual(values = c('grey77', 'red', 'green')) +
    theme_light()
}
png(file.path(projdir_deg2, 'Plots', 'volcano_plots.png'),
    width = 4900, height = 2000, res = 300)
print(wrap_plots(vp), ncol = 3)
dev.off()

fgsea_file2 <- file.path(projdir_deg2,
  paste0('fGSEA_annotation_', gmt_annotations_fgsea, '_rankby_', rankby, '.rds'))
if (do.fgsea) {
  if (!file.exists(fgsea_file2) || force) {
    fgseaResAll2 <- list()
    fgseaRanks2  <- list()
    gmt.file  <- file.path(gsea_db_dir, org, gmt_annotations_fgsea)
    pathways  <- gmtPathways(gmt.file)
    pathways  <- lapply(pathways, function(x) x[!is.na(x) & x != 'NA' & x != ''])
    fgseaResCluster2 <- list()
    for (i in unique(deg2$cluster)) {
      deg2Cluster  <- deg2[deg2$cluster == i, ]
      fgsea_ranks  <- if (rankby == 'LFC') deg2Cluster$avg_log2FC else
        -log10(deg2Cluster$p_val + 1e-300) * sign(deg2Cluster$avg_log2FC)
      fgsea_ranks  <- setNames(fgsea_ranks, deg2Cluster$gene)
      fgsea_ranks  <- fgsea_ranks[fgsea_ranks != 0]
      tryCatch({
        fgseaRes    <- fgseaMultilevel(pathways, fgsea_ranks,
                                       minSize = 15, maxSize = 500, BPPARAM = NULL)
        fgseaResCol <- collapsePathways(fgseaRes, stats = fgsea_ranks,
                                        pathway = pathways)
        fgseaResCluster2[[i]] <- fgseaRes[fgseaRes$pathway %in%
                                             fgseaResCol$mainPathways]
        fgseaRanks2[[gmt_annotations_fgsea]][[i]] <- fgsea_ranks
      }, error = function(e)
        message(paste('fGSEA failed for cluster', i, ':', e$message)))
    }
    fgseaResAll2[[gmt_annotations_fgsea]] <- fgseaResCluster2
    saveRDS(fgseaResAll2, fgsea_file2)
    saveRDS(fgseaRanks2,
      file.path(projdir_deg2,
        paste0('fgsea_ranks_', gmt_annotations_fgsea, '.rds')))
  } else {
    fgseaResAll2 <- readRDS(fgsea_file2)
    fgseaRanks2  <- readRDS(file.path(projdir_deg2,
      paste0('fgsea_ranks_', gmt_annotations_fgsea, '.rds')))
  }

  if (length(fgseaResAll2[[1]]) > 1) {
    fgseaResAll_dp2 <- lapply(fgseaResAll2, function(y)
      dotGSEA(y, padj_threshold = pvalAdjTrheshold, type = 'fgsea',
              top_pathways = top_pathways, cluster_rows = TRUE, cluster_cols = TRUE))
    lapply(seq_along(fgseaResAll_dp2), function(x) {
      if (!is.null(fgseaResAll_dp2[[x]])) {
        pdf(file.path(projdir_deg2, 'Plots',
            paste0('fGSEA_annotation_', names(fgseaResAll_dp2)[x], '_dotplots.pdf')),
            width  = 8,
            height = 3 + length(unique(fgseaResAll_dp2[[x]]$data$pathway)) / 7)
        print(fgseaResAll_dp2[[x]])
        dev.off()
      }
    })
  }
}


### ===== 14. GENE VIOLIN / BOXPLOTS =====================================

## --- Rgs4 across all cell types -----------------------------------------
gene <- 'Rgs4'
ccomp_df <- srt@meta.data
ccomp_df <- cbind(ccomp_df, as.data.frame(
  t(srt@assays$RNA@data[gene, , drop = FALSE])))

noise <- rnorm(n = nrow(ccomp_df)) / 100000
if (!all(ccomp_df[, gene] == ccomp_df[, gene][1]))
  ccomp_df[, gene] <- ccomp_df[, gene] + noise

box_p <- ggplot(ccomp_df, aes_string(x = 'condgen', y = gene)) +
  geom_violin(trim = TRUE, aes_string(fill = 'condgen'), scale = 'width') +
  gtheme +
  facet_wrap(~celltype_predicted2, ncol = 7) +
  scale_fill_manual(values = palette_condgen)

stat.test <- box_p$data %>%
  group_by(celltype_predicted2) %>%
  wilcox_test(reformulate('condgen', gene)) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = 'condgen', step.increase = 0.1)
box_p <- box_p +
  stat_pvalue_manual(stat.test, remove.bracket = FALSE,
                     bracket.nudge.y = 0.4, hide.ns = TRUE,
                     label = "p.adj.signif") +
  NoLegend()

write.csv(as.data.frame(stat.test[, c(2:9)]), paste0(gene, '_stats.csv'))
pdf(paste0('Plots/', gene, '_vlnplots.pdf'), height = 5, width = 10)
print(wrap_plots(box_p))
dev.off()

## geneDot Rgs4 by condgen ------------------------------------------------
palette_gene_expression <- viridis::inferno(100)
srt$condgen2_r <- factor(srt$condgen, levels = rev(levels(srt$condgen)))
dp_rgs4 <- geneDot(seurat_obj = srt, gene = gene,
  x = 'condgen2_r', y = NULL,
  min_expression = 0, facet_ncol = 5,
  scale.data = FALSE, x_name = '', y_name = '',
  plotcol = palette_gene_expression) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

pdf(paste0('Plots/', gene, '_expression_leg.pdf'), height = 4, width = 3)
dp_rgs4 + scale_fill_gradientn(colors = palette_gene_expression)
dev.off()

## Rgs4 per sample (boxplot) ----------------------------------------------
noise <- rnorm(n = nrow(ccomp_df)) / 100000
if (!all(ccomp_df[, gene] == ccomp_df[, gene][1]))
  ccomp_df[, gene] <- ccomp_df[, gene] + noise

box_p_sample <- ggplot(ccomp_df, aes_string(x = 'sampleID2', y = gene)) +
  geom_boxplot(aes_string(fill = 'condgen'), color = 'grey22',
               width = .5, alpha = 1, lwd = .2, outlier.shape = NA) +
  gtheme +
  scale_fill_manual(values = palette_condgen)

pdf(paste0('Plots/', gene, '_per_sample_vlnplots.pdf'), height = 4, width = 5)
print(wrap_plots(box_p_sample))
dev.off()

## geneDot Rgs4 per sample ------------------------------------------------
palette_gene_expression2 <- viridis::viridis(100)
srt$sampleID2_r <- factor(srt$sampleID2, levels = rev(levels(srt$sampleID2)))
dp_rgs4_sample <- geneDot(seurat_obj = srt, gene = gene,
  x = 'sampleID2_r', y = NULL,
  min_expression = 0, facet_ncol = 5,
  scale.data = FALSE, x_name = '', y_name = '',
  swap_axes = FALSE,
  plotcol = palette_gene_expression2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

pdf(paste0('Plots/', gene, '_per_sample_dotplot.pdf'), height = 2, width = 5)
dp_rgs4_sample
dev.off()

## Rgs4 in NP cells only --------------------------------------------------
ccomp_df_np <- srt@meta.data[srt$celltype_predicted2 == 'NP', ]
ccomp_df_np <- cbind(ccomp_df_np, as.data.frame(
  t(srt[, srt$celltype_predicted2 == 'NP']@assays$RNA@data[gene, , drop = FALSE])))

noise <- rnorm(n = nrow(ccomp_df_np)) / 100000
if (!all(ccomp_df_np[, gene] == ccomp_df_np[, gene][1]))
  ccomp_df_np[, gene] <- ccomp_df_np[, gene] + noise

box_np <- ggplot(ccomp_df_np, aes_string(x = 'sampleID2', y = gene)) +
  geom_boxplot(aes_string(fill = 'condgen'), color = 'grey22',
               width = .5, alpha = 1, lwd = .2, outlier.shape = NA) +
  gtheme +
  scale_fill_manual(values = palette_condgen)

pdf(paste0('Plots/', gene, '_per_sample_NP_boxplots.pdf'), height = 4, width = 5)
print(wrap_plots(box_np))
dev.off()

box_np_vln <- ggplot(ccomp_df_np, aes_string(x = 'sampleID2', y = gene)) +
  geom_violin(trim = TRUE, aes_string(fill = 'condgen'), scale = 'width') +
  gtheme +
  scale_fill_manual(values = palette_condgen)

pdf(paste0('Plots/', gene, '_per_sample_NP_vlnplots.pdf'), height = 4, width = 5)
print(wrap_plots(box_np_vln))
dev.off()

dp_np <- geneDot(seurat_obj = srt[, srt$celltype_predicted2 == 'NP'],
  gene = gene, x = 'sampleID2_r', y = NULL,
  min_expression = 0, facet_ncol = 5, scale.data = TRUE,
  x_name = '', y_name = '', swap_axes = FALSE,
  plotcol = palette_gene_expression2) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

pdf(paste0('Plots/', gene, '_per_sample_NP_dotplot.pdf'), height = 3, width = 5)
dp_np
dev.off()

## --- Phf24 across all cell types ----------------------------------------
gene <- 'Phf24'
ccomp_df <- srt@meta.data
ccomp_df <- cbind(ccomp_df, as.data.frame(
  t(srt@assays$RNA@data[gene, , drop = FALSE])))

noise <- rnorm(n = nrow(ccomp_df)) / 100000
if (!all(ccomp_df[, gene] == ccomp_df[, gene][1]))
  ccomp_df[, gene] <- ccomp_df[, gene] + noise

box_phf24 <- ggplot(ccomp_df, aes_string(x = 'condgen', y = gene)) +
  geom_violin(trim = TRUE, aes_string(fill = 'condgen'), scale = 'width') +
  gtheme +
  facet_wrap(~celltype_predicted2, ncol = 7) +
  scale_fill_manual(values = palette_condgen)

stat.test_phf24 <- box_phf24$data %>%
  group_by(celltype_predicted2) %>%
  wilcox_test(reformulate('condgen', gene)) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()
stat.test_phf24 <- stat.test_phf24 %>%
  add_xy_position(x = 'condgen', step.increase = 0.1)
box_phf24 <- box_phf24 +
  stat_pvalue_manual(stat.test_phf24, remove.bracket = FALSE,
                     bracket.nudge.y = 0.4, hide.ns = TRUE,
                     label = "p.adj.signif") +
  NoLegend()

write.csv(as.data.frame(stat.test_phf24[, c(2:9)]),
          paste0(gene, '_stats.csv'))
pdf(paste0('Plots/', gene, '_vlnplots.pdf'), height = 5, width = 10)
print(wrap_plots(box_phf24))
dev.off()

## per-sample Phf24 boxplot -----------------------------------------------
noise <- rnorm(n = nrow(ccomp_df)) / 100000
if (!all(ccomp_df[, gene] == ccomp_df[, gene][1]))
  ccomp_df[, gene] <- ccomp_df[, gene] + noise

box_phf24_sample <- ggplot(ccomp_df, aes_string(x = 'sampleID2', y = gene)) +
  geom_boxplot(aes_string(fill = 'condgen'), color = 'grey22',
               width = .5, alpha = 1, lwd = .2, outlier.shape = NA) +
  gtheme +
  scale_fill_manual(values = palette_condgen)

pdf(paste0('Plots/', gene, '_per_sample_vlnplots.pdf'), height = 4, width = 5)
print(wrap_plots(box_phf24_sample))
dev.off()


### ===== 15. Mrgprd / Calca VLN + DOT PLOTS ==============================

gene_pair <- c('Mrgprd', 'Calca')

srt$condition2 <- ifelse(
  srt$condition %in% c('SNI_RGS4_KO', 'SNI_RGS4_WT'), 'SNI', 'SHAM')

pdf(file.path('Plots', paste0(paste(gene_pair, collapse = '_'), '_dotplot.pdf')),
    height = 5)
DotPlot(srt[, srt$condition %in% c('SNI_RGS4_WT', 'SHAM_RGS4_WT')],
        features = gene_pair, group.by = 'condition')
DotPlot(srt[, srt$condition %in% c('SNI_RGS4_WT', 'SNI_RGS4_KO')],
        features = gene_pair, group.by = 'condition')
dev.off()

# Contrast violin plots
for (contrast in list(c('SNI_RGS4_WT', 'SNI_RGS4_KO'),
                      c('SNI_RGS4_WT', 'SHAM_RGS4_WT'))) {
  ccomp_df_c <- srt[, srt$condition %in% contrast]@meta.data
  genes_mat  <- as.data.frame(
    t(srt[, srt$condition %in% contrast]@assays$RNA@data[gene_pair, , drop = FALSE]))
  ccomp_df_c <- cbind(ccomp_df_c, genes_mat)

  boxL <- list()
  for (g in gene_pair) {
    ccomp_df1 <- ccomp_df_c[ccomp_df_c[, g] > 0, ]
    box_c <- ggplot(ccomp_df1, aes_string(x = 'condition', y = g)) +
      geom_violin(trim = TRUE, aes_string(fill = 'condition')) +
      geom_boxplot(aes_string(fill = 'condition')) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    st <- box_c$data %>%
      wilcox_test(reformulate('condition', g)) %>%
      adjust_pvalue(method = "none") %>%
      add_significance()
    st <- st %>% add_xy_position(x = 'condition', step.increase = 0.1)
    boxL[[g]] <- box_c +
      stat_pvalue_manual(st, remove.bracket = FALSE, bracket.nudge.y = 0,
                         hide.ns = FALSE, label = "p.adj.signif") +
      NoLegend()
  }
  pdf(paste0('Plots/', paste(gene_pair, collapse = '_'), '_',
             paste(contrast, collapse = '_'), '_vlnPlot.pdf'), height = 5)
  print(wrap_plots(boxL))
  dev.off()

  # Per-sample means
  ccomp_df_agg <- aggregate(ccomp_df_c[, gene_pair], drop = FALSE,
    by = as.list(ccomp_df_c[, c('sampleID', 'condition'), drop = FALSE]), mean)
  ccomp_df_agg <- na.omit(ccomp_df_agg)
  boxL2 <- list()
  for (g in gene_pair) {
    ccomp_df_a1 <- ccomp_df_agg[ccomp_df_agg[, g] > 0, ]
    box_a <- ggplot(ccomp_df_a1, aes_string(x = 'condition', y = g)) +
      geom_violin(trim = TRUE, aes_string(fill = 'condition')) +
      geom_boxplot(aes_string(fill = 'condition')) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    st_a <- box_a$data %>%
      t_test(reformulate('condition', g)) %>%
      adjust_pvalue(method = "none") %>%
      add_significance()
    st_a <- st_a %>% add_xy_position(x = 'condition', step.increase = 0.1)
    boxL2[[g]] <- box_a +
      stat_pvalue_manual(st_a, remove.bracket = FALSE, bracket.nudge.y = 0,
                         hide.ns = FALSE, label = "p.adj.signif") +
      NoLegend()
  }
  pdf(paste0('Plots/', paste(gene_pair, collapse = '_'), '_',
             paste(contrast, collapse = '_'), '_vlnPlot_per_sample.pdf'),
      height = 5)
  print(wrap_plots(boxL2))
  dev.off()
}


### ===== 16. Hdac1 DOTPLOT ===============================================

gene_hdac <- 'Hdac1'
srt$condition2b <- gsub('RGS4_', '', srt$condition)

pdf(file.path('Plots', paste0(gene_hdac, '_dotplot.pdf')), height = 5)
DotPlot(srt[, srt$condition %in% c('SNI_RGS4_WT', 'SHAM_RGS4_WT')],
        features = gene_hdac, group.by = 'celltype_predicted',
        split.by  = 'condition2b')
dev.off()


### ===== 17. Phf24 VLN WITH STAT TESTS ===================================

gene <- 'Phf24'
pdf(paste0('Plots/', gene, '_vlnplot.pdf'))
VlnPlot(srt, gene, group.by = 'condition')
dev.off()

message('Done! All visualisations written to: ', file.path(projdir, 'Plots'))
