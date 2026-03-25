########################################################################
# Self-contained analysis script – Serafini RGS4 mouse DRG study: NP sub-cluster
# Starts from NP-subset srt.rds produced by the scrna_pipeline
# Does NOT run SCENIC or WGCNA; reads pre-computed outputs when available
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

## ggplot theme ------------------------------------------------------------
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

## fp – quick FeaturePlot --------------------------------------------------
fp <- function(seurat_obj, gene, reduction = reductionName) {
  require(viridis); require(Seurat)
  p <- FeaturePlot(seurat_obj, features = gene, keep.scale = 'all',
                   order = FALSE, combine = FALSE, pt.size = .01,
                   reduction = reduction)
  lapply(p, function(x)
    x + theme_void() + NoAxes() +
      ggtitle(colnames(x$data)[4]) +
      scale_colour_gradientn(colours = viridis::turbo(100)))
}

## ModScoreCor -------------------------------------------------------------
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
  message(paste('Annotate cells based on highest module score, column:',
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
      !colnames(seurat_obj@meta.data) %in%
        paste0(names(geneset_list), '_r', cor_threshold)]
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

## Percent_Expressing (replaces scCustomize) --------------------------------
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

## geneDot -----------------------------------------------------------------
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

## cellComp ----------------------------------------------------------------
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

## NI.genes ----------------------------------------------------------------
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

## dotGSEA -----------------------------------------------------------------
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

## diffClustHeat -----------------------------------------------------------
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

# Directory containing the NP sub-cluster srt.rds downloaded from GEO
# (input only – read-only). Place srt.rds inside this directory.
srt_source_dir <- '/path/to/NP_srt_directory'

# Pre-computed SCENIC AUC matrix (optional).
# If this file is absent the SCENIC Atf6 section is skipped gracefully.
# Expected file: <srt_source_dir>/SCENIC/vg_5000_mw_tss500bp/auc_mtx.csv
scenic_auc_csv <- file.path(srt_source_dir, 'SCENIC', 'vg_5000_mw_tss500bp',
                             'auc_mtx.csv')

# Pre-computed WGCNA module assignments (optional).
# If this file is absent the WGCNA green module section is skipped gracefully.
# Expected file: <srt_source_dir>/WGCNA_.../module_assignments.csv
wgcna_modules_csv <- file.path(srt_source_dir,
  'WGCNA_sp_3_ds_1_mch_0.2max_shared_15_minModuleSize_10_metasize_20_genes_3000',
  'module_assignments.csv')

# Directory where all outputs will be written. Will be created if absent.
projdir <- '/path/to/output_directory'

# Directory containing MSigDB GMT files (mouse).
# Expected sub-structure: <gsea_db_dir>/mouse/<annotation>.gmt
# Download from https://www.gsea-msigdb.org/gsea/downloads.jsp
gsea_db_dir <- '/path/to/GSEA_gmt_directory'

# -------------------------------------------------------------------------

setwd(projdir)
dir.create(file.path(projdir, 'Plots'), showWarnings = FALSE, recursive = TRUE)

# Variables matching the scrna_pipeline processing parameters
reductionName     <- 'umap'
reductionSave     <- 'pca'
reductionGraphSnn <- 'RNA_snn'
sigPCs            <- 15
org               <- 'mouse'


### ===== 4. LOAD NP SRT.RDS ===============================================

message('Loading NP srt.rds ...')
srt <- readRDS(file.path(srt_source_dir, 'srt.rds'))


### ===== 5. METADATA SETUP ================================================

srt$Condition <- ifelse(grepl('SHAM', srt$condition), 'Sham', 'SNI')
if (!'genotype' %in% colnames(srt@meta.data))
  srt$genotype <- ifelse(grepl('KO', srt$condition), 'KO', 'WT')

srt$condgen  <- paste0(srt$Condition, '_', srt$genotype)
srt$condgen  <- gsub('_', ' ', srt$condgen)
srt$condgen  <- factor(srt$condgen,
                       levels = c('Sham WT', 'Sham KO', 'SNI WT', 'SNI KO'))
srt$condgen2 <- sub(' ', '_', as.character(srt$condgen))
srt$condgen2 <- factor(srt$condgen2,
                       levels = c('Sham_WT', 'Sham_KO', 'SNI_WT', 'SNI_KO'))

palette_condgen  <- c('Sham WT' = 'grey', 'Sham KO' = 'firebrick1',
                      'SNI WT'  = 'darkslategrey', 'SNI KO' = 'violetred4')
condition_pal    <- setNames(c('darkseagreen4', 'orange'), c('Sham', 'SNI'))
genotype_pal     <- setNames(c('grey55', 'red'), c('WT', 'KO'))

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


### ===== 6. UMAP OVERVIEW ================================================

p1 <- DimPlot(srt, reduction = reductionName,
              group.by = paste0(reductionGraphSnn, '_res.0.8'),
              label = TRUE) + NoLegend() + theme_void()
p2 <- DimPlot(srt, reduction = reductionName, group.by = 'condgen',
              cols = palette_condgen) + theme_void()
p3 <- DimPlot(srt, reduction = reductionName, group.by = 'sampleID',
              label = FALSE) + theme_void()

pdf('Plots/NP_umap_overview.pdf', height = 4, width = 12)
print(p1 | p2 | p3)
dev.off()

pdf('Plots/NP_umap_condgen.pdf', height = 4, width = 5)
print(p2)
dev.off()


### ===== 7. DEG STANDARD (per cluster markers) ===========================

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


### ===== 8. CELL COMPOSITION (condition) ==================================

metaGroupName1_cc <- paste0(reductionGraphSnn, '_res.0.8')
metaGroupName2_cc <- 'condition'

cc_box <- cellComp(
  seurat_obj  = srt,
  metaGroups  = c('sampleID', metaGroupName2_cc, metaGroupName1_cc),
  plot_as     = 'box',
  facet_ncol  = 8
) + gtheme

cc_df <- cc_box$data
stat.test <- cc_df %>%
  group_by_at(metaGroupName1_cc) %>%
  t_test(reformulate(metaGroupName2_cc, 'Freq')) %>%
  adjust_pvalue(method = "none") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = metaGroupName1_cc,
                                            step.increase = 0.1)

png(file.path('Plots', paste0('cell_composition_', metaGroupName1_cc,
    '_barboxplots2.png')), width = 2000, height = 1000, res = 300)
print(cc_box + stat_pvalue_manual(stat.test, remove.bracket = FALSE,
  bracket.nudge.y = 0, hide.ns = TRUE, label = "p.adj.signif"))
dev.off()

pdf(file.path('Plots', paste0('cell_composition_', metaGroupName1_cc,
    '_barboxplots2.pdf')), width = 10, height = 6)
print(cc_box + stat_pvalue_manual(stat.test, remove.bracket = FALSE,
  bracket.nudge.y = 0, hide.ns = TRUE, label = "p.adj.signif"))
dev.off()


### ===== 9. SCENIC ATF6 ANALYSIS =========================================

if (file.exists(scenic_auc_csv)) {
  message('Loading SCENIC auc_mtx.csv ...')
  auc_mtx <- read.csv(scenic_auc_csv, header = TRUE)
  rownames(auc_mtx) <- auc_mtx[, 1]
  auc_mtx <- auc_mtx[, -1]
  # Only add rows (cells) that are in srt
  common_cells <- intersect(rownames(auc_mtx), colnames(srt))
  if (length(common_cells) > 0) {
    srt <- AddMetaData(srt, metadata = auc_mtx[common_cells, , drop = FALSE])
  }

  # Identify Atf6 regulon column (R turns Atf6(+) -> Atf6...)
  atf6_col <- colnames(srt@meta.data)[grep('^Atf6', colnames(srt@meta.data))]
  if (length(atf6_col) > 0) {
    atf6_col <- atf6_col[1]
    message(paste('Using SCENIC column:', atf6_col))

    metaGroupNames_sc <- c('sampleID', 'condgen')

    ## UMAP feature plot
    umap_atf6 <- FeaturePlot(srt, features = atf6_col, ncol = 1,
                              combine = FALSE, pt.size = 1,
                              reduction = reductionName)
    umap_atf6[[1]] <- umap_atf6[[1]] + theme_void() + NoAxes() +
      ggtitle(atf6_col) +
      theme(plot.title = element_text(size = 8)) +
      scale_colour_gradientn(colours = viridis::viridis(100),
                             limits = c(0, max(umap_atf6[[1]]$data[, 4])))

    pdf('Plots/SCENIC_Atf6_distribution_umap.pdf', height = 2.7, width = 3)
    print(umap_atf6[[1]])
    dev.off()

    ## Boxplots per condgen (pseudo-bulk per sample)
    ccomp_df_sc <- srt@meta.data[, c(atf6_col, metaGroupNames_sc), drop = FALSE]
    ccomp_df_sc <- aggregate(ccomp_df_sc[, atf6_col, drop = FALSE],
      by = as.list(ccomp_df_sc[, metaGroupNames_sc, drop = FALSE]), mean)

    sample_rep <- table(srt@meta.data[, 'condgen'])
    keep_group <- names(sample_rep[sample_rep >= 2])

    box_sc <- ggplot(ccomp_df_sc, aes_string(x = 'condgen', y = atf6_col)) +
      geom_boxplot(aes_string(fill = 'condgen'), outlier.size = 2,
                   outlier.alpha = 0.2, notch = FALSE, lwd = .2) +
      scale_fill_manual(values = palette_condgen) +
      ggtitle(atf6_col) +
      ylab('Atf6 regulon score') +
      gtheme

    st_sc <- box_sc$data[box_sc$data$condgen %in% keep_group, ] %>%
      t_test(reformulate('condgen', atf6_col)) %>%
      adjust_pvalue(method = "none") %>%
      add_significance()
    st_sc <- st_sc %>% add_xy_position(x = 'condgen', step.increase = 0.1)
    box_sc_p <- box_sc + stat_pvalue_manual(st_sc, remove.bracket = FALSE,
      bracket.nudge.y = 0, hide.ns = TRUE, label = "p.adj.signif") + NoLegend()

    pdf('Plots/SCENIC_Atf6_distribution_boxplots.pdf', height = 2.7, width = 1.7)
    print(box_sc_p)
    dev.off()

    ## WT only boxplot (Sham WT vs SNI WT)
    ccomp_df_wt <- ccomp_df_sc[ccomp_df_sc$condgen %in% c('Sham WT', 'SNI WT'), ]
    ccomp_df_wt$condgen <- factor(as.character(ccomp_df_wt$condgen),
                                   levels = c('Sham WT', 'SNI WT'))
    box_wt <- ggplot(ccomp_df_wt, aes_string(x = 'condgen', y = atf6_col)) +
      geom_boxplot(aes_string(fill = 'condgen'), outlier.size = 2,
                   outlier.alpha = 0.2, notch = FALSE, lwd = .2) +
      scale_fill_manual(values = palette_condgen) +
      ggtitle(atf6_col) + ylab('Atf6 regulon score') + gtheme
    st_wt <- box_wt$data %>%
      t_test(reformulate('condgen', atf6_col)) %>%
      adjust_pvalue(method = "none") %>%
      add_significance()
    st_wt <- st_wt %>% add_xy_position(x = 'condgen', step.increase = 0.1)
    box_wt_p <- box_wt + stat_pvalue_manual(st_wt, remove.bracket = FALSE,
      bracket.nudge.y = 0, hide.ns = TRUE, label = "p.adj.signif") + NoLegend()

    pdf('Plots/SCENIC_Atf6_WT_only_boxplots.pdf', height = 2.7, width = 1.7)
    print(box_wt_p)
    dev.off()

    ## DEG2 – Atf6+ vs rest: SNI_RGS4_KO vs SNI_RGS4_WT ------------------
    atf6_threshold <- 0.05
    atf6_pos_col   <- paste0(atf6_col, '_pos')
    srt@meta.data[, atf6_pos_col] <- ifelse(
      srt@meta.data[, atf6_col] > atf6_threshold, atf6_col, 'rest')

    force            <- FALSE
    do.fgsea         <- TRUE
    rankby           <- 'LFC'
    logfcThreshold   <- 0
    pvalAdjTrheshold <- 0.05
    topGenes         <- 50
    pval_column      <- 'pval'
    metaGroupName1   <- atf6_pos_col
    metaGroupName2   <- 'condition'
    deg2Ident        <- c('SNI_RGS4_KO', 'SNI_RGS4_WT')
    top_pathways     <- 5
    gmt_annotations_fgsea <- 'h.all.v7.1.symbol.gmt'

    projdir_deg2 <- paste0('deg2_', metaGroupName2, '_', atf6_col,
                            '_pos_', paste(deg2Ident, collapse = '_vs_'),
                            '_lfc_', logfcThreshold)
    dir.create(file.path(projdir_deg2, 'Plots'), recursive = TRUE,
               showWarnings = FALSE)

    min.pct <- 0.1
    if (!file.exists(file.path(projdir_deg2, 'DEG2_genes.rds')) || force) {
      deg2L <- list()
      Idents(srt) <- srt@meta.data[, metaGroupName1]
      for (i in unique(srt@meta.data[, metaGroupName1])) {
        nigenes <- NI.genes(rownames(srt),
          ni.goterms = c('ribosome biogenesis','ncRNA processing',
                         'RNA processing','oxidation.reduction process'),
          ni.genes   = c('^Mt', '^Rpl', '^mt-', '^Rps'))
        srt2 <- srt[!rownames(srt) %in% nigenes, ]
        if ('layers' %in% slotNames(srt@assays$RNA)) {
          count_mat <- srt@assays$RNA@layers$counts
          rownames(count_mat) <- rownames(srt)
        } else {
          count_mat <- srt2@assays$RNA@counts
        }
        is.expressed <- apply(count_mat, 1,
                              function(x) sum(x > 0) / length(x))
        srt2 <- srt2[is.expressed >= min.pct, ]
        tryCatch({
          deg2Cluster <- FindMarkers(srt2,
            ident.1 = deg2Ident[1], ident.2 = deg2Ident[2],
            group.by = metaGroupName2, subset.ident = i,
            only.pos = FALSE, min.pct = min.pct,
            logfc.threshold = logfcThreshold, test.use = 'wilcox')
          if (nrow(deg2Cluster) < 1) stop('no genes')
          deg2Cluster$cluster <- i
          deg2Cluster$gene    <- rownames(deg2Cluster)
          deg2L[[i]] <- deg2Cluster
        }, error = function(e)
          message(paste('Skipping cluster', i, ':', e$message)))
      }
      deg2 <- do.call(rbind, deg2L)
      if (!is.null(deg2) && pval_column == 'pval')
        deg2$p_val_adj <- deg2$p_val
      write.csv(
        deg2[deg2$p_val_adj < pvalAdjTrheshold, ],
        file.path(projdir_deg2, 'DEG2_genes.csv'))
      saveRDS(deg2, file.path(projdir_deg2, 'DEG2_genes.rds'))
    } else {
      message('DEG2-Atf6 object found!')
      deg2 <- readRDS(file.path(projdir_deg2, 'DEG2_genes.rds'))
    }

    if (!is.null(deg2) && nrow(deg2) > 0) {
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
          geom_bar(aes(y = Freq), fill = 'red', color = 'black',
                   stat = "identity", alpha = 0.8) +
          geom_bar(aes(y = Freq_neg), fill = 'blue', color = 'black',
                   stat = "identity", position = "identity", alpha = 0.8) +
          ylim(c(-max(abs(na.omit(c(deg_sum$Freq_neg, deg_sum$Freq)))),
                  max(na.omit(abs(c(deg_sum$Freq_neg, deg_sum$Freq)))))) +
          labs(title = "DEG per cluster", y = "Value") +
          gtheme
        pdf(file.path(projdir_deg2, 'Plots', 'NP_Atf6_number_deg_cluster.pdf'),
            height = 4, width = 1)
        print(deg_bar)
        dev.off()
      }

      if (length(unique(deg2$cluster)) > 1 && nrow(deg2_sig) > 0) {
        dhm <- diffClustHeat(deg_genes_df = deg2, sort_by = 'p_val_adj',
          topGenes = topGenes, pvalAdjTrheshold = pvalAdjTrheshold,
          col_limit = NULL,
          plotcol = rev(RColorBrewer::brewer.pal(3, 'RdBu')),
          name   = paste(deg2Ident, collapse = '-'),
          cluster_columns = TRUE, cluster_rows = TRUE, border = TRUE)
        pdf(file.path(projdir_deg2, 'Plots',
            paste0('Top_', topGenes, '_genes_heatmap.pdf')),
            width = 3 + ncol(dhm@matrix) / 15,
            height = 3 + nrow(dhm@matrix) / 20)
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
        vp[[cl]] <- ggplot(deg2_cl, aes(x = avg_log2FC,
                                         y = -log10(p_val_adj))) +
          geom_point(size = 1, shape = 19, aes(color = Direction)) +
          geom_vline(xintercept =  0.5, linetype = "dotted",
                     color = "blue", size = 1) +
          geom_vline(xintercept = -0.5, linetype = "dotted",
                     color = "blue", size = 1) +
          geom_hline(yintercept = -log10(pvalAdjTrheshold),
                     linetype = "dotted", color = "blue", size = 1) +
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
      fgsea_file <- file.path(projdir_deg2,
        paste0('fGSEA_annotation_', gmt_annotations_fgsea,
               '_rankby_', rankby, '.rds'))
      if (do.fgsea) {
        if (!file.exists(fgsea_file) || force) {
          fgseaResAll <- list()
          gmt.file    <- file.path(gsea_db_dir, org, gmt_annotations_fgsea)
          pathways    <- gmtPathways(gmt.file)
          pathways    <- lapply(pathways,
            function(x) x[!is.na(x) & x != 'NA' & x != ''])
          fgseaResCluster <- list()
          fgseaRanks      <- list()
          for (i in unique(deg2$cluster)) {
            deg2Cluster <- deg2[deg2$cluster == i, ]
            fgsea_ranks <- if (rankby == 'LFC') deg2Cluster$avg_log2FC else
              -log10(deg2Cluster$p_val + 1e-300) * sign(deg2Cluster$avg_log2FC)
            fgsea_ranks <- setNames(fgsea_ranks, deg2Cluster$gene)
            fgsea_ranks <- fgsea_ranks[fgsea_ranks != 0]
            tryCatch({
              fgseaRes    <- fgseaMultilevel(pathways, fgsea_ranks,
                               minSize = 15, maxSize = 500, BPPARAM = NULL)
              fgseaResCol <- collapsePathways(fgseaRes, stats = fgsea_ranks,
                               pathway = pathways)
              fgseaResCluster[[i]] <- fgseaRes[
                fgseaRes$pathway %in% fgseaResCol$mainPathways]
              fgseaRanks[[gmt_annotations_fgsea]][[i]] <- fgsea_ranks
            }, error = function(e)
              message(paste('fGSEA failed for cluster', i, ':', e$message)))
          }
          fgseaResAll[[gmt_annotations_fgsea]] <- fgseaResCluster
          saveRDS(fgseaResAll, fgsea_file)
          saveRDS(fgseaRanks, file.path(projdir_deg2,
            paste0('fgsea_ranks_', gmt_annotations_fgsea, '.rds')))
        } else {
          fgseaResAll <- readRDS(fgsea_file)
          fgseaRanks  <- readRDS(file.path(projdir_deg2,
            paste0('fgsea_ranks_', gmt_annotations_fgsea, '.rds')))
        }
        if (length(fgseaResAll[[1]]) > 1) {
          fgsea_dp <- lapply(fgseaResAll, function(y)
            dotGSEA(y, padj_threshold = pvalAdjTrheshold, type = 'fgsea',
                    top_pathways = top_pathways,
                    cluster_rows = TRUE, cluster_cols = TRUE))
          lapply(seq_along(fgsea_dp), function(x) {
            if (!is.null(fgsea_dp[[x]])) {
              pdf(file.path(projdir_deg2, 'Plots',
                  paste0('fGSEA_', names(fgsea_dp)[x], '_dotplots.pdf')),
                  width  = 8,
                  height = 3 +
                    length(unique(fgsea_dp[[x]]$data$pathway)) / 7)
              print(fgsea_dp[[x]])
              dev.off()
            }
          })
        }
      }

      # Violin plots: selected genes in Atf6+ cells
      srt_atf6 <- srt[, srt@meta.data[, atf6_pos_col] == atf6_col]
      atf6_deg_genes <- deg2[deg2$cluster == atf6_col, ]
      atf6_deg_genes <- atf6_deg_genes[order(atf6_deg_genes$avg_log2FC), ]

      gene_selection_atf6 <- c('Gabbr2', 'Mapt', 'Kcnj3', 'Kcnab1',
                                'Kcnma1', 'Kcnq5')
      gene_selection_atf6 <- gene_selection_atf6[
        gene_selection_atf6 %in% rownames(srt_atf6)]

      if (length(gene_selection_atf6) > 0) {
        ccomp_df_atf6 <- srt_atf6@meta.data
        genes_mat_atf6 <- as.data.frame(
          t(srt_atf6@assays$RNA@data[gene_selection_atf6, , drop = FALSE]))
        ccomp_df_atf6 <- cbind(ccomp_df_atf6, genes_mat_atf6)
        ccomp_df_atf6 <- ccomp_df_atf6[
          ccomp_df_atf6$condgen %in% c('SNI WT', 'SNI KO'), ]
        ccomp_df_atf6$condgen <- factor(as.character(ccomp_df_atf6$condgen),
                                         levels = c('SNI WT', 'SNI KO'))
        # sanitize column names
        colnames(ccomp_df_atf6) <- gsub('-', '_', colnames(ccomp_df_atf6))
        gene_selection_atf6 <- gsub('-', '_', gene_selection_atf6)

        box_atf6 <- list()
        stat_atf6_store <- list()
        for (g in gene_selection_atf6) {
          if (!g %in% colnames(ccomp_df_atf6)) next
          ccomp2 <- ccomp_df_atf6[, c('condgen', g)]
          box_atf6[[g]] <- ggplot(ccomp2,
            aes_string(x = 'condgen', y = g)) +
            geom_violin(trim = TRUE, aes_string(fill = 'condgen'),
                        alpha = 0.7, lwd = .2) +
            gtheme +
            scale_fill_manual(values = palette_condgen)
          stat_g <- box_atf6[[g]]$data %>%
            rstatix::wilcox_test(reformulate('condgen', g)) %>%
            adjust_pvalue(method = "none") %>%
            add_significance()
          stat_atf6_store[[g]] <- stat_g
          stat_g <- stat_g %>% add_xy_position(x = 'condgen',
                                                step.increase = 0.1)
          box_atf6[[g]] <- box_atf6[[g]] +
            stat_pvalue_manual(stat_g, remove.bracket = FALSE,
                               bracket.nudge.y = 0.4, hide.ns = FALSE,
                               label = "p.adj.signif") + NoLegend()
        }

        write.csv(as.data.frame(do.call(rbind, stat_atf6_store)),
                  'SCENIC_Atf6_NP_selected_genes.csv')

        dp_atf6 <- DimPlot(srt, group.by = atf6_pos_col, combine = FALSE)
        dp_atf6[[1]] <- dp_atf6[[1]] + theme_void()
        fp_atf6 <- fp(srt, atf6_col)

        pdf(file.path('Plots', paste0('SCENIC_', atf6_col, '_module_pos_umap.pdf')),
            height = 3.5, width = 20)
        print(wrap_plots(c(dp_atf6, fp_atf6)))
        dev.off()

        pdf(file.path('Plots', paste0('SCENIC_', atf6_col,
            '_module_deg_vlnplots.pdf')), height = 4, width = 5)
        print(wrap_plots(box_atf6))
        dev.off()
      }
    }
  } else {
    message('No Atf6 column found in srt after adding auc_mtx – skipping SCENIC section')
  }
} else {
  message(paste('SCENIC auc_mtx.csv not found at', scenic_auc_csv,
                '– skipping SCENIC section'))
}


### ===== 10. WGCNA GREEN MODULE ANALYSIS ==================================

if (file.exists(wgcna_modules_csv)) {
  message('Loading WGCNA module_assignments.csv ...')
  modulesL <- read.csv(wgcna_modules_csv)
  wgcna_module_score_ngenes <- 100
  modulesL <- split(modulesL$gene_name, modulesL$module)
  modulesL <- lapply(modulesL, function(x) head(x, wgcna_module_score_ngenes))

  message('Running AddModuleScore for WGCNA modules ...')
  srt <- ModScoreCor(seurat_obj = srt, geneset_list = modulesL,
                     cor_threshold = NULL, pos_threshold = NULL,
                     listName = 'wgcna_modules',
                     outdir = file.path(projdir, 'Plots'))

  if ('green' %in% colnames(srt@meta.data)) {
    mod <- 'green'
    metaGroupNames_wg <- c('sampleID', 'condgen')

    ## UMAP feature plot
    umap_green <- FeaturePlot(srt, features = mod, ncol = 1,
                              combine = FALSE, pt.size = 1,
                              reduction = reductionName)
    umap_green[[1]] <- umap_green[[1]] + theme_void() + NoAxes() +
      ggtitle(mod) +
      theme(plot.title = element_text(size = 8)) +
      scale_colour_gradientn(colours = viridis::viridis(100),
                             limits = c(0, max(umap_green[[1]]$data[, 4])))

    pdf('Plots/WGCNA_green_distribution_umap.pdf', height = 2.7, width = 3)
    print(umap_green[[1]])
    dev.off()

    ## Boxplots per condgen (pseudo-bulk per sample)
    ccomp_df_wg <- srt@meta.data[, c(mod, metaGroupNames_wg), drop = FALSE]
    ccomp_df_wg <- aggregate(ccomp_df_wg[, mod, drop = FALSE],
      by = as.list(ccomp_df_wg[, metaGroupNames_wg, drop = FALSE]), mean)

    sample_rep_wg <- table(srt@meta.data[, 'condgen'])
    keep_group_wg <- names(sample_rep_wg[sample_rep_wg >= 2])

    box_wg <- ggplot(ccomp_df_wg, aes_string(x = 'condgen', y = mod)) +
      geom_boxplot(aes_string(fill = 'condgen'), outlier.size = 2,
                   outlier.alpha = 0.2, notch = FALSE, lwd = .2) +
      scale_fill_manual(values = palette_condgen) +
      ggtitle(mod) + ylab('green module score') + gtheme

    st_wg <- box_wg$data[box_wg$data$condgen %in% keep_group_wg, ] %>%
      t_test(reformulate('condgen', mod)) %>%
      adjust_pvalue(method = "none") %>%
      add_significance()
    st_wg <- st_wg %>% add_xy_position(x = 'condgen', step.increase = 0.1)
    box_wg_p <- box_wg + stat_pvalue_manual(st_wg, remove.bracket = FALSE,
      bracket.nudge.y = 0, hide.ns = TRUE, label = "p.adj.signif") + NoLegend()

    pdf('Plots/WGCNA_green_distribution_boxplots.pdf', height = 2.7, width = 1.7)
    print(box_wg_p)
    dev.off()

    ## DEG2 – green module+ vs rest: SNI_RGS4_KO vs SNI_RGS4_WT -----------
    green_threshold  <- 0.0
    green_pos_col    <- paste0(mod, '_pos')
    srt@meta.data[, green_pos_col] <- ifelse(
      srt@meta.data[, mod] > green_threshold, mod, 'rest')

    force            <- FALSE
    do.fgsea         <- TRUE
    rankby           <- 'LFC'
    logfcThreshold   <- 0
    pvalAdjTrheshold <- 0.05
    topGenes         <- 50
    pval_column      <- 'pval'
    metaGroupName1   <- green_pos_col
    metaGroupName2   <- 'condition'
    deg2Ident        <- c('SNI_RGS4_KO', 'SNI_RGS4_WT')
    top_pathways     <- 5
    gmt_annotations_fgsea <- 'h.all.v7.1.symbol.gmt'

    projdir_deg2_wg <- paste0('deg2_', metaGroupName2, '_', mod,
                               '_pos_', paste(deg2Ident, collapse = '_vs_'),
                               '_lfc_', logfcThreshold)
    dir.create(file.path(projdir_deg2_wg, 'Plots'), recursive = TRUE,
               showWarnings = FALSE)

    min.pct <- 0.1
    if (!file.exists(file.path(projdir_deg2_wg, 'DEG2_genes.rds')) || force) {
      deg2L_wg <- list()
      Idents(srt) <- srt@meta.data[, metaGroupName1]
      for (i in unique(srt@meta.data[, metaGroupName1])) {
        nigenes <- NI.genes(rownames(srt),
          ni.goterms = c('ribosome biogenesis', 'ncRNA processing',
                         'RNA processing', 'oxidation.reduction process'),
          ni.genes   = c('^Mt', '^Rpl', '^mt-', '^Rps'))
        srt2 <- srt[!rownames(srt) %in% nigenes, ]
        if ('layers' %in% slotNames(srt@assays$RNA)) {
          count_mat <- srt@assays$RNA@layers$counts
          rownames(count_mat) <- rownames(srt)
        } else {
          count_mat <- srt2@assays$RNA@counts
        }
        is.expressed <- apply(count_mat, 1,
                              function(x) sum(x > 0) / length(x))
        srt2 <- srt2[is.expressed >= min.pct, ]
        tryCatch({
          deg2Cluster <- FindMarkers(srt2,
            ident.1 = deg2Ident[1], ident.2 = deg2Ident[2],
            group.by = metaGroupName2, subset.ident = i,
            only.pos = FALSE, min.pct = min.pct,
            logfc.threshold = logfcThreshold, test.use = 'wilcox')
          if (nrow(deg2Cluster) < 1) stop('no genes')
          deg2Cluster$cluster <- i
          deg2Cluster$gene    <- rownames(deg2Cluster)
          deg2L_wg[[i]] <- deg2Cluster
        }, error = function(e)
          message(paste('Skipping cluster', i, ':', e$message)))
      }
      deg2_wg <- do.call(rbind, deg2L_wg)
      if (!is.null(deg2_wg) && pval_column == 'pval')
        deg2_wg$p_val_adj <- deg2_wg$p_val
      write.csv(
        deg2_wg[deg2_wg$p_val_adj < pvalAdjTrheshold, ],
        file.path(projdir_deg2_wg, 'DEG2_genes.csv'))
      saveRDS(deg2_wg, file.path(projdir_deg2_wg, 'DEG2_genes.rds'))
    } else {
      message('DEG2-WGCNA-green object found!')
      deg2_wg <- readRDS(file.path(projdir_deg2_wg, 'DEG2_genes.rds'))
    }

    if (!is.null(deg2_wg) && nrow(deg2_wg) > 0) {
      deg2_wg_sig <- deg2_wg[deg2_wg$p_val_adj < pvalAdjTrheshold, ]
      clusters_idx_wg <- unique(srt@meta.data[, metaGroupName1])

      if (nrow(deg2_wg_sig) > 0 && any(deg2_wg_sig$avg_log2FC != 0)) {
        deg_sum_pos_wg <- as.data.frame(
          table(deg2_wg_sig$cluster[deg2_wg_sig$avg_log2FC > 0]))
        rownames(deg_sum_pos_wg) <- deg_sum_pos_wg$Var1
        deg_sum_pos_wg <- deg_sum_pos_wg[clusters_idx_wg, ]
        rownames(deg_sum_pos_wg) <- clusters_idx_wg
        deg_sum_neg_wg <- as.data.frame(
          table(deg2_wg_sig$cluster[deg2_wg_sig$avg_log2FC < 0]))
        deg_sum_neg_wg$Freq <- -1 * deg_sum_neg_wg$Freq
        rownames(deg_sum_neg_wg) <- deg_sum_neg_wg$Var1
        deg_sum_neg_wg <- deg_sum_neg_wg[clusters_idx_wg, ]
        rownames(deg_sum_neg_wg) <- clusters_idx_wg
        deg_sum_wg <- cbind(deg_sum_pos_wg, deg_sum_neg_wg)
        colnames(deg_sum_wg)[3:4] <- c('Var1_neg', 'Freq_neg')

        if (mod %in% rownames(deg_sum_wg)) {
          deg_sum_wg_green <- deg_sum_wg[mod, , drop = FALSE]
          deg_bar_wg <- ggplot(deg_sum_wg_green,
            aes(x = rownames(deg_sum_wg_green))) +
            geom_bar(aes(y = Freq), fill = 'red', color = 'black',
                     stat = "identity", alpha = 0.8) +
            geom_bar(aes(y = Freq_neg), fill = 'blue', color = 'black',
                     stat = "identity", position = "identity", alpha = 0.8) +
            ylim(c(-max(abs(na.omit(c(deg_sum_wg_green$Freq_neg,
                                      deg_sum_wg_green$Freq)))),
                    max(na.omit(abs(c(deg_sum_wg_green$Freq_neg,
                                      deg_sum_wg_green$Freq)))))) +
            labs(title = "DEG green module cluster", y = "Value") +
            gtheme
          pdf(file.path(projdir_deg2_wg, 'Plots',
              'WGCNA_green_number_deg_cluster.pdf'), height = 4, width = 1)
          print(deg_bar_wg)
          dev.off()
        }
      }

      if (length(unique(deg2_wg$cluster)) > 1 && nrow(deg2_wg_sig) > 0) {
        dhm_wg <- diffClustHeat(deg_genes_df = deg2_wg, sort_by = 'p_val_adj',
          topGenes = topGenes, pvalAdjTrheshold = pvalAdjTrheshold,
          col_limit = NULL,
          plotcol = rev(RColorBrewer::brewer.pal(3, 'RdBu')),
          name   = paste(deg2Ident, collapse = '-'),
          cluster_columns = TRUE, cluster_rows = TRUE, border = TRUE)
        pdf(file.path(projdir_deg2_wg, 'Plots',
            paste0('Top_', topGenes, '_genes_heatmap.pdf')),
            width = 3 + ncol(dhm_wg@matrix) / 15,
            height = 3 + nrow(dhm_wg@matrix) / 20)
        print(dhm_wg)
        dev.off()
      }

      # Volcano plots
      vp_wg <- list()
      for (cl in unique(deg2_wg$cluster)) {
        deg2_cl         <- deg2_wg[deg2_wg$cluster == cl, ]
        deg2_cl$label   <- ''
        deg2_cl$label[abs(deg2_cl$avg_log2FC) > 0.5 &
                        deg2_cl$p_val_adj < pvalAdjTrheshold] <-
          deg2_cl$gene[abs(deg2_cl$avg_log2FC) > 0.5 &
                         deg2_cl$p_val_adj < pvalAdjTrheshold]
        deg2_cl$Direction <- ''
        deg2_cl$Direction[deg2_cl$avg_log2FC >  0.5 &
                            deg2_cl$p_val_adj < pvalAdjTrheshold] <- deg2Ident[1]
        deg2_cl$Direction[deg2_cl$avg_log2FC < -0.5 &
                            deg2_cl$p_val_adj < pvalAdjTrheshold] <- deg2Ident[2]
        vp_wg[[cl]] <- ggplot(deg2_cl, aes(x = avg_log2FC,
                                            y = -log10(p_val_adj))) +
          geom_point(size = 1, shape = 19, aes(color = Direction)) +
          geom_vline(xintercept =  0.5, linetype = "dotted",
                     color = "blue", size = 1) +
          geom_vline(xintercept = -0.5, linetype = "dotted",
                     color = "blue", size = 1) +
          geom_hline(yintercept = -log10(pvalAdjTrheshold),
                     linetype = "dotted", color = "blue", size = 1) +
          geom_text_repel(size = 2, data = deg2_cl, aes(label = label)) +
          ggtitle(cl) +
          scale_color_manual(values = c('grey77', 'red', 'green')) +
          theme_light()
      }
      png(file.path(projdir_deg2_wg, 'Plots', 'volcano_plots.png'),
          width = 4900, height = 2000, res = 300)
      print(wrap_plots(vp_wg), ncol = 3)
      dev.off()

      # Violin plots: selected genes in green module+ cells
      srt_green <- srt[, srt@meta.data[, green_pos_col] == mod]
      gene_selection_wg <- c('Syn3', 'Scn8a', 'Actg1', 'Tubb2b')
      gene_selection_wg <- gene_selection_wg[
        gene_selection_wg %in% rownames(srt_green)]

      if (length(gene_selection_wg) > 0) {
        ccomp_df_wg2 <- srt_green@meta.data
        genes_mat_wg <- as.data.frame(
          t(srt_green@assays$RNA@data[gene_selection_wg, , drop = FALSE]))
        ccomp_df_wg2 <- cbind(ccomp_df_wg2, genes_mat_wg)
        ccomp_df_wg2 <- ccomp_df_wg2[
          ccomp_df_wg2$condgen %in% c('SNI WT', 'SNI KO'), ]
        ccomp_df_wg2$condgen <- factor(as.character(ccomp_df_wg2$condgen),
                                        levels = c('SNI WT', 'SNI KO'))

        box_wg2 <- list()
        stat_wg_store <- list()
        for (g in gene_selection_wg) {
          if (!g %in% colnames(ccomp_df_wg2)) next
          ccomp_wg_g <- ccomp_df_wg2[, c('condgen', g)]
          box_wg2[[g]] <- ggplot(ccomp_wg_g,
            aes_string(x = 'condgen', y = g)) +
            geom_violin(trim = TRUE, aes_string(fill = 'condgen'),
                        alpha = 0.7, lwd = .2) +
            gtheme +
            scale_fill_manual(values = palette_condgen)
          stat_g_wg <- box_wg2[[g]]$data %>%
            rstatix::wilcox_test(reformulate('condgen', g)) %>%
            adjust_pvalue(method = "none") %>%
            add_significance()
          stat_wg_store[[g]] <- stat_g_wg
          stat_g_wg <- stat_g_wg %>%
            add_xy_position(x = 'condgen', step.increase = 0.1)
          box_wg2[[g]] <- box_wg2[[g]] +
            stat_pvalue_manual(stat_g_wg, remove.bracket = FALSE,
                               bracket.nudge.y = 0.4, hide.ns = FALSE,
                               label = "p.adj.signif") + NoLegend()
        }

        write.csv(as.data.frame(do.call(rbind, stat_wg_store)),
                  'WGCNA_green_NP_selected_genes.csv')

        dp_green <- DimPlot(srt, group.by = green_pos_col, combine = FALSE)
        dp_green[[1]] <- dp_green[[1]] + theme_void()
        fp_green <- fp(srt, mod)

        pdf(file.path('Plots', paste0('WGCNA_', mod,
            '_module_pos_umap.pdf')), height = 3.5, width = 20)
        print(wrap_plots(c(dp_green, fp_green)))
        dev.off()

        pdf(file.path('Plots', paste0('WGCNA_', mod,
            '_module_deg_vlnplots.pdf')), height = 2, width = 7)
        print(wrap_plots(box_wg2, ncol = 4))
        dev.off()
      }

      # fGSEA for WGCNA green
      fgsea_file_wg <- file.path(projdir_deg2_wg,
        paste0('fGSEA_annotation_', gmt_annotations_fgsea,
               '_rankby_', rankby, '.rds'))
      if (do.fgsea) {
        if (!file.exists(fgsea_file_wg) || force) {
          fgseaResAll_wg <- list()
          gmt.file  <- file.path(gsea_db_dir, org, gmt_annotations_fgsea)
          pathways  <- gmtPathways(gmt.file)
          pathways  <- lapply(pathways,
            function(x) x[!is.na(x) & x != 'NA' & x != ''])
          fgseaResCluster_wg <- list()
          fgseaRanks_wg      <- list()
          for (i in unique(deg2_wg$cluster)) {
            deg2Cluster_wg <- deg2_wg[deg2_wg$cluster == i, ]
            fgsea_ranks_wg <- if (rankby == 'LFC')
              deg2Cluster_wg$avg_log2FC else
              -log10(deg2Cluster_wg$p_val + 1e-300) *
                sign(deg2Cluster_wg$avg_log2FC)
            fgsea_ranks_wg <- setNames(fgsea_ranks_wg, deg2Cluster_wg$gene)
            fgsea_ranks_wg <- fgsea_ranks_wg[fgsea_ranks_wg != 0]
            tryCatch({
              fgseaRes_wg    <- fgseaMultilevel(pathways, fgsea_ranks_wg,
                                 minSize = 15, maxSize = 500, BPPARAM = NULL)
              fgseaResCol_wg <- collapsePathways(fgseaRes_wg,
                                 stats = fgsea_ranks_wg, pathway = pathways)
              fgseaResCluster_wg[[i]] <- fgseaRes_wg[
                fgseaRes_wg$pathway %in% fgseaResCol_wg$mainPathways]
              fgseaRanks_wg[[gmt_annotations_fgsea]][[i]] <- fgsea_ranks_wg
            }, error = function(e)
              message(paste('fGSEA-green failed for cluster', i, ':',
                            e$message)))
          }
          fgseaResAll_wg[[gmt_annotations_fgsea]] <- fgseaResCluster_wg
          saveRDS(fgseaResAll_wg, fgsea_file_wg)
          saveRDS(fgseaRanks_wg, file.path(projdir_deg2_wg,
            paste0('fgsea_ranks_', gmt_annotations_fgsea, '.rds')))
        } else {
          fgseaResAll_wg <- readRDS(fgsea_file_wg)
        }
        if (length(fgseaResAll_wg[[1]]) > 1) {
          fgsea_dp_wg <- lapply(fgseaResAll_wg, function(y)
            dotGSEA(y, padj_threshold = pvalAdjTrheshold, type = 'fgsea',
                    top_pathways = top_pathways,
                    cluster_rows = TRUE, cluster_cols = TRUE))
          lapply(seq_along(fgsea_dp_wg), function(x) {
            if (!is.null(fgsea_dp_wg[[x]])) {
              pdf(file.path(projdir_deg2_wg, 'Plots',
                  paste0('fGSEA_', names(fgsea_dp_wg)[x], '_dotplots.pdf')),
                  width  = 8,
                  height = 3 +
                    length(unique(fgsea_dp_wg[[x]]$data$pathway)) / 7)
              print(fgsea_dp_wg[[x]])
              dev.off()
            }
          })
        }
      }
    }
  } else {
    message('No "green" column found after WGCNA module scoring – skipping')
  }
} else {
  message(paste('WGCNA module_assignments.csv not found at', wgcna_modules_csv,
                '– skipping WGCNA section'))
}


### ===== 11. MRGPRD / CALCA VIOLIN + DOT PLOTS ===========================

gene_pair <- c('Mrgprd', 'Calca')
gene_pair <- gene_pair[gene_pair %in% rownames(srt)]

if (length(gene_pair) > 0) {
  contrasts <- list(c('SNI_RGS4_WT', 'SNI_RGS4_KO'),
                    c('SNI_RGS4_WT', 'SHAM_RGS4_WT'))

  # Dotplots
  pdf(file.path('Plots',
      paste0(paste(gene_pair, collapse = '_'), '_dotplot.pdf')),
      height = 5)
  tryCatch(
    print(DotPlot(srt[, srt$condition %in% c('SNI_RGS4_WT', 'SHAM_RGS4_WT')],
                  features = gene_pair, group.by = 'condition')),
    error = function(e) message(e$message))
  tryCatch(
    print(DotPlot(srt[, srt$condition %in% c('SNI_RGS4_WT', 'SNI_RGS4_KO')],
                  features = gene_pair, group.by = 'condition')),
    error = function(e) message(e$message))
  dev.off()

  # Per-cell violin plots per contrast
  boxL2_all <- list()
  for (contrast in contrasts) {
    srt_sub    <- srt[, srt$condition %in% contrast]
    ccomp_df_c <- srt_sub@meta.data
    genes_mat_c <- as.data.frame(
      t(srt_sub@assays$RNA@data[gene_pair, , drop = FALSE]))
    ccomp_df_c  <- cbind(ccomp_df_c, genes_mat_c)
    ccomp_df_c$condition <- factor(ccomp_df_c$condition, levels = contrast)

    boxL <- list()
    for (g in gene_pair) {
      ccomp1 <- ccomp_df_c[ccomp_df_c[, g] > 0, ]
      box_c  <- ggplot(ccomp1, aes_string(x = 'condition', y = g)) +
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
        stat_pvalue_manual(st, remove.bracket = FALSE,
                           bracket.nudge.y = 0, hide.ns = FALSE,
                           label = "p.adj.signif") + NoLegend()
    }
    boxL2_all[[paste(contrast, collapse = '_')]] <- boxL

    pdf(file.path('Plots', paste0(paste(gene_pair, collapse = '_'), '_',
        paste(contrast, collapse = '_'), '_vlnPlot.pdf')), height = 5)
    print(wrap_plots(boxL))
    dev.off()
  }

  # Per-sample mean plots per contrast
  for (contrast in contrasts) {
    srt_sub    <- srt[, srt$condition %in% contrast]
    ccomp_df_c <- srt_sub@meta.data
    genes_mat_c <- as.data.frame(
      t(srt_sub@assays$RNA@data[gene_pair, , drop = FALSE]))
    ccomp_df_c  <- cbind(ccomp_df_c, genes_mat_c)

    ccomp_agg <- aggregate(ccomp_df_c[, gene_pair, drop = FALSE],
      by = as.list(ccomp_df_c[, c('sampleID', 'condition'), drop = FALSE]),
      mean)
    ccomp_agg <- na.omit(ccomp_agg)
    ccomp_agg$condition <- factor(ccomp_agg$condition, levels = contrast)

    boxL_s <- list()
    for (g in gene_pair) {
      ccomp_s <- ccomp_agg[ccomp_agg[, g] > 0, ]
      box_s   <- ggplot(ccomp_s, aes_string(x = 'condition', y = g)) +
        geom_violin(trim = TRUE, aes_string(fill = 'condition')) +
        geom_boxplot(aes_string(fill = 'condition')) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      st_s <- box_s$data %>%
        t_test(reformulate('condition', g)) %>%
        adjust_pvalue(method = "none") %>%
        add_significance()
      st_s <- st_s %>% add_xy_position(x = 'condition', step.increase = 0.1)
      boxL_s[[g]] <- box_s +
        stat_pvalue_manual(st_s, remove.bracket = FALSE,
                           bracket.nudge.y = 0, hide.ns = FALSE,
                           label = "p.adj.signif") + NoLegend()
    }

    pdf(file.path('Plots', paste0(paste(gene_pair, collapse = '_'), '_',
        paste(contrast, collapse = '_'), '_vlnPlot_per_sample.pdf')),
        height = 5)
    print(wrap_plots(boxL_s))
    dev.off()
  }
}

message('NP script completed successfully.')
