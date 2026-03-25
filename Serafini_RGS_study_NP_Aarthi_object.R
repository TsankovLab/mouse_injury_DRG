use UGER
conda activate scrnatools 
R

#data.dir = '/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_1_v1/filtered_feature_bc_matrix/'
#samples = c(paste0('GFP_',1:3), paste0('CB2_PTENL_',1:3), paste0('CB2_PTENL_C124S_',1:3))
#samples_path1 = paste0('/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_1_v1/',samples,'/gex/')
#samples_path2 = paste0('/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_2_v1/',samples,'/gex/')
samples_path = c(
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-SNI-Rgs4KO-1/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-SNI-Rgs4KO-2/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-SNI-Rgs4KO-3/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-SNI-Rgs4KO-4/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-SNI-Rgs4WT-1/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-SNI-Rgs4WT-2/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-SNI-Rgs4WT-3/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-SNI-Rgs4WT-4/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-Sham-Rgs4WT-1/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-Sham-Rgs4WT-2/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-Sham-Rgs4WT-3/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-Sham-Rgs4WT-4/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-Sham-Rgs4KO-1/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-Sham-Rgs4KO-2/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-Sham-Rgs4KO-3/outs/raw_feature_bc_matrix/',
'/broad/hptmp/bgiotti/aarthi/ramaka02.dmz.hpc.mssm.edu/mouse_pain_snRNA/cellranger_210329/DRG-Sham-Rgs4KO-4/outs/raw_feature_bc_matrix/')

#meta = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/KaLung_collab/metadata_scRNA-Seq_IECs_Hdac4.csv')
meta = data.frame (sampleID = c(
'DRG-SNI-Rgs4KO-1',
'DRG-SNI-Rgs4KO-2',
'DRG-SNI-Rgs4KO-3',
'DRG-SNI-Rgs4KO-4',
'DRG-SNI-Rgs4WT-1',
'DRG-SNI-Rgs4WT-2',
'DRG-SNI-Rgs4WT-3',
'DRG-SNI-Rgs4WT-4',
'DRG-Sham-Rgs4WT-1',
'DRG-Sham-Rgs4WT-2',
'DRG-Sham-Rgs4WT-3',
'DRG-Sham-Rgs4WT-4',
'DRG-Sham-Rgs4KO-1',
'DRG-Sham-Rgs4KO-2',
'DRG-Sham-Rgs4KO-3',
'DRG-Sham-Rgs4KO-4'),
 condition = c(
  'SNI_RGS4_KO','SNI_RGS4_KO',
  'SNI_RGS4_KO','SNI_RGS4_KO',
  'SNI_RGS4_WT','SNI_RGS4_WT',
  'SNI_RGS4_WT','SNI_RGS4_WT',
  'SHAM_RGS4_WT','SHAM_RGS4_WT',
  'SHAM_RGS4_WT','SHAM_RGS4_WT',
  'SHAM_RGS4_KO','SHAM_RGS4_KO',
'SHAM_RGS4_KO','SHAM_RGS4_KO'))

##### Import annotated seurat object from Aarthi  #####
#srt2 = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Serafini_prj/GSE154659_C57_Raw_counts.rds')
rerun=FALSE
if (rerun)
  {
  srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Serafini_prj/mouse_pain_study2_annotated.rds')
  srt = SplitObject(srt, split.by = "orig.ident")
  srt = srt[meta$sampleID]
  }

#meta = meta[!duplicated(meta$assignment),]
#samples_path = '/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_1_v1/raw_feature_bc_matrix/'
# Set project directory
proj_name = 'nerve_injury_RGS4_aarthi'

# Set cellranger_to_seurat parameters
cr_to_seurat = list(
  run_cellbender=F,
  run_cellbender = FALSE, # if this is set to TRUE then input raw (and not filtered) cellranger count matrices
  cellbender_samples = NULL,
  #cellbender_parameters = list(p811 = c(expected_cells = 12000,total_droplets_included = 100000, low_count_threshold = 5)),
  cellbender_parameters = NULL,
  org = 'mouse',
  datatype = 'RNA',
  cr_output = 'filtered',
  samples_path = samples_path,
  meta = meta, 
  is.hashed = FALSE,
  srt_provided = NULL
  )
# Set QC parameters
qc_params = list(
  filtering = 'hard', # 'emptyDrops' or 'hard' filtering
  nFeat = 0,#400, # Number of features per cells. default 400
  nCounts = 0, #1000, # Number of UMI per cell. Default 800
  pchM = 100, # Percent mitochondrial genes. Default 25 
  remove.samples = NULL, # Remove bad samples. Takes vector of sampleIDs 
  processInd = FALSE # run preprocessing per sample before running it on the merged data
  )
### Data processing and clustering variables ###
harmony_params = list(
  batch = 'no'
  )

data_processing_param = list(
  variablefeatures = 'seurat',
  nfeat = 2000, # number of variable genes to consider for dimentionality reduction
  sigPCs = 15,
  vars_to_regress='nCount_RNA',
  metaGroupNames = c('sampleID','condition'),
  res = c(0.2, 0.8, 2, 5,10,15) # denovo cluster resolutions 
  )

# Initiate pipeline
force = FALSE # re run pipeline from the beginning no matter if objects are found
subclustername = 'NP'
projdir_init = paste0 ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Serafini_prj/",proj_name,"_analysis")
scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'
source (file.path(scrna_pipeline_dir,'master_scrna.R'))

srt$Condition = ifelse (grepl ('SHAM',srt$condition), 'Sham','SNI')
srt$condgen = paste0 (srt$Condition, '_', srt$genotype)
srt$condgen = gsub ('_',' ', srt$condgen)
srt$condgen = factor (srt$condgen, levels = c('Sham WT','Sham KO','SNI WT','SNI KO'))
srt$condgen2 = paste0 (srt$Condition, '_', srt$genotype)
srt$condgen2 = factor (srt$condgen2, levels = c('Sham_WT','Sham_KO','SNI_WT','SNI_KO'))

names(palette_condgen) = gsub ('_',' ', names(palette_condgen))
palette_condgen = c('Sham WT' = 'grey','Sham KO' = 'firebrick1', 'SNI WT'= 'darkslategrey', 'SNI KO' = 'violetred4')

gtheme = theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),          
      axis.line = element_blank(),  # Remove axis lines
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )

# # DEG ####
enricher_universe = 'all'
logfcThreshold = .25
pvalAdjTrheshold = 0.01
metaGroupName = paste0(reductionGraphSnn,'_res.0.8')
top_pathways = 5
force = FALSE
top_genes = 5
source (paste0(scrna_pipeline_dir,'DEG_standard.R'))

metaGroupName1 = paste0(reductionGraphSnn,'_res.0.8')
metaGroupName2 = 'condition'
cc_box = cellComp (
  seurat_obj = srt, 
  metaGroups = c('sampleID',metaGroupName2,metaGroupName1),
  plot_as = 'box',
  #pal = genotype_pal,
  facet_ncol = 8,
  ) + 
  gtheme

cc_df = cc_box$data
stat.test <- cc_df %>%
  group_by_at (metaGroupName1) %>%
  t_test (reformulate (metaGroupName2, 'Freq')) %>%
  adjust_pvalue (method = "none") %>%
  add_significance ()
stat.test <- stat.test %>% add_xy_position (x = metaGroupName1, step.increase=0.1)

png (paste0 ('Plots/cell_composition_',metaGroupName1,'_barboxplots2.png'), width=2000, height=1000, res=300)
(cc_box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
   bracket.nudge.y = 0, hide.ns = TRUE,
    label = "p.adj.signif"))# + plot_layout (widths= c(1,1,6))
dev.off()

pdf (paste0 ('Plots/cell_composition_',metaGroupName1,'_barboxplots2.pdf'), width=10, height=6)
(cc_box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
   bracket.nudge.y = 0, hide.ns = TRUE,
    label = "p.adj.signif"))# + plot_layout (widths= c(1,1,6))
dev.off()



### Run WGCNA ####
# Set variables
force=FALSE # force re-running WGCNA 
do.fgsea=TRUE
powerTable=FALSE # plot the powerTable plot. Can take a while to generate
do.plots=TRUE
softPower=3 # Set the softPower
deepSplit=1 # Set this to increase decrease the number of modules idendified. 1-4
mergeCutHeight = 0.20 # Height below which two modules are merged
metacells_k = 20 # number of cells to create pseudobulks
max_shared = 15
metacells_groups = c('RNA_snn_res.0.2') # set metagroup in which to find metacells (usually your clustering / celltypes)
metaGroupNames = c('sampleID','celltype_predicted','condition') # set of metagroups for generating the boxplots
minModuleSize = 10
#genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat = 10000)) # Genes to use to compute WGCNA
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=3000))
enricher_universe = genes.keep # Genes to use as background for pathway enrichments analysis
source (file.path(scrna_pipeline_dir,'scWGCNA.R'))


### Run SCENIC ####
force = FALSE
motif_window = 'tss500bp'#'10kbp'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
source (paste0(scrna_pipeline_dir, 'SCENIC.R'))

# Run SCENIC plots ####
motif_window = 'tss500bp'#'10kbp'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
metaGroupNames = c('sampleID','celltype_predicted','condgen')
source (paste0(scrna_pipeline_dir, 'SCENIC_plots.R'))

# ccomp_df = srt@meta.data[,SCENIC_mods, drop=FALSE]
# ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,metaGroupNames,drop=F]), mean)
# ccomp_df = ccomp_df[ccomp_df$condgen %in% c('SNI_KO','SNI_WT'),]
# ccomp_df$condgen = as.factor (as.character (ccomp_df$condgen))
# box_p = lapply (seq(SCENIC_mods[[1]]), function(x) 
#     {
#     sample_rep = srt@meta.data
#     sample_rep = table (sample_rep[,metaGroupNames[2]], sample_rep[,metaGroupNames[3]])
#     keep_group = apply (sample_rep , 1, function(y) !any (y < 2))
#     keep_group = names(keep_group[keep_group])
#       if (length(keep_group) > 0 & metaGroupNames[1] != metaGroupNames[3])
#       {
#       box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods[x])) +
#           #geom_violin (trim=TRUE) +
#           geom_boxplot (aes_string(fill = metaGroupNames[3])) +
#           #geom_jitter (color="black", size=0.4, alpha=0.9) +
#           scale_fill_manual (values= palette_condgen) + 
#           ggtitle (SCENIC_mods[x]) + 
#           gthem
  
#       stat.test = box$data[box$data[,metaGroupNames[2]] %in% keep_group,] %>%
#         group_by_at (metaGroupNames[2]) %>%
#         t_test(reformulate (metaGroupNames[3], SCENIC_mods[x])) %>%
#         adjust_pvalue (method = "none") %>%
#         add_significance ()
#         box$data[,metaGroupNames[3]] = as.factor (as.character(box$data[,metaGroupNames[3]]))
#         stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
#         box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
#         bracket.nudge.y = 0, hide.ns = T,
#         label = "p.adj.signif") + NoLegend()
#       } else {
#       ccomp_df = srt@meta.data[,c(SCENIC_mods, metaGroupNames), drop=FALSE]
#       #ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)    
#       box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods[x])) +
#         geom_violin (trim=TRUE, aes_string (fill = metaGroupNames[3])) +
#         geom_boxplot (aes_string(fill = metaGroupNames[3])) +
#         #geom_bar (stats='identity') +
#         #geom_jitter (color="black", size=0.4, alpha=0.9) +
#         #scale_fill_manual (values= module_pal) + 
#         ggtitle (SCENIC_mods[x]) + 
#         gtheme

#         stat.test = box$data %>%
#         group_by_at (metaGroupNames[2]) %>%
#         wilcox_test(reformulate (metaGroupNames[3], SCENIC_mods[x])) %>%
#         adjust_pvalue (method = "none") %>%
#         add_significance ()
#         stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
        
#         box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
#         bracket.nudge.y = 0, hide.ns = T,
#         label = "p.adj.signif") + NoLegend()
#       }})

#   if (length(unique(ccomp_df[,metaGroupNames[2]])) > 1) box_p = lapply (seq_along(SCENIC_mods), function(x) box_p[[x]] + facet_wrap (as.formula(paste("~", metaGroupNames[2]))) + theme_classic() + 
#           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + NoLegend()) 
#   if (exists ('module_pal')) box_p = lapply (seq_along(SCENIC_mods), function(x) box_p[[x]] + scale_fill_manual (values= module_pal)) 
      
      
#     png (file.path ('Plots',paste0('SCENIC_module_scores_',paste(metaGroupNames,collapse='_'),'_umap.png')), width = 500, height = 1000, pointsize=10, res = 300, type="cairo")
#     box_p
#     dev.off()
      
#     pdf (file.path ('Plots',paste0('SCENIC_module_scores_',paste(metaGroupNames,collapse='_'),'_umap.pdf')), width = 5, height = 3)
#     box_p
#     dev.off()
      
  

  # Make UMAP of ATF6 SCENIC score
  metaGroupNames = c('sampleID','celltype_predicted','condgen')
  SCENIC_mods = colnames(srt@meta.data)[grepl ('SCENIC', colnames(srt@meta.data))]
  SCENIC_mods = 'SCENIC_Atf6'
umap_p1 = FeaturePlot (srt, 
    features = SCENIC_mods, 
    ncol=4, 
    combine = FALSE, 
    pt.size = 1, 
    reduction = reductionName)
for(i in 1:length(umap_p1)) 
  {
  umap_p1[[i]] = umap_p1[[i]] + 
    theme_void() + 
    NoAxes() +
    ggtitle (SCENIC_mods[[i]]) + 
    theme(plot.title = element_text(size = 8)) +
    #scale_colour_gradientn (colours = rev(brewer.pal(n = 11, name = "RdBu")),limits=c(-max (abs(p[[i]]$data[,4])), max (abs(p[[i]]$data[,4]))))
    scale_colour_gradientn (colours = viridis::viridis(100),limits=c(0, max (umap_p1[[i]]$data[,4])))
  }

pdf ('Plots/SCENIC_Atf6_distribution_umap2.pdf',height=2.7,width = 3)  
umap_p1
dev.off()


### Generate boxplots per meta groups
ccomp_df = srt@meta.data[,SCENIC_mods, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,metaGroupNames,drop=F]), mean)
#ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 

sample_rep = srt@meta.data
sample_rep = table (sample_rep[,metaGroupNames[2]], sample_rep[,metaGroupNames[3]])
keep_group = apply (sample_rep , 1, function(y) !any (y < 2))
keep_group = names(keep_group[keep_group])

box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods)) +
    #geom_violin (trim=TRUE) +
    geom_boxplot (aes_string(fill = metaGroupNames[3]),outlier.size=2,outlier.alpha = 0.2, notch=FALSE, lwd=.2) +
    #geom_jitter (color="black", size=0.4, alpha=0.9) +
    scale_fill_manual (values= palette_condgen) + 
    #scale_fill_manual (values= module_pal) + 
    ggtitle (SCENIC_mods) +
    ylab ('Atf6 regulon score') + 
   gtheme
  
stat.test = box$data[box$data[,metaGroupNames[2]] %in% keep_group,] %>%
  group_by_at (metaGroupNames[2]) %>%
  t_test(reformulate (metaGroupNames[3], SCENIC_mods)) %>%
  adjust_pvalue (method = "none") %>%
  add_significance ()
  stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
box_p = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
  bracket.nudge.y = 0, hide.ns = T,
  label = "p.adj.signif") + NoLegend()
        
pdf ('Plots/SCENIC_Atf6_distribution_boxplots.pdf',height=2.7,width = 1.7)  
box_p
dev.off()

ccomp_df = srt@meta.data[,SCENIC_mods, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,metaGroupNames,drop=F]), mean)
#ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 
ccomp_df = ccomp_df[ccomp_df$condgen %in% c('Sham WT','SNI WT'),]
ccomp_df$condgen = factor (as.character(ccomp_df$condgen), levels = c('Sham WT','SNI WT'))
box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods)) +
    #geom_violin (trim=TRUE) +
    geom_boxplot (aes_string(fill = metaGroupNames[3]),outlier.size=2,outlier.alpha = 0.2, notch=FALSE, lwd=.2) +
    #geom_jitter (color="black", size=0.4, alpha=0.9) +
    scale_fill_manual (values= palette_condgen) + 
    #scale_fill_manual (values= module_pal) + 
    ggtitle (SCENIC_mods) +
    ylab ('Atf6 regulon score') + 
   gtheme

stat.test = box$data[box$data[,metaGroupNames[2]] %in% keep_group,] %>%
  group_by_at (metaGroupNames[2]) %>%
  t_test(reformulate (metaGroupNames[3], SCENIC_mods)) %>%
  adjust_pvalue (method = "none") %>%
  add_significance ()
  stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
box_p = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
  bracket.nudge.y = 0, hide.ns = T,
  label = "p.adj.signif") + NoLegend()

pdf ('Plots/SCENIC_Atf6_distribution_GRANT_boxplots.pdf',height=2.7,width = 1.7)  
box_p
dev.off()

srtWT = RunUMAP (srt[, srt$condgen %in% c('Sham WT','SNI WT')], reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName)
umap_p1 = FeaturePlot (srtWT, 
    features = 'Atf6...', 
    ncol=4, 
    combine = FALSE, 
    pt.size = .01, 
    reduction = reductionName)
for(i in 1:length(umap_p1)) 
  {
  umap_p1[[i]] = umap_p1[[i]] + 
    theme_void() + 
    NoAxes() +
    #ggtitle (SCENIC_mods[[i]]) + 
    theme(plot.title = element_text(size = 8)) +
    #scale_colour_gradientn (colours = rev(brewer.pal(n = 11, name = "RdBu")),limits=c(-max (abs(p[[i]]$data[,4])), max (abs(p[[i]]$data[,4]))))
    scale_colour_gradientn (colours = viridis::viridis(100),limits=c(0, max (umap_p1[[i]]$data[,4])))
  }

pdf (file.path('Plots','label_transfer_celltype_WT_GRANT_umap.pdf'),2,2)
umap_p1
dev.off()


# # Compare WGCNA with SCENIC modules ####
# library (hdWGCNA)
# srt_wgcna = readRDS ('WGCNA_sp_3_ds_1_mch_0.2max_shared_15_minModuleSize_10_metasize_20_genes_3000/srt_wgcna.rds')
# modules = GetModules (srt_wgcna)
# wgcna_mod = split (modules$gene_name, modules$module)
# scenic_mod = lapply (tgenes, function(x) x$gene)

# heat2 = ovmat (scenic_mod)

# heat = ovmat (c(wgcna_mod, scenic_mod), ov_threshold=.3)
# pdf (paste0(projdir, 'Plots/heatmap_wgcna_SCENIC_mod_overlap.pdf'),10,10)
# heat
# heat2
# dev.off()


# Identify DEG in SCENIC ####
motif_window = 'tss500bp'
auc_mtx <- read.csv(paste0('SCENIC/', paste0('vg_',5000,'_mw_',motif_window), '/auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
srt = AddMetaData (srt, metadata = auc_mtx)


# srt_wgcna = readRDS (paste0('WGCNA_sp_3_ds_1_mch_0.2max_shared_15_minModuleSize_10_metasize_20_genes_3000/srt_wgcna.rds'))
# modules = GetModules (srt_wgcna)
# srt = AddMetaData (srt, metadata = srt_wgcna@meta.data[, unique (as.character(modules$module))])

module_thresholds = c(Atf6... = .05)#, green = .0)#, brown=.1, blue = .4, yellow=.2, pink=.4)
mod = 'Atf6...'
metaGroupNames = c('sampleID','celltype_predicted','condition')
  
### Run DEG between condtitions in tan module  + cells
srt@meta.data[,paste0(mod,'_pos')] = ifelse (srt@meta.data[, mod] > module_thresholds[mod], mod, 'rest')
table (srt@meta.data[,paste0(mod,'_pos')], srt$sampleID)

### DEG on genotype ####
  force = FALSE
  do.fgsea = TRUE
  rankby = 'LFC'
  logfcThreshold = 0
  pvalAdjTrheshold = 0.05
  topGenes = 50
  addGene=NULL
  pval_column = 'pval' # specify this if adjusted pvalues are all non-significant
  FeatureSets = list (
  all = NULL)
  #metaGroupName1 = paste0(reductionGraphSnn,'_res.0.8')
  metaGroupName1 = paste0(mod,'_pos')
  metaGroupName2 = 'condition'
  deg2Ident = c('SNI_RGS4_KO','SNI_RGS4_WT')
  source ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scripts/scrna_pipeline/DEG2.R')

write.csv (deg2[deg2$p_val_adj < pvalAdjTrheshold,], file.path(projdir_deg2,'DEG2_genes.csv')) 
deg2 = deg2 [deg2$p_val < pvalAdjTrheshold,]

  
mod_genes = deg2[deg2$cluster == mod,]
mod_genes = mod_genes[order(mod_genes$avg_log2FC),]
#srt_tmp = srt[, srt$RNA_snn_res.0.2 == 5]
srt_tmp = srt[, srt@meta.data[,paste0(mod,'_pos')] == mod]
srt_tmp = ScaleData (srt_tmp, features = mod_genes$gene)

Idents (srt_tmp) = srt_tmp@meta.data[,metaGroupNames[1]]
heat1 = DoHeatmap (srt_tmp, features = mod_genes$gene, raster=F) + ggtitle (paste0('SCENIC ',mod,' module'))

clusters_idx = unique(srt@meta.data[,metaGroupName1])
deg2_sig = deg2[deg2$p_val < pvalAdjTrheshold,]
deg_sum_pos = as.data.frame (table (deg2_sig$cluster[deg2_sig$avg_log2FC > 0]))
rownames(deg_sum_pos) = deg_sum_pos$Var1
deg_sum_pos = deg_sum_pos[clusters_idx,]
rownames (deg_sum_pos) = clusters_idx
deg_sum_neg = as.data.frame (table (deg2_sig$cluster[deg2_sig$avg_log2FC < 0]))
deg_sum_neg$Freq = -1 * deg_sum_neg$Freq
rownames(deg_sum_neg) = deg_sum_neg$Var1
deg_sum_neg = deg_sum_neg[clusters_idx,]
rownames(deg_sum_neg) = clusters_idx
deg_sum = cbind (deg_sum_pos, deg_sum_neg)
colnames(deg_sum)[3:4] = c('Var1_neg','Freq_neg')
deg_sum = deg_sum[mod,]
deg_bar = ggplot(deg_sum, aes(x = rownames(deg_sum))) +
  geom_bar(aes(y = Freq), fill = 'red',color='black', stat = "identity", alpha=0.8) +
  geom_bar(aes(y = Freq_neg), fill = 'blue',color='black', stat = "identity", position = "identity", alpha=0.8) +
  ylim (c(-max(abs(na.omit(c(deg_sum$Freq_neg,deg_sum$Freq)))), 
    max(na.omit(abs(c(deg_sum$Freq_neg,deg_sum$Freq)))))) +
  #scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, by = 10)) +
  labs(title = "DEG per cluster", y = "Value") +
  gtheme

pdf (paste0('Plots/NP_Atf6_number_deg_cluster2.pdf'),4,width=1)
print (deg_bar)
dev.off()
  



# Plot violin plots of selected genes ####
ccomp_df = srt_tmp@meta.data
meta_modules_names = mod_genes$gene
metaGroupNames = c('condgen',meta_modules_names)
ccomp_df = cbind (ccomp_df, as.data.frame (t(srt_tmp@assays$RNA@data[mod_genes$gene,])))
ccomp_df = ccomp_df[,metaGroupNames]
genes_names_starting_with_numbers = meta_modules_names[grep ('^\\d', meta_modules_names)]
meta_modules_names = meta_modules_names[!meta_modules_names %in% genes_names_starting_with_numbers]
  
ccomp_df = ccomp_df[ccomp_df[,'condgen'] %in% c('SNI_WT','SNI_KO'),]
ccomp_df$condgen = factor (ccomp_df$condgen, levels = c('SNI_WT','SNI_KO'))

gene_selection = c('Gabbr2','Mapt','Kcnj3','Kcnab1','Kcnma1','Kcnq5','condgen')
palette_condgen2 = c('SNI KO' = 'grey', 'SNI WT' = 'orange')
ccomp_df = ccomp_df[,gene_selection]
meta_modules_names = gene_selection[-length(gene_selection)]
box_p = list()
stat.store = list()
for (x in meta_modules_names)
  {
  #module_pal = setNames (c('grey55','red'), levels (ccomp_df$condgen))
  max_score = max (ccomp_df[,x]) + 1
  min_score = min (ccomp_df[,x])
  ccomp_df2 = ccomp_df[,c(metaGroupNames[1], x)]
  if (grepl ('\\-',colnames(ccomp_df2)[2])) {
    colnames(ccomp_df2)[2] = sub ('\\-','_',colnames(ccomp_df2)[2])
    x = sub ('\\-','_',colnames(ccomp_df2)[2])
    }
  box_p[[x]] = ggplot (ccomp_df2, aes_string (x= metaGroupNames[1], y= x)) +
      geom_violin (trim=TRUE,aes_string(fill = metaGroupNames[1]),outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
      gtheme +
      scale_fill_manual (values= palette_condgen)
      
  stat.test = box_p[[x]]$data %>%
  #group_by (metaGroupNames[2]) %>%
  rstatix::wilcox_test(reformulate (metaGroupNames[1], x)) %>%
  adjust_pvalue (method = "none") %>%
  add_significance ()
  stat.store[[x]] = stat.test
  stat.test = stat.test %>% add_xy_position (x = metaGroupNames[1], step.increase=0.1)
  box_p[[x]] = box_p[[x]] + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
      bracket.nudge.y = 0.4, hide.ns = F,
      label = "p.adj.signif") + NoLegend()
  }

stat.store = as.data.frame (do.call (rbind, stat.store))
write.csv (stat.store, 'SCENIC_Atf6_NP_selected_genes.csv')

          
dp1 = DimPlot (srt, group.by= paste0(mod,'_pos'), combine=F)
dp1[[1]] = dp1[[1]] + theme_void()
dp2 = fp (srt, mod)
pdf (paste0('Plots/SCENIC_',mod,'_module_pos_umap.pdf'),height=3.5,20)
print(wrap_plots (c(dp1, dp2)))
dev.off()
pdf (paste0('Plots/SCENIC_',mod,'_module_deg_heatmaps_vlnplots.pdf'),height=4,width=5)
#print (heat1)
print (wrap_plots (box_p))
dev.off()




pdf (paste0(projdir,'Plots/Atf6.pdf'))
plotEnrichment(pathways[['GO_synapse_assembly']],
           fgseaRanks[[3]][[1]])
dev.off()

# Pathway analysis of Atf6 SCENIC module ####
scenic_mod = lapply (tgenes, function(x) x$gene)
enricher_universe = rownames (srt)

gmt_annotations = c(
'c2.cp.kegg.v7.1.symbol.gmt',
'c2.cp.reactome.v7.1.symbol.gmt',
'c5.bp.v7.1.symbol.gmt'
#'h.all.v7.1.symbol.gmt'
)

if (!file.exists (paste0(projdir, 'Pathway_Enrichment_clusters.rds')))
  {
  # GSEA analysis on DEG per cluster
  EnrichRResAll = list()
  for (ann in gmt_annotations)
    {
    gmt.file = paste0 ('/ahg/regevdata/projects/ICA_Lung/Bruno/DBs/GSEA_gs/',org,'/',ann)
    pathways = read.gmt (gmt.file)
    #pathways = gmtPathways (gmt.file)
    message (paste('Compute enrichment per cluster using annotation:', ann))
    EnrichRResCluster = list()
    for (i in  seq_along(scenic_mod))
      {
      sig_genes = scenic_mod[[i]]
      if (!all (sig_genes %in% pathways$gene)) next
      egmt <- enricher(sig_genes, TERM2GENE=pathways, universe = enricher_universe)
      EnrichRResCluster[[i]] = egmt@result
      }
    EnrichRResAll[[ann]] = EnrichRResCluster    
    names (EnrichRResAll[[ann]]) = names(scenic_mod)
    }
  saveRDS (EnrichRResAll, paste0(projdir, 'Pathway_Enrichment_clusters.rds')) 
  } else {
  EnrichRResAll = readRDS (paste0(projdir, 'Pathway_Enrichment_clusters.rds'))
  }

lapply (seq_along(EnrichRResAll), function(x) write.csv (EnrichRResAll[[x]], paste0(projdir,'Pathway_enrichment_ann_',x)))  







gene = c('Mrgprd','Calca')

pdf (paste0('Plots/',gene,'_vln_plot.pdf'),width=15)
VlnPlot (srt, features = gene, group.by = 'celltype_predicted2',split.by = 'condition2')
dev.off()


pdf (paste0('Plots/',gene,'_vln_plot.pdf'))
VlnPlot (srt[,srt$condition %in% c('SNI_RGS4_WT','SHAM_RGS4_WT')], features = gene,split.by = 'condition')
VlnPlot (srt[,srt$condition %in% c('SNI_RGS4_WT','SNI_RGS4_KO')], features = gene,split.by = 'condition')
dev.off()



pdf (paste0('Plots/',paste(gene,collapse='_'),'_dotplot.pdf'),height=5)
DotPlot (srt[,srt$condition %in% c('SNI_RGS4_WT','SHAM_RGS4_WT')], features = gene, group.by = 'condition')
DotPlot (srt[,srt$condition %in% c('SNI_RGS4_WT','SNI_RGS4_KO')], features = gene, group.by='condition')
dev.off()


contrast = c('SNI_RGS4_WT','SNI_RGS4_KO')
contrast = c('SNI_RGS4_WT','SHAM_RGS4_WT')
ccomp_df = srt[,srt$condition %in% contrast]@meta.data
genes_mat = as.data.frame (t(srt[,srt$condition %in% contrast]@assays$RNA@data[gene, , drop=F]))
ccomp_df = cbind (ccomp_df, genes_mat)

boxL2=list()

for (contrast in contrasts)
  {
  ccomp_df = srt[,srt$condition %in% contrast]@meta.data
  genes_mat = as.data.frame (t(srt[,srt$condition %in% contrast]@assays$RNA@data[gene, , drop=F]))
  ccomp_df = cbind (ccomp_df, genes_mat)
  
  boxL = list()
  for (i in gene)
    {
    #ccomp_df1 = ccomp_df[ccomp_df[,i] > 0,]  
    box = ggplot (ccomp_df, aes_string (x= 'condition', y= i)) +
            geom_violin (trim=TRUE, aes_string (fill = 'condition')) +
            geom_boxplot (aes_string(fill = 'condition')) +
            #geom_bar (stats='identity') +
            #geom_jitter (color="black", size=0.4, alpha=0.9) +
            gtheme

            stat.test = box$data %>%
            t_test(reformulate ('condition', i)) %>%
            adjust_pvalue (method = "none") %>%
            add_significance ()
            stat.test = stat.test %>% add_xy_position (x = 'condition', step.increase=0.1)
            
    boxL[[i]] = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
            bracket.nudge.y = 0, hide.ns = F,
            label = "p.adj.signif") + NoLegend()
    } 
  boxL2[[paste(contrast, collapse='_')]] = boxL
  }


pdf (paste0('Plots/',paste(gene,collapse='_'),'_', paste(contrast, collapse='_'),'_vlnPlot.pdf'),height=8)
wrap_plots (c(boxL2[[1]],boxL2[[2]]))
dev.off()


# Distributions per sample
gene = c('Mrgprd','Calca')
contrasts = list (c('SNI_RGS4_WT','SNI_RGS4_KO'), c('SNI_RGS4_WT','SHAM_RGS4_WT'))
boxL2=list()
for (contrast in contrasts)
{
ccomp_df = srt[,srt$condition %in% contrast]@meta.data
genes_mat = as.data.frame (t(srt[,srt$condition %in% contrast]@assays$RNA@data[gene, , drop=F]))
ccomp_df = cbind (ccomp_df, genes_mat)

boxL = list()
ccomp_df = aggregate (ccomp_df[,gene],drop=F, by=as.list(ccomp_df[,c('sampleID','condition'),drop=F]), mean)
ccomp_df = na.omit (ccomp_df)
for (i in gene)
  {
  #ccomp_df1 = ccomp_df[ccomp_df[,i] > 0,]  
  box = ggplot (ccomp_df1, aes_string (x= 'condition', y= i)) +
          geom_violin (trim=TRUE, aes_string (fill = 'condition')) +
          geom_boxplot (aes_string(fill = 'condition')) +
          #geom_bar (stats='identity') +
          #geom_jitter (color="black", size=0.4, alpha=0.9) +
          theme_classic() + 
          #scale_fill_manual (values= module_pal) + 
          #ggtitle (SCENIC_mods[x]) + 
          theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
          
          stat.test = box$data %>%
          t_test(reformulate ('condition', i)) %>%
          adjust_pvalue (method = "none") %>%
          add_significance ()
          stat.test = stat.test %>% add_xy_position (x = 'condition', step.increase=0.1)
          
  boxL[[i]] = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
          bracket.nudge.y = 0, hide.ns = F,
          label = "p.adj.signif") + NoLegend()
  } 
boxL2[[paste(contrast, collapse='_')]] = boxL
}

pdf (paste0('Plots/',paste(gene,collapse='_'),'_vlnPlot_per_sample.pdf'),height=8)
wrap_plots (c(boxL2[[1]],boxL2[[2]]))
dev.off()








#### Do the same with green module from WGCNA ####
modulesL = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Serafini_prj/nerve_injury_RGS4_aarthi_analysis/_cellranger_filtered_Filter_0_0_100/no_harmony_nCount_RNA_regressed/NP_subset/no_harmony_nCount_RNA_regressed/WGCNA_sp_3_ds_1_mch_0.2max_shared_15_minModuleSize_10_metasize_20_genes_3000/module_assignments.csv')
wgcna_module_score_ngenes = 100
modulesL = split (modulesL$gene_name, modulesL$module)
modulesL = lapply (modulesL, function(x) head (x, wgcna_module_score_ngenes))
#names (modulesL) = names(split (modules, modules$module))
message ('Run AddModuleScore for each WGCNA module')
srt = ModScoreCor (
      seurat_obj = srt, 
      geneset_list = modulesL, 
      cor_threshold = NULL, 
      pos_threshold = NULL, # threshold for fetal_pval2
      listName = 'wgcna_modules', outdir = paste0(projdir,'Plots/'))

  # Make UMAP of green WGCNA score
  metaGroupNames = c('sampleID','celltype_predicted','condgen')
  SCENIC_mods = 'green'
umap_p1 = FeaturePlot (srt, 
    features = SCENIC_mods, 
    ncol=4, 
    combine = FALSE, 
    pt.size = 1, 
    reduction = reductionName)
for(i in 1:length(umap_p1)) 
  {
  umap_p1[[i]] = umap_p1[[i]] + 
    theme_void() + 
    NoAxes() +
    ggtitle (SCENIC_mods[[i]]) + 
    theme(plot.title = element_text(size = 8)) +
    #scale_colour_gradientn (colours = rev(brewer.pal(n = 11, name = "RdBu")),limits=c(-max (abs(p[[i]]$data[,4])), max (abs(p[[i]]$data[,4]))))
    scale_colour_gradientn (colours = viridis::viridis(100),limits=c(0, max (umap_p1[[i]]$data[,4])))
  }

pdf (file.path('Plots','WGCNA_green_distribution_umap3.pdf'),height=2.7,width = 3)  
umap_p1
dev.off()


### Generate boxplots per meta groups
ccomp_df = srt@meta.data[,SCENIC_mods, drop=FALSE]
ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,metaGroupNames,drop=F]), mean)
#ccomp_df = cbind (ccomp_df, srt@meta.data[,metaGroupNames]) 


sample_rep = srt@meta.data
sample_rep = table (sample_rep[,metaGroupNames[2]], sample_rep[,metaGroupNames[3]])
keep_group = apply (sample_rep , 1, function(y) !any (y < 2))
keep_group = names(keep_group[keep_group])
    
box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods)) +
    #geom_violin (trim=TRUE) +
    geom_boxplot (aes_string(fill = metaGroupNames[3]),outlier.size=2,outlier.alpha = 0.2, notch=FALSE, lwd=.2) +
    #geom_jitter (color="black", size=0.4, alpha=0.9) +
    scale_fill_manual (values= palette_condgen) + 
    #scale_fill_manual (values= module_pal) + 
    ggtitle (SCENIC_mods) +
    ylab ('green module score') + 
   gtheme
  
stat.test = box$data[box$data[,metaGroupNames[2]] %in% keep_group,] %>%
  group_by_at (metaGroupNames[2]) %>%
  t_test(reformulate (metaGroupNames[3], SCENIC_mods)) %>%
  adjust_pvalue (method = "none") %>%
  add_significance ()
  stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
box = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
  bracket.nudge.y = 0, hide.ns = T,
  label = "p.adj.signif") + NoLegend()

pdf (file.path('Plots','WGCNA_green_distribution_boxplots3.pdf'),height=2.7,width = 1.7)  
box_p
dev.off()





library (hdWGCNA)
srt_wgcna = readRDS (paste0('WGCNA_sp_3_ds_1_mch_0.2max_shared_15_minModuleSize_10_metasize_20_genes_3000/srt_wgcna.rds'))
modules = GetModules (srt_wgcna)
srt = AddMetaData (srt, metadata = srt_wgcna@meta.data[, unique (as.character(modules$module))])

module_thresholds = c(green = .0)#, brown=.1, blue = .4, yellow=.2, pink=.4)
mod = 'green'
metaGroupNames = c('sampleID','celltype_predicted','condition')
  
### Run DEG between condtitions in tan module  + cells
srt@meta.data[,paste0(mod,'_pos')] = ifelse (srt@meta.data[, mod] > module_thresholds[mod], mod, 'rest')
table (srt@meta.data[,paste0(mod,'_pos')], srt$sampleID)

### DEG on genotype ####
  force = FALSE
  do.fgsea = TRUE
  rankby = 'LFC'
  logfcThreshold = 0
  pvalAdjTrheshold = 0.05
  topGenes = 50
  addGene=NULL
  pval_column = 'pval' # specify this if adjusted pvalues are all non-significant
  FeatureSets = list (
  all = NULL#,  
  #Il1_pathway = Il1_pathway, 
  #TNF_pathway = TNF_pathway
  )
  #metaGroupName1 = paste0(reductionGraphSnn,'_res.0.8')
  metaGroupName1 = paste0(mod,'_pos')
  metaGroupName2 = 'condition'
  deg2Ident = c('SNI_RGS4_KO','SNI_RGS4_WT')
  source (file.path(scrna_pipeline_dir,'DEG2.R'))

write.csv (deg2[deg2$p_val_adj < pvalAdjTrheshold,], paste0(projdir_deg2,'DEG2_genes.csv'))
deg2 = deg2 [deg2$p_val < pvalAdjTrheshold,]


mod_genes = deg2[deg2$cluster == mod,]
mod_genes = mod_genes[order(mod_genes$avg_log2FC),]
#srt_tmp = srt[, srt$RNA_snn_res.0.2 == 5]
srt_tmp = srt[, srt@meta.data[,paste0(mod,'_pos')] == mod]
srt_tmp = ScaleData (srt_tmp, features = mod_genes$gene)

Idents (srt_tmp) = srt_tmp@meta.data[,metaGroupNames[1]]
heat1 = DoHeatmap (srt_tmp, features = mod_genes$gene, raster=F) + ggtitle (paste0('SCENIC ',mod,' module'))

clusters_idx = unique(srt@meta.data[,metaGroupName1])
deg2_sig = deg2[deg2$p_val < pvalAdjTrheshold,]
deg_sum_pos = as.data.frame (table (deg2_sig$cluster[deg2_sig$avg_log2FC > 0]))
rownames(deg_sum_pos) = deg_sum_pos$Var1
deg_sum_pos = deg_sum_pos[clusters_idx,]
rownames (deg_sum_pos) = clusters_idx
deg_sum_neg = as.data.frame (table (deg2_sig$cluster[deg2_sig$avg_log2FC < 0]))
deg_sum_neg$Freq = -1 * deg_sum_neg$Freq
rownames(deg_sum_neg) = deg_sum_neg$Var1
deg_sum_neg = deg_sum_neg[clusters_idx,]
rownames(deg_sum_neg) = clusters_idx
deg_sum = cbind (deg_sum_pos, deg_sum_neg)
colnames(deg_sum)[3:4] = c('Var1_neg','Freq_neg')
deg_sum = deg_sum['green',]
deg_bar = ggplot(deg_sum, aes(x = rownames(deg_sum))) +
  geom_bar(aes(y = Freq), fill = 'red',color='black', stat = "identity", alpha=0.8) +
  geom_bar(aes(y = Freq_neg), fill = 'blue',color='black', stat = "identity", position = "identity", alpha=0.8) +
  ylim (c(-max(abs(na.omit(c(deg_sum$Freq_neg,deg_sum$Freq)))), 
    max(na.omit(abs(c(deg_sum$Freq_neg,deg_sum$Freq)))))) +
  #scale_y_continuous(limits = c(-60, 60), breaks = seq(-60, 60, by = 10)) +
  labs(title = "DEG per cluster", y = "Value") +
  gtheme

pdf (paste0('Plots/NP_WGCNA_green_module_number_deg_cluster3.pdf'),4,width=1)
print (deg_bar)
dev.off()
  



# Plot violin plots of selected genes ####
ccomp_df = srt_tmp@meta.data
meta_modules_names = mod_genes$gene
metaGroupNames = c('condgen',meta_modules_names)
ccomp_df = cbind (ccomp_df, as.data.frame (t(srt_tmp@assays$RNA@data[mod_genes$gene,])))
ccomp_df = ccomp_df[,metaGroupNames]
genes_names_starting_with_numbers = meta_modules_names[grep ('^\\d', meta_modules_names)]
meta_modules_names = meta_modules_names[!meta_modules_names %in% genes_names_starting_with_numbers]
  
ccomp_df = ccomp_df[ccomp_df[,'condgen'] %in% c('SNI WT','SNI KO'),]
ccomp_df$condgen = factor (ccomp_df$condgen, levels = c('SNI WT','SNI KO'))

gene_selection = c('Syn3','Scn8a','Actg1','Tubb2b', 'condgen')
ccomp_df = ccomp_df[,gene_selection]
meta_modules_names = gene_selection[-length(gene_selection)]
box_p = list()
stat.store = list()
for (x in meta_modules_names)
  {
  #module_pal = setNames (c('grey55','red'), levels (ccomp_df$condgen))
  max_score = max (ccomp_df[,x]) + 1
  min_score = min (ccomp_df[,x])
  ccomp_df2 = ccomp_df[,c(metaGroupNames[1], x)]
  if (grepl ('\\-',colnames(ccomp_df2)[2])) {
    colnames(ccomp_df2)[2] = sub ('\\-','_',colnames(ccomp_df2)[2])
    x = sub ('\\-','_',colnames(ccomp_df2)[2])
    }
  box_p[[x]] = ggplot (ccomp_df2, aes_string (x= metaGroupNames[1], y= x)) +
      geom_violin (trim=TRUE,aes_string(fill = metaGroupNames[1]),outlier.size=2,outlier.alpha = 0.2, notch=FALSE,alpha = 0.7, lwd=.2) +
      gtheme +
      scale_fill_manual (values= palette_condgen)
      
  stat.test = box_p[[x]]$data %>%
  #group_by (metaGroupNames[2]) %>%
  wilcox_test(reformulate (metaGroupNames[1], x)) %>%
  adjust_pvalue (method = "none") %>%
  add_significance ()
  stat.store[[x]] = stat.test
  stat.test = stat.test %>% add_xy_position (x = metaGroupNames[1], step.increase=0.1)
  box_p[[x]] = box_p[[x]] + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
      bracket.nudge.y = 0.4, hide.ns = F,
      label = "p.adj.signif") + NoLegend()
  }

          
dp1 = DimPlot (srt, group.by= paste0(mod,'_pos'), combine=F)
dp1[[1]] = dp1[[1]] + theme_void()
dp2 = fp (srt, mod)
pdf (paste0('Plots/WGCNA_',mod,'_module_pos_umap.pdf'),height=3.5,20)
print(wrap_plots (c(dp1, dp2)))
dev.off()
pdf (paste0('Plots/WGCNA_',mod,'_module_deg_heatmaps_vlnplots.pdf'),height=2,width=7)
#print (heat1)
print (wrap_plots (box_p, ncol=4))
dev.off()


stat.store = as.data.frame (do.call (rbind, stat.store))
write.csv (stat.store, 'WGCNA_green_NP_selected_genes.csv')
