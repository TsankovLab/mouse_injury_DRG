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
  variablefeatures = 'seurat', # options are 'scran', 'seurat','manual'
  nfeat = 2000, # number of variable genes to consider for dimentionality reduction
  sigPCs = 15,
  vars_to_regress='nCount_RNA',
  metaGroupNames = c('sampleID','condition'),
  res = c(0.2, 0.8, 2, 5,10,15) # denovo cluster resolutions 
  )

# Initiate pipeline
force = FALSE # re run pipeline from the beginning no matter if objects are found
subclustername = NULL
projdir_init = paste0 ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Serafini_prj/",proj_name,"_analysis/")
scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'
source (file.path(scrna_pipeline_dir,'master_scrna.R'))

# gtheme = theme(
#       axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),          
#       axis.line = element_blank(),  # Remove axis lines
#       panel.background = element_blank()#,
#     #panel.border = element_blank(),
#     #panel.grid.major = element_blank(),
#     #panel.grid.minor = element_blank()
#   )


if (!file.exists (paste0(projdir,'GSE154659_srt.rds')))
  {
  counts_mat = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Serafini_prj/GSE154659_C57_Raw_counts.rds')
  ref = CreateSeuratObject (counts_mat)
  ref = NormalizeData (object = ref, normalization.method = "LogNormalize", scale.factor = 10000)
  ref = FindVariableFeatures (ref, selection.method = "vst", nfeat = 2000)
  ref = ScaleData (ref, features = VariableFeatures (object=ref))
  ref = RunPCA (ref, features = VariableFeatures (object = ref), npcs = ifelse(ncol(srt) <= 30,ncol(ref)-1,30), ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
  ref = RunUMAP (object = ref, reduction = reductionSave, dims = 1:15)
  ref$celltype = sapply (colnames(counts_mat), function(x) unlist(strsplit(x, '_'))[6])
  saveRDS (ref, paste0(projdir,'GSE154659_srt.rds'))
  } else {
  ref = readRDS (paste0(projdir,'GSE154659_srt.rds'))  
  }
  
  # Label transfer ####
if (!file.exists('prediction_scores_celltype.csv'))
  {
  anchors = FindTransferAnchors (reference = ref, query = srt, dims = 1:15, npcs= 15)
  predictionsMaj = TransferData (anchorset = anchors, refdata = ref$celltype,
  dims = 1:15)
  write.csv (predictionsMaj, 'prediction_scores_celltype.csv')
  } else {
  predictionsMaj = read.csv ('prediction_scores_celltype.csv')
  srt$celltype_predicted = predictionsMaj$predicted.id
  }
  
# Stacked barplot of condition proportions across celltypes ####
srt$condgen2 = sub ('_', ' ', srt$condgen)
metadata = srt@meta.data
metadata = metadata[metadata$celltype_predicted %in% names(table (metadata$celltype_predicted) [table (metadata$celltype_predicted) > 1]),]
bp = cellComp(
  seurat_obj = metadata,
  ptable_factor = 1,
  metaGroups = c('celltype_predicted','condgen2'),
  plot_as = 'bar',
  pal = palette_condgen
  ) + gtheme

pdf (file.path('Plots','cell_proportions_celltypes_conditions.pdf'), width=6, height=4)
bp
dev.off()



srt = RunTSNE (object = srt, reduction = reductionSave, dims = 1:sigPCs)

celltype_ann = c('predicted.subclass','predicted.class','predicted.cluster','predicted.cross_species_cluster')

dp = lapply (celltype_ann, function(x) DimPlot (srt, group.by = x, reduction = reductionName))
dp1 = lapply (celltype_ann, function(x) DimPlot (srt, group.by = x))
dp2 = lapply (celltype_ann, function(x) DimPlot (srt, group.by = x, reduction = 'tsne'))
dp3 = lapply (celltype_ann, function(x) DimPlot (srt, group.by = x, reduction= 'tsne'))
pdf (paste0(projdir,'Plots/celltypes_umaps.pdf'),20,10)
wrap_plots (dp)
wrap_plots (dp2)
dev.off()

srt$celltype_predicted2 = srt$celltype_predicted
srt$celltype_predicted2[srt$RNA_snn_res.0.2 == 7] = 'Macrophage'
srt$celltype_predicted2[srt$RNA_snn_res.0.2 == 5] = 'Endothelial'
srt$celltype_predicted2[srt$RNA_snn_res.0.2 == 2] = 'Fibroblast'
srt$celltype_predicted2[srt$RNA_snn_res.0.2 == 8] = 'Pericyte'
srt$celltype_predicted2[srt$celltype_predicted2 %in% c('NF2','NF3')] = 'NF2/NF3'

srt$Condition = ifelse (grepl ('SHAM',srt$condition), 'Sham','SNI')
srt$condgen = paste0 (srt$Condition, '_', srt$genotype)
srt$condgen = gsub ('_',' ', srt$condgen)
srt$condgen = factor (srt$condgen, levels = c('Sham WT','Sham KO','SNI WT','SNI KO'))

table (srt$condgen)
palette_condgen = setNames (as.character(paletteer::paletteer_d("LaCroixColoR::CranRaspberry")[c(1,2,5,4)]), c('Sham_WT','SNI_WT','Sham_KO','SNI_KO'))
names(palette_condgen) = gsub ('_',' ', names(palette_condgen))
palette_condgen = c('Sham WT' = 'grey','Sham KO' = 'firebrick1', 'SNI WT'= 'darkslategrey', 'SNI KO' = 'violetred4')

store_pal (list (palette_condgen = palette_condgen))
condition_pal = setNames (c('darkseagreen4','orange'), c('Sham','SNI'))
genotype_pal = setNames (c('grey55','red'), c('WT','KO'))


p1 = DimPlot (srt, group.by = 'celltype_predicted2', label=F, label.size = 3) + theme_void() 
p2 = DimPlot (srt, group.by = 'celltype_predicted2', label=T, label.size = 3) + theme_void() 
p2 = DimPlot (srt, group.by = 'Condition', label=F, label.size = 3, cols = condition_pal) + NoLegend()  + theme_void()
p3 =DimPlot (srt, group.by = 'genotype', label=F, label.size = 3, cols = genotype_pal) + NoLegend()  + theme_void()
p4 =DimPlot (srt, group.by = 'condgen', label=F, label.size = 3, cols = palette_condgen) + theme_void() +
scale_fill_manual(values = palette_condgen)

png (paste0(projdir, 'Plots/label_transfer_celltype_umap.png'),width = 3200,1000, res=300)
wrap_plots (p1,p2,p4)
dev.off()

# Rerun UMAP only for WTs ####
#srtWT = RunUMAP()
srtWT = RunUMAP (srt[, srt$condgen %in% c('Sham_WT','SNI_WT')], reduction = reductionSave, dims = 1:sigPCs, reduction.name = reductionName)
p1 = DimPlot (srtWT, group.by = 'celltype_predicted2', label=T, label.size = 3) + theme_void()
pdf (file.path('Plots','label_transfer_celltype_WT_umap.pdf'),4,4)
p1
dev.off()



# QC PLOTS ####
srt$sampleID2 = gsub ('DRG-','',srt$sampleID)
srt$sampleID2 = gsub ('Rgs4','',srt$sampleID2)

srt$nFeature_RNAL = log10 (srt$nFeature_RNA)
srt$nCount_RNAL = log10 (srt$nCount_RNA)
palette_condgen2 = srt$condgen[!duplicated (srt$sampleID2)]
names (palette_condgen2) = srt$sampleID2[!duplicated (srt$sampleID2)]
palette_condgen2 = palette_condgen[palette_condgen2]
names (palette_condgen2) = srt$sampleID2[!duplicated (srt$sampleID2)]
srt$sampleID2 = factor (srt$sampleID2, levels = c('Sham-WT-1','Sham-WT-2','Sham-WT-3','Sham-WT-4','Sham-KO-1','Sham-KO-2','Sham-KO-3','Sham-KO-4','SNI-WT-1','SNI-WT-2','SNI-WT-3','SNI-WT-4','SNI-KO-1','SNI-KO-2','SNI-KO-3','SNI-KO-4'))
vln_p = VlnPlot (srt, features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"), combine=F, group.by = 'sampleID2',pt.size = 0, ncol = 3, col=palette_condgen2)
vln_p = lapply (vln_p, function (x) x + theme(axis.text=element_text(size=9)) + NoLegend())
ccomp_df = as.data.frame (table(srt$sampleID2))
cc_p2 = ggplot (ccomp_df, aes (x= Var1, y= Freq, fill = Var1)) +
        geom_bar (position="stack", stat="identity") +
        ggtitle (paste('Tot cells',ncol(srt))) + theme_classic() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        scale_fill_manual (values = palette_condgen2) +theme(axis.text=element_text(size=9)) + NoLegend()

png (paste0("Plots/QC_nFeat_nCount_m.percent_vlnPlot.png"), 4200, 1000, res=300)
print (cc_p2 | vln_p[[1]] | vln_p[[2]] | vln_p[[3]]) + plot_layout (widths=c(2,2,2,2))
dev.off()

### Check RGS4 expression ####
pdf ('Plots/RGS4_expression.pdf')
VlnPlot (srt, 'Rgs4', group.by = 'condgen', col = palette_condgen) + scale_fill_manual (values = palette_condgen)
dev.off()

# # DEG ####
enricher_universe = 'all'
logfcThreshold = .25
pvalAdjTrheshold = 0.01
metaGroupName = paste0(reductionGraphSnn,'_res.0.8')
top_pathways = 5
force = FALSE
top_genes = 5
source (paste0(scrna_pipeline_dir,'DEG_standard.R'))

### Run WGCNA ####
# Set variables
force=FALSE # force re-running WGCNA 
do.fgsea=TRUE
powerTable=FALSE # plot the powerTable plot. Can take a while to generate
do.plots=TRUE
softPower=6 # Set the softPower
deepSplit=4 # Set this to increase decrease the number of modules idendified. 1-4
mergeCutHeight = 0.20 # Height below which two modules are merged
metacells_k = 20 # number of cells to create pseudobulks
max_shared = 15
metacells_groups = c('RNA_snn_res.0.2') # set metagroup in which to find metacells (usually your clustering / celltypes)
metaGroupNames = c('sampleID','celltype_predicted','condition') # set of metagroups for generating the boxplots
minModuleSize = 10
#genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat = 10000)) # Genes to use to compute WGCNA
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=3000))
enricher_universe = genes.keep # Genes to use as background for pathway enrichments analysis
source (paste0(scrna_pipeline_dir,'scWGCNA.R'))

### Run SCENIC ####
force = FALSE
motif_window = 'tss500bp'#'10kbp'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
source (paste0(scrna_pipeline_dir, 'SCENIC.R'))

# Run SCENIC plots ####
scenic_name = ''
motif_window = 'tss500bp'#'10kbp'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
metaGroupNames = c('sampleID','celltype_predicted','condition')
source (file.path(scrna_pipeline_dir, 'SCENIC_plots.R'))

metaGroupNames = c('sampleID','celltype_predicted','condgen')
ccomp_df = srt@meta.data
scenic_tf = c('Ets1','Fos','Gtf2ird1','Hdx','Mafa','Ppara','Prrxl1','Stat1')
scenic_tf = colnames (srt@meta.data)[grep ('SCENIC_',colnames (srt@meta.data))]
ccomp_df$condgen = gsub ('_',' ',ccomp_df$condgen)
ccomp_df$condgen = factor (ccomp_df$condgen, levels = c('Sham WT','Sham KO','SNI WT','SNI KO'))
SCENIC_mods2 = scenic_tf
box_p = lapply (seq(SCENIC_mods2), function(x) 
    {
    sample_rep = srt@meta.data
    sample_rep = table (sample_rep[,metaGroupNames[2]], sample_rep[,metaGroupNames[3]])
    keep_group = apply (sample_rep , 1, function(y) !any (y < 2))
    keep_group = names(keep_group[keep_group])
      if (length(keep_group) > 0 & metaGroupNames[1] != metaGroupNames[3])
      {
      box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods2[x])) +
          #geom_violin (trim=TRUE) +
          geom_boxplot (aes_string(fill = metaGroupNames[3])) +
          #geom_jitter (color="black", size=0.4, alpha=0.9) +
          theme_classic() + 
          scale_fill_manual (values= palette_condgen) + 
          ggtitle (SCENIC_mods2[x]) + 
          theme (axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          labs(x = "cond - geno", y = "module score")
  
      stat.test = box$data[box$data[,metaGroupNames[2]] %in% keep_group,] %>%
        group_by_at (metaGroupNames[2]) %>%
        t_test(reformulate (metaGroupNames[3], SCENIC_mods2[x])) %>%
        adjust_pvalue (method = "none") %>%
        add_significance ()
        stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
        box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
        bracket.nudge.y = 0, hide.ns = T,
        label = "p.adj.signif") + NoLegend()
      } else {
      ccomp_df = srt@meta.data[,c(SCENIC_mods2, metaGroupNames), drop=FALSE]
      #ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)    
      box = ggplot (ccomp_df, aes_string (x= metaGroupNames[3], y= SCENIC_mods2[x])) +
        geom_violin (trim=TRUE, aes_string (fill = metaGroupNames[3])) +
        geom_boxplot (aes_string(fill = metaGroupNames[3])) +
        #geom_bar (stats='identity') +
        #geom_jitter (color="black", size=0.4, alpha=0.9) +
        theme_classic() + 
        scale_fill_manual (values= palette_condgen) + 
        ggtitle (SCENIC_mods2[x]) + 
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(x = "cond - geno", y = "module score")
        
        stat.test = box$data %>%
        group_by_at (metaGroupNames[2]) %>%
        wilcox_test(reformulate (metaGroupNames[3], SCENIC_mods2[x])) %>%
        adjust_pvalue (method = "none") %>%
        add_significance ()
        stat.test = stat.test %>% add_xy_position (x = metaGroupNames[3], step.increase=0.1)
        
        box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
        bracket.nudge.y = 0, hide.ns = T,
        label = "p.adj.signif") + NoLegend()
      }})

  if (length(unique(ccomp_df[,metaGroupNames[2]])) > 1) box_p = lapply (seq_along(SCENIC_mods2), function(x) box_p[[x]] + facet_wrap (as.formula(paste("~", metaGroupNames[2]))) + theme_classic() + 
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + NoLegend()) 
      
    png (file.path (projdir_SC_run, 'Plots',paste0('SCENIC_module_scores_',paste(metaGroupNames,collapse='_'),'_umap.png')), width = 5000, height = 6000, pointsize=10, res = 300, type="cairo")
    print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))) / 
      wrap_plots (box_p, ncol= ifelse (length(box_p) > 5,ceiling(length(box_p)/2),length(box_p))) + plot_layout (heights=c(1,2)))
    dev.off()
      
    pdf (file.path (projdir_SC_run, 'Plots',paste0('SCENIC_module_scores_',paste(metaGroupNames,collapse='_'),'_umap.pdf')), width = 80+length(umap_p1), height = 28)
    print (wrap_plots (umap_p1, ncol = ifelse (length(umap_p1) > 8,ceiling(length(umap_p1)/2),length(umap_p1))) / 
      wrap_plots (box_p, ncol= ifelse (length(box_p) > 5,ceiling(length(box_p)/2),length(box_p))) + plot_layout (heights=c(1,2)))
    dev.off()

    png (file.path (projdir_SC_run, 'Plots',paste0('SCENIC_module_scores_',paste(metaGroupNames,collapse='_'),'_boxplots.png')), width = 8000, height = 6000, pointsize=10, res = 300, type="cairo")
    wrap_plots (box_p, ncol= ifelse (length(box_p) > 5,ceiling(length(box_p)/2),length(box_p))) + plot_layout (heights=c(1,1))
    dev.off()
          
  

# Subset to NP ####
force=FALSE
subclustername = 'NP'
#metaGroupName='RNA_snn_res.0.8'
#metaGroupSelection=c(2,3,4,6,12,15,16,18,20,21)
metaGroupName=paste0('celltype_predicted')
metaGroupSelection=c('NP')
exclude=FALSE
source (paste0(scrna_pipeline_dir,'subcluster.R'))


# Subset to NP ####
force=FALSE
subclustername = 'Schwann'
#metaGroupName='RNA_snn_res.0.8'
#metaGroupSelection=c(2,3,4,6,12,15,16,18,20,21)
metaGroupName=paste0('celltype_predicted')
metaGroupSelection=c('Schwann')
exclude=FALSE
source (paste0(scrna_pipeline_dir,'subcluster.R'))



### DEG on genotype ####
srt$condition2 = ifelse (srt$condition %in% c('SNI_RGS4_KO','SNI_RGS4_WT'), 'SNI','SHAM')
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
metaGroupName1 = 'celltype_predicted2'
metaGroupName2 = 'condgen'
deg2Ident = c('SNI','SHAM')
deg2Ident = c('Sham_WT','Sham_KO')
source (file.path (scrna_pipeline_dir,'DEG2.R'))


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


boxL = list()
for (i in gene)
  {
  ccomp_df1 = ccomp_df[ccomp_df[,i] > 0,]  
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
          wilcox_test(reformulate ('condition', i)) %>%
          adjust_pvalue (method = "none") %>%
          add_significance ()
          stat.test = stat.test %>% add_xy_position (x = 'condition', step.increase=0.1)
          
  boxL[[i]] = box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
          bracket.nudge.y = 0, hide.ns = F,
          label = "p.adj.signif") + NoLegend()
  }


pdf (paste0('Plots/',paste(gene,collapse='_'),'_', paste(contrast, collapse='_'),'_vlnPlot.pdf'),height=5)
wrap_plots (boxL)
dev.off()


contrast = c('SNI_RGS4_WT','SNI_RGS4_KO')
contrast = c('SNI_RGS4_WT','SHAM_RGS4_WT')
ccomp_df = srt[,srt$condition %in% contrast]@meta.data
genes_mat = as.data.frame (t(srt[,srt$condition %in% contrast]@assays$RNA@data[gene, , drop=F]))
ccomp_df = cbind (ccomp_df, genes_mat)

boxL = list()
ccomp_df = aggregate (ccomp_df[,gene],drop=F, by=as.list(ccomp_df[,c('sampleID','condition'),drop=F]), mean)
ccomp_df = na.omit (ccomp_df)
for (i in gene)
  {
  ccomp_df1 = ccomp_df[ccomp_df[,i] > 0,]  
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


pdf (paste0('Plots/',paste(gene,collapse='_'),'_', paste(contrast, collapse='_'),'_vlnPlot_per_sample.pdf'),height=5)
wrap_plots (boxL)
dev.off()



### DEG on genotype ####
srt$project = 'RGS4'
force = TRUE
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
metaGroupName1 = 'project'
metaGroupName2 = 'condition'
deg2Ident = c('SNI_RGS4_KO','SNI_RGS4_WT')
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/DEG2.R')




### vlnplots ####
gene = 'Phf24'
pdf (paste0('Plots/',gene,'_vlnplot.pdf'))
VlnPlot (srt, gene, group.by = 'condition')
dev.off()









  # Distributions per sample
gene = c('Phf24')
gene = c('Rgs4')

# Plot violin plots of selected genes ####
ccomp_df = srt@meta.data
meta_modules_names = gene
metaGroupNames = c('condgen','celltype_predicted2', meta_modules_names,'sampleID2')
ccomp_df = cbind (ccomp_df, as.data.frame (t(srt@assays$RNA@data[gene,, drop=F])))
ccomp_df = ccomp_df[,metaGroupNames]

noise <- rnorm(n = length(x = ccomp_df[, gene])) / 100000  
if (all(ccomp_df[, gene] == ccomp_df[, gene][1])) { 
  warning(paste0("All cells have the same value of ", gene, ".")) 
} else{ 
  ccomp_df[, gene] <- ccomp_df[, gene] + noise
}
#ccomp_df = ccomp_df[ccomp_df[,'condgen'] %in% c('SNI WT','SNI KO'),]
#ccomp_df$condgen = factor (ccomp_df$condgen, levels = c('SNI WT','SNI KO'))
box_p = ggplot (ccomp_df, aes_string (x= 'condgen', y= gene)) +
#geom_jitter() + 
    geom_violin (trim=TRUE,aes_string (fill = 'condgen'),
      scale = 'width',
      )  + 
    gtheme +
    facet_wrap (~celltype_predicted2, ncol=7) + 
    scale_fill_manual (values= palette_condgen)
      
  stat.test = box_p$data %>%
  group_by (celltype_predicted2) %>%
  wilcox_test(reformulate (metaGroupNames[1], gene)) %>%
  adjust_pvalue (method = "none") %>%
  add_significance ()
  stat.test = stat.test %>% add_xy_position (x = 'condgen', step.increase=0.1)
  box_p = box_p + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
      bracket.nudge.y = 0.4, hide.ns = T,
      label = "p.adj.signif") + NoLegend()

write.csv (as.data.frame (stat.test[,c(2,3,4,5,6,7,8,9)]), paste0(gene,'_stats.csv'))

pdf (paste0('Plots/',gene,'_vlnplots.pdf'),height=5,width=10)
#print (heat1)
print (wrap_plots (box_p))
dev.off()

palette_gene_expression = viridis::inferno (100)
srt$condgen2 = factor (srt$condgen, levels = rev (levels (srt$condgen)))
dp = geneDot (
  seurat_obj = srt,
  #gene = top_tfs2, 
  gene = gene,
  x = 'condgen2', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  y = 'condgen2',
  z = NULL, 
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=F,
  x_name ='',
  y_name = '',
  plotcol = palette_gene_expression) +
    theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) 
dp + scale_fill_gradientn (colors = palette_gene_expression)

pdf (paste0('Plots/',gene,'_expression_leg.pdf'), height=4, width=3)
dp
dev.off()

# per sample

noise <- rnorm(n = length(x = ccomp_df[, gene])) / 100000  
if (all(ccomp_df[, gene] == ccomp_df[, gene][1])) { 
  warning(paste0("All cells have the same value of ", gene, ".")) 
} else{ 
  ccomp_df[, gene] <- ccomp_df[, gene] + noise
}
box_p = ggplot (ccomp_df, aes_string (x= 'sampleID2', y= gene)) +
#geom_jitter() + 
    # geom_violin (trim=TRUE,aes_string (fill = 'condgen'),
    #   scale = 'width',
    #   )  + 
geom_boxplot (aes_string(fill='condgen'),color = 'grey22', width=.5, alpha = 1, lwd=.2, outlier.shape = NA) +
    gtheme +
    #facet_wrap (~celltype_predicted2, ncol=7) + 
    scale_fill_manual (values= palette_condgen)
      
write.csv (as.data.frame (stat.test[,c(2,3,4,5,6,7,8,9)]), paste0(gene,'_stats.csv'))

pdf (paste0('Plots/',gene,'_per_sample_vlnplots.pdf'),height=4,width=5)
#print (heat1)
print (wrap_plots (box_p))
dev.off()

palette_gene_expression = viridis::viridis (100)
srt$sampleID2 = factor (srt$sampleID2, levels = rev (levels(srt$sampleID2)))
gene = c('Rgs4')
dp = geneDot (
  seurat_obj = srt,
  #gene = top_tfs2, 
  gene = gene,
  x = 'sampleID2', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  y = NULL, 
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=F,
  x_name ='',
  y_name = '',
  swap_axes = F,
  plotcol = palette_gene_expression) +
    theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) 
#dp + scale_fill_gradientn (colors = palette_gene_expression)

pdf (paste0('Plots/',gene,'_per_sample_dotplot.pdf'),height=2,width=5)
dp
dev.off()








# check only in NP cells
gene = 'Rgs4'
ccomp_df = srt@meta.data[srt$celltype_predicted2 == 'NP',]
meta_modules_names = gene
metaGroupNames = c('condgen',meta_modules_names,'sampleID2')
ccomp_df = cbind (ccomp_df, as.data.frame (t(srt[,srt$celltype_predicted2 == 'NP']@assays$RNA@data[gene,, drop=F])))
ccomp_df = ccomp_df[,metaGroupNames]

noise <- rnorm(n = length(x = ccomp_df[, gene])) / 100000  
if (all(ccomp_df[, gene] == ccomp_df[, gene][1])) { 
  warning(paste0("All cells have the same value of ", gene, ".")) 
} else{ 
  ccomp_df[, gene] <- ccomp_df[, gene] + noise
}
box_p = ggplot (ccomp_df, aes_string (x= 'sampleID2', y= gene)) +
#geom_jitter() + 
    # geom_violin (trim=TRUE,aes_string (fill = 'condgen'),
    #   scale = 'width',
    #   )  + 
geom_boxplot (aes_string(fill='condgen'),color = 'grey22', width=.5, alpha = 1, lwd=.2, outlier.shape = NA) +
    gtheme +
    #facet_wrap (~celltype_predicted2, ncol=7) + 
    scale_fill_manual (values= palette_condgen)
      
write.csv (as.data.frame (stat.test[,c(2,3,4,5,6,7,8,9)]), paste0(gene,'_stats.csv'))

pdf (paste0('Plots/',gene,'_per_sample_NP_boxplots.pdf'),height=4,width=5)
#print (heat1)
print (wrap_plots (box_p))
dev.off()

box_p = ggplot (ccomp_df, aes_string (x= 'sampleID2', y= gene)) +
# geom_jitter() + 
    geom_violin (trim=TRUE,aes_string (fill = 'condgen'),
      scale = 'width',
      )  + 
# geom_boxplot (aes_string(fill='condgen'),color = 'grey22', width=.5, alpha = 1, lwd=.2, outlier.shape = NA) +
    gtheme +
    #facet_wrap (~celltype_predicted2, ncol=7) + 
    scale_fill_manual (values= palette_condgen)
      
write.csv (as.data.frame (stat.test[,c(2,3,4,5,6,7,8,9)]), paste0(gene,'_stats.csv'))

pdf (paste0('Plots/',gene,'_per_sample_NP_vlnplots.pdf'),height=4,width=5)
#print (heat1)
print (wrap_plots (box_p))
dev.off()

palette_gene_expression = viridis::viridis (100)
srt$sampleID2 = factor (srt$sampleID2, levels = rev (levels(srt$sampleID2)))
gene = c('Rgs4')
dp = geneDot (
  seurat_obj = srt[,srt$celltype_predicted2 == 'NP'],
  #gene = top_tfs2, 
  gene = gene,
  x = 'sampleID2', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  y = NULL, 
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data=TRUE,
  x_name ='',
  y_name = '',
  swap_axes = F,
  plotcol = palette_gene_expression) +
    theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) 
#dp + scale_fill_gradientn (colors = palette_gene_expression)

pdf (paste0('Plots/',gene,'_per_sample_NP_dotplot.pdf'),height=3,width=5)
dp
dev.off()



gene = 'Hdac1'
srt$condition2 = gsub ('RGS4_','',srt$condition)
pdf (file.path('Plots',paste0(paste(gene,collapse='_'),'_dotplot.pdf')),height=5)
DotPlot (srt[,srt$condition %in% c('SNI_RGS4_WT','SHAM_RGS4_WT')], features = gene, group.by = 'celltype_predicted', split.by = 'condition2')

# dp = geneDot (
#   seurat_obj = srt,
#   #gene = top_tfs2, 
#   gene = gene,
#   x = 'condgen2', # if multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
#   #y = 'condgen2',
#   #z = NULL, 
#   min_expression = 0,
#   facet_ncol = 5,
#   lim_expression = NULL,
#   scale.data=F,
#   x_name ='',
#   y_name = '',
#   plotcol = palette_gene_expression) +
#     theme (axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) 
# dp + scale_fill_gradientn (colors = palette_gene_expression)

#DotPlot (srt[,srt$condition %in% c('SNI_RGS4_WT','SNI_RGS4_KO')], features = 'Hdac1', group.by='condition')
dev.off()












#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(10:30)
k_selections = c(10:30)
cores= 100

# Run cNMF only in samples matching scATAC-seq samples 
  
cnmf_name = paste0('main')
#cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
#dir.create (file.path(cnmf_out,'Plots'), recursive=T)

  
### RUN consensus NMF ####
force=F
source (file.path (scrna_pipeline_dir,'cnmf_prepare_inputs.R'))     


### Import and format spectra files ####
cnmf_name = 'main'
k_selection=30
nfeat=2000
source (file.path (scrna_pipeline_dir,'cnmf_format_spectra_files.R')) 
cnmf_spectra_unique

# Check module including Lrfn5
mod = 'cNMF18'
gene = 'Lrfn5'
module = cnmf_spectra_unique[[mod]]
module_mat = srt@assays$RNA@data[rownames(srt) %in% module, ]
module_mat = t(scale(t(module_mat)))

palette_expression = colorRamp2(c(-2,0,2), c('blue','white','red'))
ha = HeatmapAnnotation (ct = srt$celltype_predicted2)

pdf()
hm = draw (Heatmap (module_mat, top_annotation = ha,
  row_names_gp = gpar(fontsize = 6), row_km = 5,
  column_names_gp = gpar(fontsize = 0), clustering_distance_rows = 'pearson',
  col = palette_expression))
dev.off()
pdf (file.path ('Plots',paste0(mod,'_heatmap.pdf')), height=5)
hm
#VlnPlot (srt, features = 'Bmp6', group.by = 'Hughes_celltype')
dev.off()

gene_groups = sapply (row_order (hm), function(x) rownames (module_mat)[x])
gene_coexp = gene_groups[unlist(sapply(gene_groups, function(x) gene %in% x))]

genes = gene_coexp[[1]]

names(palette_condgen) = gsub(' ','_',names(palette_condgen))
metagroupnames = c('sampleID','celltype_predicted2','condgen')
ccomp_df = srt@meta.data[,metagroupnames]
exp_df = as.data.frame (t(srt@assays$RNA@data[genes,,drop=F]))
ccomp_df = aggregate (exp_df, 
  by = as.list(ccomp_df), 
  mean)  
head (ccomp_df)
#ccomp_df = ccomp_df[ccomp_df$condgen %in% c('Sham_WT','SNI_WT'),]
box = list()  
for (gene in genes)
    {
    box[[gene]] = ggplot (ccomp_df, aes_string (x= metagroupnames[2], y= gene)) +
    geom_boxplot (aes_string (fill = metagroupnames[3]),
    linewidth = .2,
    width=.8,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7
     ) +
    # geom_violin (aes_string (fill = metagroupnames[3]),
    #   trim=TRUE,size=2,
    #   width=1,
    #   scale='width',
    #   linewidth = .2, alpha=0.7) +
    theme_classic() +
    scale_fill_manual (values= palette_condgen) +
    ggtitle (gene) +
    theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
    #NoLegend()
    
    stat.test = box[[gene]]$data %>%
  group_by_at(metagroupnames[2])  %>% #filter(
    #sum(condgen == "SNI_WT") >= 2,
    #sum(condgen == "Sham_WT") >= 2) %>%

  t_test(
    formula = reformulate(metagroupnames[3], response = gene),
    ref.group = "Sham_WT"
  ) %>% 
  adjust_pvalue (method = "none") %>%
  add_significance () %>%
  #ungroup() %>%
  add_xy_position(
    data = box[[gene]]$data,
    formula = reformulate(metagroupnames[3], response = gene),
    x = metagroupnames[2],     # <-- x = RNA_snn_res.0.8 (the cluster)
    #group.by = metagroupnames[2],
    step.increase = 0.01
  )

    box[[gene]] = box[[gene]] + stat_pvalue_manual (
      stat.test, remove.bracket=FALSE,
    bracket.nudge.y = 0, hide.ns = T,
    x = metagroupnames[2],
    label = "p.adj.signif") 
    }

pdf (file.path ('Plots',paste0(gene,'_to_test_boxplots.pdf')), width=20, height=8)
wrap_plots(box)
dev.off()

pdf (file.path ('Plots','clustering_dimplot.pdf'), 4,4)
DimPlot (srt, group.by = 'celltype_predicted2')
dev.off()




cnmf_spectra_unique_l = lapply (cnmf_spectra_unique, function(x) head(x,50))

srt = ModScoreCor (
      seurat_obj = srt, 
      geneset_list = cnmf_spectra_unique_l, 
      cor_threshold = NULL, 
      pos_threshold = NULL, # threshold for fetal_pval2
      listName = 'module', outdir = paste0(projdir,'Plots/'))

meta_modules_names = names(cnmf_spectra_unique_l)
metaGroupNames = c('sampleID','condition','celltype2')
  umap_df = data.frame (srt[[reductionName]]@cell.embeddings, srt@meta.data[,c(meta_modules_names,metaGroupNames)])
  umap_p1 = lapply (meta_modules_names, function(x) ggplot(data = umap_df) + 
  geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
  scale_colour_gradientn (colours = rev(brewer.pal (n = 11, name = "RdBu")),limits=c(-max (abs(umap_df[,x])), max (abs(umap_df[,x])))) +
  ggtitle (x) + 
  #facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + 
  theme_classic() +
  theme_void())

pdf (file.path ('Plots','cNMF_modules_umaps.pdf'), width=12)
wrap_plots (umap_p1)
dev.off()

ccomp_df = srt@meta.data
ccomp_df = aggregate (ccomp_df[,meta_modules_names], 
  by = list(
  sample = ccomp_df$sampleID, 
  condition = ccomp_df$condition, 
  celltype2 = ccomp_df$celltype2), mean)

#ccomp_df$treatment = factor (ccomp_df$treatment, levels = c('pre','post'))
ccomp_df = gather (ccomp_df, module, score, which (colnames(ccomp_df) %in% meta_modules_names))
ccomp_df$module = factor (ccomp_df$module, levels = meta_modules_names)
bp = ggplot (ccomp_df, aes (x= celltype2, y= score)) +
      #geom_violin (trim=TRUE, aes (fill = treatment), alpha=.6) +
      #geom_violin (aes_string(fill = metaGroupNames[3])) +
      #geom_point (aes (x = condition, y = score), position='identity', alpha=.7, color="blue", size=1.2) +
      geom_boxplot(width=0.5, aes (fill = condition, color = condition), alpha=.6) +
      #scale_fill_manual (values= palette_treatment) + 
      #scale_color_manual (values= palette_treatment) + 
      #geom_line (data = ccomp_df[ccomp_df$paired2 == 'paired',], aes(x = condition, y = score, group = patient), color='blue',linewidth=.2, alpha=.7) +
      gtheme + 
      NoLegend() + 
      facet_wrap (~module, scales = 'free_y',strip.position = "left") +
      theme(strip.background = element_rect(fill = "white", color = "black")) 

stat.test2 = bp$data %>% group_by (module) %>%
  rstatix::t_test(reformulate('condition', 'score'), paired=F) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

stat.test2 = stat.test2 %>% add_xy_position (x = "condition", step.increase=0.005)
bp = bp + stat_pvalue_manual (stat.test2, remove.bracket=FALSE,
   bracket.nudge.y = 0, hide.ns = TRUE, position = position_nudge (y = -0.01),
    label = "p.adj.signif") 

pdf (file.path ('Plots','cNMF_modules_boxplots.pdf'),width=18, height=5)
wrap_plots (bp)
dev.off()


#ccomp_df$treatment = factor (ccomp_df$treatment, levels = c('pre','post'))


ccomp_df = srt@meta.data
ccomp_df$neurons = 'Neurons'
ccomp_df = aggregate (ccomp_df[,meta_modules_names], 
  by = list(
  sample = ccomp_df$sampleID, 
  condition = ccomp_df$condition, 
  celltype2 = ccomp_df$neurons), mean)

ccomp_df = gather (ccomp_df, module, score, which (colnames(ccomp_df) %in% meta_modules_names))
ccomp_df$module = factor (ccomp_df$module, levels = meta_modules_names)
bp = ggplot (ccomp_df, aes (x= celltype2, y= score)) +
      #geom_violin (trim=TRUE, aes (fill = treatment), alpha=.6) +
      #geom_violin (aes_string(fill = metaGroupNames[3])) +
      #geom_point (aes (x = condition, y = score), position='identity', alpha=.7, color="blue", size=1.2) +
      geom_boxplot(width=0.5, aes (fill = condition, color = condition), alpha=.6) +
      #scale_fill_manual (values= palette_treatment) + 
      #scale_color_manual (values= palette_treatment) + 
      #geom_line (data = ccomp_df[ccomp_df$paired2 == 'paired',], aes(x = condition, y = score, group = patient), color='blue',linewidth=.2, alpha=.7) +
      gtheme + 
      #NoLegend() + 
      facet_wrap (~module, scales = 'free_y',strip.position = "left") +
      theme(strip.background = element_rect(fill = "white", color = "black")) 

stat.test2 = bp$data %>% group_by (module) %>% 
  rstatix::t_test(reformulate('condition', 'score'), paired=F) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

stat.test2 = stat.test2 %>% add_xy_position (x = "condition", step.increase=0.005)
bp = bp + stat_pvalue_manual (stat.test2, remove.bracket=FALSE,
   bracket.nudge.y = 0, hide.ns = TRUE, position = position_nudge (y = -0.01),
    label = "p.adj.signif") 

pdf (file.path ('Plots','cNMF_modules2_boxplots.pdf'),width=18, height=5)
wrap_plots (bp)
dev.off()


# Check in which module is Mef2c and Itgb5
sapply (cnmf_spectra_unique, function(x) 'Mef2c' %in% x)
sapply (cnmf_spectra_unique, function(x) 'Itgb5' %in% x)






### Run cellphoneDB ####
force=T
metaGroupName1 = 'sampleID'
metaGroupName2 = 'celltype_predicted2'
#metaGroupName2 = 'Hughes_celltype'
metaGroupName3 = 'condgen'
source (file.path (scrna_pipeline_dir, 'cellphoneDB_run.R'))

### generate plots from CPDB output ####
compGroups = c('SNI_WT','Sham_WT')
compGroups = c('SNI_WT','Sham_WT')
compGroups = c('SNI_KO','SNI_WT')
compGroups = c('Sham_KO','Sham_WT')

lr_filter = NULL
#celltype_filter = c('Bcells','DCs','Hypoxic_mac','Macrophages_1','Macrophages_2','Mast_cells','T_cells')
force=F
celltype_filter = c('T_cells')
celltype_filter = c('Myeloid')
celltype_filter = NULL
top_lr_filter = NULL
minCells = 20
pValThreshold = -log10(0.05)
source (file.path (scrna_pipeline_dir, 'cellphoneDB_downstream.R'))


meta_modules_names = c('Bmp6','Acvr1')
metaGroupName = 'celltype_predicted2'
ccomp_df = srt@meta.data
metaGroupNames = c('condgen','celltype_predicted2',meta_modules_names)
ccomp_df = cbind (ccomp_df, as.data.frame (t(srt@assays$RNA@data[meta_modules_names,,drop=F])))
ccomp_df = ccomp_df[,metaGroupNames]
meta_modules_names = meta_modules_names[!meta_modules_names %in% '5730522E02Rik']
box = lapply (meta_modules_names, function(x) 
  {
    box = ggplot (ccomp_df, aes_string (x= metaGroupNames[2], y= x)) +
        geom_violin (trim=TRUE, aes_string (fill = metaGroupNames[1])) +
        #geom_bar (stats='identity') +
        #geom_jitter (color="black", size=0.4, alpha=0.9) +
        theme_classic() + 
        #scale_fill_manual (values= module_pal) + 
        ggtitle (x) + 
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        stat.test = box$data %>%
        dplyr::group_by_at (metaGroupNames[2]) %>%
        rstatix::wilcox_test (reformulate ('condgen', x)) %>%
        adjust_pvalue (method = "none") %>%
        add_significance ()
        stat.test = stat.test %>% add_xy_position (x = metaGroupNames[1], step.increase=0.01)
        box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
        bracket.nudge.y = 0, hide.ns = F,
        label = "p.adj.signif") + NoLegend()
      })
          
dp1 = DimPlot (srt, group.by= metaGroupName, combine=F) 
dp2 = fp (srt_wgcna, meta_modules_names) 
v1 = VlnPlot(srt, 'Mef2c', group.by = 'celltype2', split.by = 'condgen')
pdf (file.path(projdir,'Plots',paste0(meta_modules_names,'_deg_heatmaps_vlnplots.pdf')),height=6,8)
wrap_plots (c(dp1, dp2))
wrap_plots (box)
v1
dev.off()









#### plot genes of interest
genes = c ('Ntng2','Snap25','Sirt1','Slc22a6','Plk2','Terf2')
genes = c ('Th','Ddc','Dbh','Comt','Maoa','Maob','Slc6a3','Drd1','Drd2','Drd3',
  'Drd4','Drd5','Nln','Nts','Ntsr1','Ntsr2')
genes = genes[genes %in% rownames(srt)]

names(palette_condgen) = gsub(' ','_',names(palette_condgen))
metagroupnames = c('sampleID','celltype_predicted2','condgen')
ccomp_df = srt@meta.data[,metagroupnames]
exp_df = as.data.frame (t(srt@assays$RNA@data[genes,]))
ccomp_df = aggregate (exp_df, 
  by = as.list(ccomp_df), 
  mean)  
head (ccomp_df)
ccomp_df = ccomp_df[ccomp_df$condgen %in% c('Sham_WT','SNI_WT'),]
box = list()  
for (gene in genes)
    {
    box[[gene]] = ggplot (ccomp_df, aes_string (x= metagroupnames[2], y= gene)) +
    geom_boxplot (aes_string (fill = metagroupnames[3]),
    linewidth = .2,
    width=.8,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7
     ) +
    # geom_violin (aes_string (fill = metagroupnames[3]),
    #   trim=TRUE,size=2,
    #   width=1,
    #   scale='width',
    #   linewidth = .2, alpha=0.7) +
    theme_classic() +
    scale_fill_manual (values= palette_condgen) +
    ggtitle (gene) +
    theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    NoLegend()
    
    stat.test = box[[gene]]$data %>%
  group_by_at(metagroupnames[2])  %>% filter(
    sum(condgen == "SNI_WT") >= 2,
    sum(condgen == "Sham_WT") >= 2
  ) %>%

  t_test(
    formula = reformulate(metagroupnames[3], response = gene),
    ref.group = "Sham_WT"
  ) %>% 
  adjust_pvalue (method = "none") %>%
  add_significance () %>%
  #ungroup() %>%
  add_xy_position(
    data = box[[gene]]$data,
    formula = reformulate(metagroupnames[3], response = gene),
    x = metagroupnames[2],     # <-- x = RNA_snn_res.0.8 (the cluster)
    #group.by = metagroupnames[2],
    step.increase = 0.01
  )

    box[[gene]] = box[[gene]] + stat_pvalue_manual (
      stat.test, remove.bracket=FALSE,
    bracket.nudge.y = 0, hide.ns = T,
    x = metagroupnames[2],
    label = "p.adj.signif") 
    }

pdf (file.path ('Plots','genes_to_test_boxplots.pdf'), width=12, height=12)
wrap_plots(box)
dev.off()

pdf (file.path ('Plots','clustering_dimplot.pdf'))
DimPlot (srt, group.by = 'celltype_predicted2')
dev.off()


genes = c ('Lrfn5')
genes = genes[genes %in% rownames(srt)]

names(palette_condgen) = gsub(' ','_',names(palette_condgen))
metagroupnames = c('sampleID','celltype_predicted2','condgen')
ccomp_df = srt@meta.data[,metagroupnames]
exp_df = as.data.frame (t(srt@assays$RNA@data[genes,,drop=F]))
ccomp_df = aggregate (exp_df, 
  by = as.list(ccomp_df), 
  mean)  
head (ccomp_df)
#ccomp_df = ccomp_df[ccomp_df$condgen %in% c('Sham_WT','SNI_WT'),]
box = list()  
for (gene in genes)
    {
    box[[gene]] = ggplot (ccomp_df, aes_string (x= metagroupnames[2], y= gene)) +
    geom_boxplot (aes_string (fill = metagroupnames[3]),
    linewidth = .2,
    width=.8,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7
     ) +
    # geom_violin (aes_string (fill = metagroupnames[3]),
    #   trim=TRUE,size=2,
    #   width=1,
    #   scale='width',
    #   linewidth = .2, alpha=0.7) +
    theme_classic() +
    scale_fill_manual (values= palette_condgen) +
    ggtitle (gene) +
    theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
    #NoLegend()
    
    stat.test = box[[gene]]$data %>%
  group_by_at(metagroupnames[2])  %>% #filter(
    #sum(condgen == "SNI_WT") >= 2,
    #sum(condgen == "Sham_WT") >= 2) %>%

  t_test(
    formula = reformulate(metagroupnames[3], response = gene),
    ref.group = "Sham_WT"
  ) %>% 
  adjust_pvalue (method = "none") %>%
  add_significance () %>%
  #ungroup() %>%
  add_xy_position(
    data = box[[gene]]$data,
    formula = reformulate(metagroupnames[3], response = gene),
    x = metagroupnames[2],     # <-- x = RNA_snn_res.0.8 (the cluster)
    #group.by = metagroupnames[2],
    step.increase = 0.01
  )

    box[[gene]] = box[[gene]] + stat_pvalue_manual (
      stat.test, remove.bracket=FALSE,
    bracket.nudge.y = 0, hide.ns = T,
    x = metagroupnames[2],
    label = "p.adj.signif") 
    }

pdf (file.path ('Plots',paste0(gene,'_to_test_boxplots.pdf')), width=7, height=4)
wrap_plots(box)
dev.off()

pdf (file.path ('Plots','clustering_dimplot.pdf'), 4,4)
DimPlot (srt, group.by = 'celltype_predicted2')
dev.off()


### Analysis for Joseph Grieco ####

# Check Rest expression / regulon / cNMF module
gene = 'Rest'
# Run SCENIC plots ####
scenic_name = ''
motif_window = 'tss500bp'#'10kbp'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
metaGroupNames = c('sampleID','celltype_predicted','condition')
source (file.path(scrna_pipeline_dir, 'SCENIC_plots.R'))
grep (gene, colnames (auc_mtx))

# check in cNMF
nfeat = 7000
force=F
k_list = c(10:30)
k_selections = c(10:30)
cores= 100
cnmf_name = paste0('main')

### RUN consensus NMF ####
force=F
source (file.path (scrna_pipeline_dir,'cnmf_prepare_inputs.R'))     


### Import and format spectra files ####
cnmf_name = 'customized_vf_main'
k_selection=30
nfeat=10000 # Found with FindVariableFeatures
source (file.path (scrna_pipeline_dir,'cnmf_format_spectra_files.R')) 
cnmf_spectra_unique
sapply (cnmf_spectra_unique, function(x) 'Rest' %in% x)


mod = 'cNMF25'
module = cnmf_spectra_unique[[mod]]
module_mat = srt@assays$RNA@data[rownames(srt) %in% module, ]
module_mat = t(scale(t(module_mat)))

#palette_expression = colorRamp2(c(-2,0,2), c('blue','white','red'))
ha = HeatmapAnnotation (ct = srt$celltype_predicted2)
pdf (file.path ('Plots',paste0(mod,'_heatmap.pdf')), height=9)
Heatmap (module_mat, top_annotation = ha,
  row_names_gp = gpar(fontsize = 3), row_km = 5,
  column_names_gp = gpar(fontsize = 0), clustering_distance_rows = 'pearson')#,
  #col = palette_expression)
#VlnPlot (srt, features = 'Bmp6', group.by = 'Hughes_celltype')
dev.off()

# ct = unique(srt$Hughes_celltype)[grep ('neur',unique(srt$Hughes_celltype))]
# ct = 'Endothelial'
# ha = HeatmapAnnotation (
#   celltype = srt$Hughes_celltype[srt$Hughes_celltype %in% ct],
#   condition = srt$condition[srt$Hughes_celltype %in% ct],
#   sample = srt$sampleID[srt$Hughes_celltype %in% ct])

# pdf (file.path ('Plots',paste0(mod,'_only_',ct,'_heatmap.pdf')), height=13)
# Heatmap (module_mat[,srt$Hughes_celltype %in% ct], 
#   top_annotation = ha,
#   row_names_gp = gpar(fontsize = 5), row_km = 5,
#   column_names_gp = gpar(fontsize = 0), 
#   clustering_distance_rows = 'pearson'#,
#   #col = palette_expression
#   )
# #VlnPlot (srt, features = 'Bmp6', group.by = 'Hughes_celltype')
# dev.off()



#### Check Rest expression ####
gene = 'Rest'
gd = geneDot (
  seurat_obj = srt,
  #mat_norm = srt@assays$RNA@data,
  #mat_counts = srt@assays$RNA@counts,
  gene = gene, 
  x = 'celltype_predicted2', # Vector of metagroup1 of equal length to mat columns, If multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  y = 'condgen', # vector of metagroup2 of equal length to mat columns
  x_name = gene,
  assay='RNA',
  y_name = 'clusters',
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data = TRUE, # scale data when multiple genes are given
  plotcol = viridis::viridis(100),
  include_NA = TRUE,
  swap_axes = FALSE,
  returnDF = FALSE
  #... # arguments to pass to facet wrap
  )

pdf (file.path ('Plots',paste0(gene, '_dotplot.pdf')), height=3, width=5)
gd
dev.off()

modules = list(
  REST = c('Rest, Rcor1, Rcor2, Rcor3, Hdac1, Hdac2, Kdm1a, Ehmt2, Ehmt1, Sirt1, Sirt2'),
  ER = c('Atf6, Hspa5, Ddit3, Xbp1, Ern1, Atf4, Eif2ak3'),
  Ubi = c('Uchl1, Trim2, Trim9, Fbxo2, Fbxo3, Rnf41, Ube2n, Ube2l3, Ube3a, Nedd4, Nedd4l, Huwe1, Vcp, Derl1, Sel1l, Syvn1, Herpud1, Os9, Ubxn1, Pten'),
  Endo = c('Rab5a, Rab7a, Rab11a, Rab6a, Arf1, Arf4, Arl1, Tmed10, Tmed2, Tgoln1, Ap1b1, Cltc, Dnm1, Eea1, Vps35, Vps26a, Snx1, Snx2, Gga1, Gga2, Stx5, Stx6, Stx16, Yipf5, Sec22b')
    )

modules = lapply (modules, function(x) unlist(strsplit (x, ', ')))

ha = HeatmapAnnotation (ct = srt$celltype_predicted2, cond = srt$condgen)
pdf (file.path ('Plots',paste0(mod,'_heatmap.pdf')), height=3)
for (mod in names(modules))
  {
  module = modules[[mod]]
  module_mat = srt@assays$RNA@data[rownames(srt) %in% module, ]
  module_mat = t(scale(t(module_mat)))
  
draw(Heatmap (module_mat, top_annotation = ha,
    row_names_gp = gpar (fontsize = 6), #row_km = 5,
    column_names_gp = gpar (fontsize = 0), clustering_distance_rows = 'pearson',
  name = mod))
  }
dev.off()

# Cross-reference with cNMF and see if there is any enrichment for any
sapply (cnmf_spectra_unique, length)

ovmat = function (ovlist, df=FALSE, ov_threshold=0.5, compare_lists= NULL, palette=NULL)
  {
  if (!df) require (ComplexHeatmap) 
  comp_mat = sapply (ovlist, function(x) sapply(ovlist, function (y) sum (unique(x) %in% unique(y)) / min(c(length(unique(x)),length(unique(y))))))
  if (!is.null (compare_lists)) comp_mat = comp_mat[compare_lists[[1]],compare_lists[[2]]]
  if (df) 
    {
    return (comp_mat)
    } else {
    if (is.null(palette)) palette = viridis::mako(100)
    return (Heatmap (
      comp_mat, 
      col=palette,
      cell_fun = function(j, i, x, y, width, height, fill) {
      if (comp_mat[i,j] > ov_threshold) grid.text (sprintf("%.1f", comp_mat[i, j]), x, y, gp = gpar (fontsize=5))      
      }))      
    }
  }

ov = ovmat (ovlist = c(cnmf_spectra_unique, modules), compare_lists = list(names(cnmf_spectra_unique), names(modules)))
pdf (file.path ('Plots','cnmf_Joseph_programs_ov_heatmap.pdf'))
ov
dev.off()

ov_df = as.data.frame (do.call (cbind, lapply(cnmf_spectra_unique, 
  function(x) unlist(sapply (modules, function(y) sum(y%in%x))))))
ov_df$module = rownames(ov_df)
ov_df = gather (ov_df, cnmf, hit,1:(ncol(ov_df)-1))
ov_df$module = gsub ('\\d','',ov_df$module)

bp = ggplot (ov_df) + geom_bar (
  aes(x = cnmf, y= hit), 
  position = 'stack',
  stat = 'identity') + facet_wrap (~module)
pdf (file.path ('Plots','overlap_barplot.pdf'))
bp
dev.off()



### Run cellphoneDB ####
force=F
metaGroupName1 = 'sampleID'
metaGroupName2 = 'celltype_predicted2'
#metaGroupName2 = 'Hughes_celltype'
metaGroupName3 = 'condgen'
source (file.path (scrna_pipeline_dir, 'cellphoneDB_run.R'))

### generate plots from CPDB output ####
compGroups = c('SNI_WT','Sham_WT')
compGroups = c('SNI_WT','Sham_WT')
compGroups = c('SNI_KO','SNI_WT')
compGroups = c('Sham_KO','Sham_WT')
compGroups = c('SNI_WT','Sham_WT')

lr_filter = NULL
#celltype_filter = c('Bcells','DCs','Hypoxic_mac','Macrophages_1','Macrophages_2','Mast_cells','T_cells')
force=F
#celltype_filter = c('Schwann')
#celltype_filter = c('Myeloid')
celltype_filter = NULL
top_lr_filter = NULL
minCells = 20
pValThreshold = -log10(0.05)
source (file.path (scrna_pipeline_dir, 'cellphoneDB_downstream.R'))

gene = c('Adgrl2','Tenm3','Lrrtm4','Nrxn2','Lrfn5','Ptprf','Ptprs','Tgfb3','Tgfbr3','Btc','Erbb3')
srt_sub = srt[,srt$celltype_predicted2 %in% c('Schwann','NF1','NF2/NF3','NP','cLTMR1') & srt$condgen %in% c(compGroups)]
srt_sub$celltype_predicted2 = gsub ('\\/','_',srt_sub$celltype_predicted2)
gd = geneDot (
  seurat_obj = srt_sub,
  #mat_norm = srt@assays$RNA@data,
  #mat_counts = srt@assays$RNA@counts,
  gene = gene, 
  x = 'celltype_predicted2', # Vector of metagroup1 of equal length to mat columns, If multiple genes are specified this is ignored and genes would make the x axis of dotplot instead
  y = 'condgen', # vector of metagroup2 of equal length to mat columns
  x_name = gene,
  assay='RNA',
  y_name = 'clusters',
  min_expression = 0,
  facet_ncol = 5,
  lim_expression = NULL,
  scale.data = TRUE, # scale data when multiple genes are given
  plotcol = viridis::viridis(100),
  include_NA = TRUE,
  swap_axes = FALSE,
  returnDF = FALSE
  #... # arguments to pass to facet wrap
  )

pdf (file.path ('Plots',paste0(gene, '_dotplot.pdf')), height=6, width=7)
gd
dev.off()
