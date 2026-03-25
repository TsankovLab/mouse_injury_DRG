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
  nfeat = 2000, # number of variable genes to consider for dimentionality reduction
  sigPCs = 15,
  vars_to_regress='nCount_RNA',
  metaGroupNames = c('sampleID','condition'),
  res = c(0.2, 0.8, 2, 5,10,15) # denovo cluster resolutions 
  )

# Initiate pipeline
force = FALSE # re run pipeline from the beginning no matter if objects are found
subclustername = 'Schwann'
projdir = paste0 ("/ahg/regevdata/projects/ICA_Lung/Bruno/Serafini_prj/",proj_name,"_analysis/")
scrna_pipeline_dir = '/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/'
source (paste0(scrna_pipeline_dir,'master_scrna.R'))

### DEG on genotype ####
#srt$condition2 = ifelse (srt$condition %in% c('SNI_RGS4_KO','SNI_RGS4_WT'), 'SNI','SHAM')
force = FALSE
do.fgsea = TRUE
rankby = 'LFC'
logfcThreshold = 0
pvalAdjTrheshold = 0.05
topGenes = 50
addGene=NULL
#pval_column = 'pval' # specify this if adjusted pvalues are all non-significant
FeatureSets = list (
all = NULL#,  
#Il1_pathway = Il1_pathway, 
#TNF_pathway = TNF_pathway
)
#metaGroupName1 = paste0(reductionGraphSnn,'_res.0.8')
metaGroupName1 = 'celltype_predicted'
metaGroupName2 = 'condition'
deg2Ident = c('SHAM_RGS4_WT','SHAM_RGS4_KO')
deg2Ident = c('SNI_RGS4_WT','SNI_RGS4_KO')
source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/DEG2.R')

write.csv (deg2[deg2$p_val < 0.05,], paste0(paste(deg2Ident, collapse='_vs_'),'_lfc_0_DEG_genes_padj.csv'))
