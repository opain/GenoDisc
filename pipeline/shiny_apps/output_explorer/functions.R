# Function to check email has valid format
is_valid_email_format <- function(email) {
  grepl("^[a-zA-Z0-9.]+@[a-zA-Z0-9.]+\\.[a-zA-Z]{2,}$", email)
}

# Create names vectors for selecting panels
fusion_panels<-value<-c(
  'Adipose_Subcutaneous',
  'Adipose_Visceral_Omentum',
  'Adrenal_Gland',
  'Artery_Aorta',
  'Artery_Coronary',
  'Artery_Tibial',
  'Brain_Amygdala',
  'Brain_Anterior_cingulate_cortex_BA24',
  'Brain_Caudate_basal_ganglia',
  'Brain_Cerebellar_Hemisphere',
  'Brain_Cerebellum',
  'Brain_Cortex',
  'Brain_Frontal_Cortex_BA9',
  'Brain_Hippocampus',
  'Brain_Hypothalamus',
  'Brain_Nucleus_accumbens_basal_ganglia',
  'Brain_Putamen_basal_ganglia',
  'Brain_Spinal_cord_cervical_c-1',
  'Brain_Substantia_nigra',
  'Breast_Mammary_Tissue',
  'Cells_EBV-transformed_lymphocytes',
  'Cells_Transformed_fibroblasts',
  'CMC.BRAIN.RNASEQ',
  'CMC.BRAIN.RNASEQ_SPLICING',
  'Colon_Sigmoid',
  'Colon_Transverse',
  'Esophagus_Gastroesophageal_Junction',
  'Esophagus_Mucosa',
  'Esophagus_Muscularis',
  'Heart_Atrial_Appendage',
  'Heart_Left_Ventricle',
  'Liver',
  'Lung',
  'METSIM.ADIPOSE.RNASEQ',
  'Minor_Salivary_Gland',
  'Muscle_Skeletal',
  'Nerve_Tibial',
  'NTR.BLOOD.RNAARR',
  'Ovary',
  'Pancreas',
  'Pituitary',
  'Prostate',
  'Skin_Not_Sun_Exposed_Suprapubic',
  'Skin_Sun_Exposed_Lower_leg',
  'Small_Intestine_Terminal_Ileum',
  'Spleen',
  'Stomach',
  'Testis',
  'Thyroid',
  'Uterus',
  'Vagina',
  'Whole_Blood',
  'YFS.BLOOD.RNAARR'
)

fusion_panels<-data.frame(original=fusion_panels,
                          clean=fusion_panels)

gtex_fusion_panels<-fusion_panels[!grepl('CMC|YFS|NTR|METSIM', fusion_panels$original),]
nongtex_fusion_panels<-fusion_panels[grepl('CMC|YFS|NTR|METSIM', fusion_panels$original),]

gtex_fusion_panels$clean<-gsub('_',' ',gtex_fusion_panels$clean)
gtex_fusion_panels$type<-gsub(' .*','',gtex_fusion_panels$clean)
gtex_fusion_panel_types<-unique(gtex_fusion_panels$type[duplicated(gtex_fusion_panels$type)])
gtex_fusion_panels$clean[gtex_fusion_panels$type %in% gtex_fusion_panel_types]<-gsub("^.*?\\s",'',gtex_fusion_panels$clean[gtex_fusion_panels$type %in% gtex_fusion_panel_types])
gtex_fusion_panels$clean<-tolower(gtex_fusion_panels$clean)
gtex_fusion_panels$clean<-gsub('ebv','EBV', gtex_fusion_panels$clean)
gtex_fusion_panels$clean<-gsub('ba24','BA24', gtex_fusion_panels$clean)
gtex_fusion_panels$clean<-gsub('ba9','BA9', gtex_fusion_panels$clean)
gtex_fusion_panels$clean <- paste0(toupper(substr(gtex_fusion_panels$clean, 1, 1)), substr(gtex_fusion_panels$clean, 2, nchar(gtex_fusion_panels$clean)))
gtex_fusion_panels$clean<-paste0(gtex_fusion_panels$clean, " (GTEx)")
gtex_fusion_panels$clean[gtex_fusion_panels$type %in% gtex_fusion_panel_types]<-paste0(gtex_fusion_panels$type,' - ', gtex_fusion_panels$clean)[gtex_fusion_panels$type %in% gtex_fusion_panel_types]
gtex_fusion_panels$type<-NULL

nongtex_fusion_panels$clean[nongtex_fusion_panels$original == 'CMC.BRAIN.RNASEQ_SPLICING']<-"Brain - DLPFC - Splice (CMC)"
nongtex_fusion_panels$clean[nongtex_fusion_panels$original == 'CMC.BRAIN.RNASEQ']<-"Brain - DLPFC (CMC)"
nongtex_fusion_panels$clean[nongtex_fusion_panels$original == 'METSIM.ADIPOSE.RNASEQ']<-"Adipose (METSIM)"
nongtex_fusion_panels$clean[nongtex_fusion_panels$original == 'NTR.BLOOD.RNAARR']<-"Blood (NTR)"
nongtex_fusion_panels$clean[nongtex_fusion_panels$original == 'YFS.BLOOD.RNAARR']<-"Blood (YFS)"

fusion_twas_panel_names <- rbind(gtex_fusion_panels, nongtex_fusion_panels)
names(fusion_twas_panel_names)<-c('original','clean')

non_fusion_panels<-data.frame(original='psychencode',
                              clean='Brain - DLPFC (PsychENCODE)')

twas_panel_names<-rbind(fusion_twas_panel_names, non_fusion_panels)
twas_panel_names <- setNames(twas_panel_names$original, twas_panel_names$clean)

pwas_panel_names <- setNames(c('rosmap','banner'), c("Brain - DLPFC (ROSMAP)","Brain - DLPFC (Banner)"))

smr_expr_panel_names <- data.frame(
  original = c(
    'smr_expression_panel_psychencode',
    'smr_expression_panel_metabrain_basalganglia',
    'smr_expression_panel_metabrain_cerebellum',
    'smr_expression_panel_metabrain_cortex',
    'smr_expression_panel_metabrain_hippocampus',
    'smr_expression_panel_metabrain_spinalcord',
    'smr_expression_panel_eqtlgen'
  ),
  clean = c(
    "Brain - DLPFC (PsychENCODE)",
    "Brain - Basalganglia (MetaBrain)",
    "Brain - Cerebellum (MetaBrain)",
    "Brain - Cortex (MetaBrain)",
    "Brain - Hippocampus (MetaBrain)",
    "Brain - Spinalcord (MetaBrain)",
    "Blood (eQTLGen)"
  )
)

smr_expr_panel_names <- setNames(smr_expr_panel_names$original, smr_expr_panel_names$clean)

smr_protein_panel_names <- data.frame(
  original = c(
    'smr_protein_panel_rosmap'
  ),
  clean = c(
    "Brain - DLPFC (ROSMAP)"
  )
)
smr_protein_panel_names <- setNames(smr_protein_panel_names$original, smr_protein_panel_names$clean)

# Create a named vector for the populations available
populations<-data.frame(original=c('EUR'),
                        cleaned=c('European'))
populations <- setNames(populations$original, populations$clean)



# Read a GWAS QC statistic (lambda_gc, max_chi2, n_sig_snp) from a gwas_qc block.
# These moved from the FOCUS munge log to the sumstat_cleaner log when the FOCUS
# munge step was dropped, so results packaged before that change carry them under
# focus_dat instead of cleaner_dat.
gd_qc_stat <- function(gwas_qc, field) {
  val <- gwas_qc$cleaner_dat$val[[field]]
  if (length(val) != 1 || is.na(val)) val <- gwas_qc$focus_dat$val[[field]]
  if (length(val) != 1) return(NA_real_)
  val
}
