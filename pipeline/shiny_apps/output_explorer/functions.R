# Create list showing possible column names in sumstats
ss_head_dict<-list(
  SNP=c(
    "SNP",
    "PREDICTOR",
    "SNPID",
    "MARKERNAME",
    "MARKER_NAME",
    "SNPTESTID",
    "ID_DBSNP49",
    "RSID",
    "ID",
    "RS_NUMBER",
    "MARKER",
    "RS",
    "RSNUMBER",
    "RS_NUMBERS",
    "SNP.NAME",
    "SNP ID",
    "SNP_ID",
    "LOCATIONALID",
    "ASSAY_NAME",
    "VARIANT_ID",
    "HM_RSID"
  ),
  A1 = c(
    "A1",
    "ALLELE1",
    "ALLELE_1",
    "INC_ALLELE",
    "EA",
    "A1_EFFECT",
    "REF",
    "EFFECT_ALLELE",
    "RISK_ALLELE",
    "EFFECTALLELE",
    "EFFECT_ALL",
    "REFERENCE_ALLELE",
    "REF_ALLELE",
    "REFERENCEALLELE",
    "EA",
    "ALLELE_1",
    "INC_ALLELE",
    "ALLELE1",
    "A",
    "A_1",
    "CODED_ALLELE",
    "TESTED_ALLELE",
    "HM_EFFECT_ALLELE"
  ),
  A2 = c(
    "A2",
    "ALLELE2",
    "ALLELE_2",
    "OTHER_ALLELE",
    "NON_EFFECT_ALLELE",
    "DEC_ALLELE",
    "OA",
    "NEA",
    "ALT",
    "A2_OTHER",
    "NONREF_ALLELE",
    "NEFFECT_ALLELE",
    "NEFFECTALLELE",
    "NONEFFECT_ALLELE",
    "OTHER_ALL",
    "OTHERALLELE",
    "NONEFFECTALLELE",
    "ALLELE0",
    "ALLELE_0",
    "ALT_ALLELE",
    "A_0",
    "NONCODED_ALLELE",
    "HM_OTHER_ALLELE"
  ),
  BETA = c(
    "BETA",
    "B",
    "EFFECT_BETA",
    "EFFECT",
    "EFFECTS",
    "SIGNED_SUMSTAT",
    "EST",
    "GWAS_BETA",
    "EFFECT_A1",
    "EFFECTA1",
    "EFFECT_NW",
    "HM_BETA"
  ), 
  OR = c(
    "OR",
    "LOG_ODDS",
    "OR",
    "ODDS-RATIO",
    "ODDS_RATIO",
    "ODDSRATIO",
    "OR(MINALLELE)",
    "OR.LOGISTIC",
    "OR_RAN",
    "OR(A1)",
    "HM_ODDS_RATIO"
  ),
  SE = c(
    "SE",
    "STDER",
    "STDERR",
    "STD",
    "STANDARD_ERROR",
    "OR_SE",
    "STANDARDERROR",
    "STDERR_NW",
    "META.SE",
    "SE_DGC",
    "SE.2GC"
  ),
  Z = c(
    "Z",
    "ZSCORE",
    "Z-SCORE",
    "ZSTAT",
    "ZSTATISTIC",
    "GC_ZSCORE",
    "BETAZSCALE"
  ),
  INFO = c(
    "INFO",
    "IMPINFO",
    "IMPQUALITY",
    "INFO.PLINK",
    "INFO_UKBB",
    "INFO_UKB"
  ),
  P = c(
    "P",
    "PVALUE",
    "PVAL",
    "P_VALUE",
    "GC_PVALUE",
    "WALD_P",
    "P.VAL",
    "GWAS_P",
    "P-VALUE",
    "P-VAL",
    "FREQUENTIST_ADD_PVALUE",
    "P.VALUE",
    "P_VAL",
    "SCAN-P",
    "P.LMM",
    "META.PVAL",
    "P_RAN",
    "P.ADD",
    "P_BOLT_LMM"
  ),
  N = c(
    "N",
    "WEIGHT",
    "NCOMPLETESAMPLES",
    "TOTALSAMPLESIZE",
    "TOTALN",
    "TOTAL_N",
    "N_COMPLETE_SAMPLES",
    "N_TOTAL",
    "N_SAMPLES",
    "N_ANALYZED",
    "NSAMPLES",
    "SAMPLESIZE",
    "SAMPLE_SIZE",
    "TOTAL_SAMPLE_SIZE",
    "TOTALSAMPLESIZE"
  ),
  N_CAS = c(
    "N_CAS",
    "NCASE",
    "N_CASE",
    "N_CASES",
    "NCAS",
    "NCA",
    "NCASES",
    "CASES",
    "CASES_N"
  ),
  N_CON = c(
    "N_CON",
    "NCONTROL",
    "N_CONTROL",
    "N_CONTROLS",
    "NCON",
    "NCO",
    "N_CON",
    "NCONTROLS",
    "CONTROLS",
    "CONTROLS_N"
  ),
  NEF = c(
    "NEF",
    "NEFF",
    "NEFFECTIVE",
    "NE"
  ),
  FRQ = c(
    "FRQ",
    "FREQ",
    "MAF",
    "AF",
    "CEUAF",
    "FREQ1",
    "EAF",
    "FREQ1.HAPMAP",
    "FREQALLELE1HAPMAPCEU",
    "FREQ.ALLELE1.HAPMAPCEU",
    "EFFECT_ALLELE_FREQ",
    "FREQ.A1",
    "F_A",
    "F_U",
    "FREQ_A",
    "FREQ_U",
    "MA_FREQ",
    "MAF_NW",
    "FREQ_A1",
    "A1FREQ",
    "CODED_ALLELE_FREQUENCY",
    "FREQ_TESTED_ALLELE_IN_HRS",
    "EAF_HRC",
    "EAF_UKB",
    "EFFECT_ALLELE_FREQUENCY",
    "HM_EFFECT_ALLELE_FREQUENCY"
  ),
  FRQ_A = c(
    "FREQ_A",
    "F_A"
  ),
  FRQ_U = c(
    "FREQ_U",
    "F_U"
  ),
  CHR = c(
    "CHR",
    "CH",
    "CHROMOSOME",
    "CHROM",
    "CHR_BUILD38",
    "CHR_BUILD37",
    "CHR_BUILD36",
    "CHR_B38",
    "CHR_B37",
    "CHR_B36",
    "CHR_ID",
    "SCAFFOLD",
    "HG19CHR",
    "CHR.HG19",
    "CHR_HG19",
    "HG18CHR",
    "CHR.HG18",
    "CHR_HG18",
    "CHR_BP_HG19B37",
    "HG19CHRC",
    "HM_CHROM"
  ),
  BP = c(
    "ORIGBP",
    "BP",
    "POS",
    "POSITION",
    "LOCATION",
    "PHYSPOS",
    "GENPOS",
    "CHR_POSITION",
    "POS_B38",
    "POS_BUILD38",
    "POS_B37",
    "POS_BUILD37",
    "BP_HG19B37",
    "POS_B36",
    "POS_BUILD36",
    "POS.HG19",
    "POS.HG18",
    "POS_HG19",
    "POS_HG18",
    "BP_HG19",
    "BP_HG18",
    "BP.GRCH38",
    "BP.GRCH37",
    "POSITION(HG19)",
    "POSITION(HG18)",
    "POS(B38)",
    "POS(B37)",
    "BASE_PAIR_LOCATION",
    "HM_POS"
  )
)

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


