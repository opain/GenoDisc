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
        "HM_BETA",
        "LOG_ODDS"
    ), 
    OR = c(
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
#    NEF = c(
#        "NEF",
#        "NEFF",
#        "NEFFECTIVE",
#        "NE"
#    ),
    FREQ = c(
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

head_interp<-function(sub_ss){

    # Read in the header and interpret column names
    sub_header<-names(sub_ss)

    # Remove columns that are all NA
    sub_ss_comp<-sub_ss[,apply(sub_ss, 2, function(x) !(all(is.na(x)))), with=F]
    sub_header_comp<-names(sub_ss_comp)

    int_header <- sub_header_comp
    for(i in names(ss_head_dict)){
        int_header[toupper(int_header) %in% ss_head_dict[[i]]] <- i
    }
    int_header[!(toupper(int_header) %in% unlist(ss_head_dict))]<-NA

    # Show original and interpreted header
    header_interp <- data.frame(Original = sub_header_comp,
                                Interpreted = int_header)

    # Show columns that are ignored due to be irrelevant, duplicated, or all NA
    header_interp$Keep<-!(
        !(int_header %in% names(ss_head_dict)) | 
        duplicated(int_header)
    )

    # Insert reason it was ignored
    header_interp$Reason<-NA
    header_interp$Reason[!(int_header %in% names(ss_head_dict))]<-'Not recognised'
    header_interp$Reason[duplicated(int_header)]<-'Duplicated'

    # Show columns ignored due to missingness
    for(i in sub_header[!(sub_header %in% header_interp$Original)]){
    header_interp<-rbind(header_interp, data.frame(Original = i,
                                                    Interpreted = NA,
                                                    Keep=F,
                                                    Reason = 'First 1000 rows NA'))
    }

    header_interp[is.na(header_interp)]<-'NA'
    header_interp$Keep<-as.character(header_interp$Keep)

    # Insert description of each column after interpretation
    header_labels<-data.frame(Interpreted=c('SNP','CHR','BP','A1','A2','P','OR','BETA','Z','SE','N','N_CAS','N_CON','NEF','FREQ','FRQ_A','FREQ_U','INFO'),
                            Description=c("RSID for variant",
                                        "Chromosome number",
                                        "Base pair position",
                                        "Allele 1 (effect allele)",
                                        "Allele 2",
                                        "P-value of association",
                                        "Odds ratio effect size",
                                        "BETA effect size",
                                        "Z-score",
                                        "Standard error of log(OR) or BETA",
                                        "Total sample size",
                                        "Number of cases",
                                        "Number of controls",
                                        "Effective sample size",
                                        "Allele frequency",
                                        "Allele frequency in cases",
                                        "Allele frequency in controls",
                                        "Imputation quality"))

    header_interp<-merge(header_interp, header_labels, by='Interpreted', all.x=T)
    header_interp$Keep<-factor(header_interp$Keep, levels=c('TRUE','FALSE'))
    header_interp<-header_interp[order(header_interp$Keep),]
    header_interp<-header_interp[,c('Original','Interpreted','Keep','Reason','Description')]

    return(header_interp)
}

format_header<-function(sumstats, log_file = NULL){

  header_interpretation<-head_interp(sumstats[1:1000,])
  header_interpretation<-header_interpretation[match(names(sumstats), header_interpretation$Original),]

  log_add(log_file = log_file, message = '---------------')
  
  if(is.null(log_file)){
    print.data.frame(header_interpretation, row.names=F, quote=F, right=F)
  } else {
    sink(file = log_file, append = T)
    print.data.frame(header_interpretation, row.names=F, quote=F, right=F)
    sink()    
  }

  log_add(log_file = log_file, message = '---------------')

  # Throw an error if required columns are not available
  error<-F
  messages <- NULL
  if(!all(c('A1','A2') %in% header_interpretation$Interpreted)){
    messages <- c(messages, "Error: Missing A1 and A2 information.")
    error<-T
  }
  if(!any(c('OR','BETA','Z') %in% header_interpretation$Interpreted)){
    messages <- c(messages, "Error: Effect size (BETA, OR, Z) must be present.")
    error<-T
  }
  if(!('SNP' %in% header_interpretation$Interpreted) & !(all(c('CHR','BP') %in% header_interpretation$Interpreted))){
    messages <- c(messages, "Error: Either SNP, or CHR and BP, must be present.")
    error<-T
  }
  if(any(c('N_CAS','N_CON') %in% header_interpretation$Interpreted) & sum(c('N_CAS','N_CON') %in% header_interpretation$Interpreted) == 1 & !('N' %in% names(GWAS))){
    messages <- c(messages, "Error: Either N, or both N_CAS and N_CON must be present.")
    error<-T
  }
  if(any(c('FRQ_A','FRQ_U') %in% header_interpretation$Interpreted) & sum(c('FRQ_A','FRQ_U') %in% header_interpretation$Interpreted) == 1){
    messages <- c(messages, "Warning: Both FRQ_A and FRQ_U must be present for either to be considered.")
  }

  if(!is.null(messages)){
    log_add(log_file = log_file, message = messages)

    if(error){
      stop(paste(messages, sep='\n'))
    }
  }

  # Remove columns that are not interpreted
  sumstats <- sumstats[, as.logical(header_interpretation$Keep), with=F]

  # Update sumstat header
  names(sumstats)<-header_interpretation$Interpreted[as.logical(header_interpretation$Keep)]

  return(sumstats)
}

# Calculate average allele frequency across cases and controls
calc_mean_freq<-function(sumstats, sampling = NA, log_file = NULL){
  if(all(c('FRQ_A','FRQ_U') %in% names(GWAS))){
    if(all(c('FRQ_A','FRQ_U','N_CAS','N_CON') %in% names(GWAS))){
        FREQ<- ((sumstats$FRQ_A*sumstats$N_CAS) + (sumstats$FRQ_U*sumstats$N_CON))/(sumstats$N_CAS + sumstats$N_CON)
        log_add(log_file = log_file, message = 'Average allele frequency calculated weighted by N_CAS and N_CON.')
    }
    if(all(c('FRQ_A','FRQ_U') %in% names(sumstats)) & any(!(c('FRQ_A','FRQ_U') %in% names(sumstats))) & !(is.na(sampling))){
        FREQ<- (sumstats$FRQ_A*sampling) + (sumstats$FRQ_U*(1-sampling))
        log_add(log_file = log_file, message = 'Average allele frequency calculated weighted by sampling fraction.')
    }
    if(all(c('FRQ_A','FRQ_U') %in% names(sumstats)) & any(!(c('FRQ_A','FRQ_U') %in% names(sumstats))) & is.na(sampling)){
        FREQ<- mean(sumstats$FRQ_A, sumstats$FRQ_U)
        log_add(log_file = log_file, message = 'Average allele frequency calculated assuming N_CAS == N_CON.')
    }
  } else {
    FREQ<-NULL
  }

  return(FREQ)
}

# Create IUPAC code function
snp_iupac<-function(x=NA, y=NA){
  if(length(x) != length(y)){
    print('x and y are different lengths')
  } else {
    iupac<-rep(NA, length(x))
    iupac[x == 'A' & y =='T' | x == 'T' & y =='A']<-'W'
    iupac[x == 'C' & y =='G' | x == 'G' & y =='C']<-'S'
    iupac[x == 'A' & y =='G' | x == 'G' & y =='A']<-'R'
    iupac[x == 'C' & y =='T' | x == 'T' & y =='C']<-'Y'
    iupac[x == 'G' & y =='T' | x == 'T' & y =='G']<-'K'
    iupac[x == 'A' & y =='C' | x == 'C' & y =='A']<-'M'
    return(iupac)
  }
}

# Determine of target using reference CHR and BP information
detect_build<-function(ref, targ, log_file = NULL){
    if(!(all(c('CHR','BP') %in% names(targ)))){
        stop('CHR and BP columns must be present in targ.\n')
    }
    if(!(all(c("CHR","A1","A2") %in% names(targ)) & any(grepl('BP_GRCh', names(ref))))){
        stop('CHR, A1, A2 and BP coordinates must be present in ref.\n')
    }

    # Insert IUPAC codes if missing
    if(!('IUPAC' %in% names(targ))){
        targ$IUPAC <- snp_iupac(targ$A1, targ$A2)
    }
    if(!('IUPAC' %in% names(ref))){
        ref$IUPAC <- snp_iupac(ref$A1, ref$A2)
    }

    builds<-as.numeric(gsub('BP_GRCh','',names(ref)[grepl('BP_GRCh', names(ref))]))

    # Check target-ref condordance of BP across builds
    ref$CHR<-as.character(ref$CHR)
    targ$CHR<-as.character(targ$CHR)

    matched<-list()
    target_build<-NA
    for(build_i in builds){
        matched<-merge(targ, ref, by.x=c('CHR','BP','IUPAC'), by.y=c('CHR',paste0('BP_GRCh', build_i),'IUPAC'))
        overlap<-nrow(matched)/nrow(targ)
        if(overlap > 0.7){
            target_build<-paste0('GRCh',build_i)
        }
        log_add(log_file = log_file, message = paste0('GRCh',build_i,' match: ',round(overlap*100, 2),'%'))
    }

    return(target_build)
}

# Create function to change allele to complement
snp_allele_comp<-function(x=NA){
  x_new<-x
  x_new[x == 'A']<-'T'
  x_new[x == 'T']<-'A'
  x_new[x == 'G']<-'C'
  x_new[x == 'C']<-'G'
  x_new[!(x %in% c('A','T','G','C'))]<-NA
  return(x_new)
}

# Create effective sample size calculator
neff<-function(ncas, ncon){
  return(4/((1/ncas)+(1/ncon)))
}

# Retain only non-ambiguous SNPs
remove_ambig<-function(dat){
  if(!('IUPAC' %in% names(dat))){
    if(!(all(c('A1','A2') %in% names(dat)))){
      stop('Either IUPAC, or A1 and A2 must be present')
    } else {
      iupac_codes <- snp_iupac(dat$A1, dat$A2)
      subset_condition <- iupac_codes %in% c('R', 'Y', 'K', 'M')
    }
  } else {
    subset_condition <- dat$IUPAC %in% c('R', 'Y', 'K', 'M')
  }

  return(dat[subset_condition, ])
}

# Identifies IUPAC codes in targ that are on opposite strand to ref
detect_strand_flip<-function(targ, ref){
    flipped<-(
        (targ == 'R' & ref == 'Y') | 
        (targ == 'Y' & ref == 'R') | 
        (targ == 'K' & ref == 'M') |
        (targ == 'M' & ref == 'K'))

    return(flipped)
}

ref_harmonise<-function(targ, ref_rds, log_file = NULL){

  ref_22<-readRDS(file = paste0(ref_rds, 22, '.rds'))
  if(!(all(c('CHR','A1','A2','REF.FRQ','IUPAC') %in% names(ref_22)) & any(grepl('BP_GRCh', names(ref_22))))){
      print(names(ref_22))
      stop('CHR, A1, A2, REF.FRQ, IUPAC and BP coordinates must be present in ref_rds files.\n')
  }

  # Check whether CHR and BP information are present
  chr_bp_avail<-sum(c('CHR','BP') %in% names(targ)) == 2 
  
  # Check whether RSIDs are available for majority of SNPs in GWAS
  rsid_avail<-(sum(grepl('rs', targ$SNP)) > 0.9*length(targ$SNP))

  targ_matched<-NULL
  flip_logical_all<-NULL
  chrs<-c(1:22)

  if(chr_bp_avail){
    log_add(log_file = log_file, message = 'Merging sumstats with reference using CHR, BP, A1, and A2')

    ###
    # Determine build
    ###

    target_build <- detect_build( ref = ref_22, 
                                  targ = targ[targ$CHR == 22,],
                                  log_file = log_file)
    
    if(!is.na(target_build)){
      for(i in chrs){
        print(i)
        
        # Read reference data
        ref_i<-readRDS(file = paste0(ref_rds,i,'.rds'))
        
        # Retain only non-ambiguous SNPs
        ref_i<-remove_ambig(ref_i)

        # Rename columns prior to merging with target
        names(ref_i)<-paste0('REF.',names(ref_i))
        names(ref_i)[names(ref_i) == 'REF.REF.FRQ']<-'REF.FREQ'
        ref_i<-ref_i[, c('REF.CHR','REF.SNP',paste0('REF.BP_',target_build),'REF.A1','REF.A2','REF.IUPAC','REF.FREQ'), with=F]
        
        # Subset chromosome i from target
        targ_i<-targ[targ$CHR == i,]

        # Merge target and reference by BP
        ref_target<-merge(targ_i, ref_i, by.x='BP', by.y=paste0('REF.BP_',target_build))
        
        # Identify targ-ref strand flips, and flip target
        flip_logical<-detect_strand_flip(targ = ref_target$IUPAC, ref = ref_target$REF.IUPAC)
        flip_logical_all<-c(flip_logical_all, flip_logical)

        flipped<-ref_target[flip_logical,]
        flipped$A1<-snp_allele_comp(flipped$A1)
        flipped$A2<-snp_allele_comp(flipped$A2)      
        flipped$IUPAC<-snp_iupac(flipped$A1, flipped$A2)
        
        # Identify SNPs that have matched IUPAC
        matched<-ref_target[ref_target$IUPAC == ref_target$REF.IUPAC,]
        matched<-rbind(matched, flipped)
        
        # Flip REF.FREQ if alleles are swapped
        matched$REF.FREQ[matched$A1 != matched$REF.A1]<-1-matched$REF.FREQ[matched$A1 != matched$REF.A1]
        
        # Retain reference SNP and REF.FREQ data
        matched<-matched[, names(matched) %in% c('CHR','BP','REF.SNP','A1','A2','BETA','SE','OR','Z','FREQ','REF.FREQ','N','INFO','P'), with=F]
        names(matched)[names(matched) == 'REF.SNP']<-'SNP'

        targ_matched<-rbind(targ_matched, matched)
      }
    }

    if(is.na(target_build) & !(rsid_avail)){
      log_add(log_file = log_file, message = 'Error: Target build could not be determined and SNP IDs unavailable.')
      stop('Target build could not be determined and SNP IDs unavailable.\n')
    }
    if(is.na(target_build) & rsid_avail){
      log_add(log_file = log_file, message = 'Target build could not be determined from CHR and BP data.')
    }
  }

  if(is.na(target_build) & rsid_avail){
    log_add(log_file = log_file, message = 'Using SNP, A1 and A2 to merge with the reference.')
    
    for(i in chrs){
      print(i)
      
      # Read reference data
      ref_i<-readRDS(file = paste0(ref_rds,i,'.rds'))
      
      # Retain only non-ambiguous SNPs
      ref_i<-remove_ambig(ref_i)

      # Rename columns prior to merging with target
      names(ref_i)<-paste0('REF.',names(ref_i))
      names(ref_i)[names(ref_i) == 'REF.REF.FRQ']<-'REF.FREQ'
      ref_i<-ref_i[, c('REF.CHR','REF.SNP',paste0('REF.BP_',target_build),'REF.A1','REF.A2','REF.IUPAC','REF.FREQ'), with=F]
    
      # Merge target and reference by SNP ID
      ref_target<-merge(targ, tmp, by='SNP')
      
      # Identify targ-ref strand flips, and flip target
      flip_logical<-detect_strand_flip(targ = ref_target$IUPAC, ref = ref_target$REF.IUPAC)
      flip_logical_all<-c(flip_logical_all, flip_logical)

      flipped<-ref_target[flip_logical,]
      flipped$A1<-snp_allele_comp(flipped$A1)
      flipped$A2<-snp_allele_comp(flipped$A2)      
      flipped$IUPAC<-snp_iupac(flipped$A1, flipped$A2)
      
      # Identify SNPs that have matched IUPAC
      matched<-ref_target[ref_target$IUPAC == ref_target$REF.IUPAC,]
      matched<-rbind(matched, flipped)
      
      # Flip REF.FREQ if alleles are swapped
      matched$REF.FREQ[matched$A1 != matched$REF.A1]<-1-matched$REF.FREQ[matched$A1 != matched$REF.A1]
      
      # Retain reference CHR and BP_GRCh37 data
      matched<-matched[, names(matched) %in% c('REF.CHR','REF.BP_GRCh37','SNP','A1','A2','BETA','SE','OR','Z','FREQ','REF.FREQ','N','INFO','P'), with=F]
      names(matched)[names(matched) == 'REF.CHR']<-'CHR'
      names(matched)[names(matched) == 'REF.BP_GRCh37']<-'BP'

      targ_matched<-rbind(targ_matched, matched)
    }
  }

  log_add(log_file = log_file, message = paste0('After matching variants to the reference ,',nrow(targ_matched),' variants remain.'))
  log_add(log_file = log_file, message = paste0(sum(flip_logical_all), ' variants were flipped to match reference.'))

  return(targ_matched)
}

# Remove variants with INFO less than threshold
filter_info<-function(targ, thresh, log_file = NULL){
  if(sum(names(targ) == 'INFO') == 1){
    targ<-targ[targ$INFO >= thresh,]
    log_add(log_file = log_file, message = paste0('After removal of SNPs with INFO < ',thresh,', ',nrow(targ),' variants remain.'))
  } else {
    log_add(log_file = log_file, message = 'INFO column is not present.')
  }
  return(targ)
}

# Remove variants with MAF less than threshold
# If ref = T, then REF.FREQ column is used
filter_maf<-function(targ, thresh, ref = F, log_file = NULL){
  if(ref == F){
    if(any(names(targ) == 'FREQ')){
        targ<-targ[targ$FREQ >= thresh & targ$FREQ <= (1-thresh),]
        log_add(log_file = log_file, message = paste0('After removal of SNPs with reported MAF < ',thresh,', ',nrow(targ),' variants remain.'))    
    } else {
        log_add(log_file = log_file, message = 'FREQ  column is not present.')
    }
  } else {
    if(any(names(targ) == 'REF.FREQ')){
        targ<-targ[targ$REF.FREQ >= thresh & targ$REF.FREQ <= (1-thresh),]
        log_add(log_file = log_file, message = paste0('After removal of SNPs with reference MAF < ',thresh,', ',nrow(targ),' variants remain.'))    
    } else {
        log_add(log_file = log_file, message = 'REF.FREQ  column is not present.')
    }    
  }

  return(targ)
}

# Remove variants with reported and reference allele frequency difference greater than threshold
discord_maf<-function(targ, thresh, log_file, plot_file = NA){
  if(sum(names(targ) == 'FREQ') == 1){
    targ$diff<-abs(targ$FREQ-targ$REF.FREQ)
    
    if(!is.na(plot_file)){
      png(plot_file, unit='px', res=300, width=1200, height=1200)
        plot(targ$REF.FREQ[targ$diff > thresh],targ$FREQ[targ$diff > thresh], xlim=c(0,1), ylim=c(0,1), xlab='Reference Allele Frequency', ylab='Sumstat Allele Frequency')
        abline(coef = c(0,1))
      dev.off()
    }

    targ<-targ[targ$diff < thresh,]
    targ$diff<-NULL
    
    log_add(log_file = log_file, message = paste0('After removal of SNPs with absolute MAF difference of < ',thresh,', ',nrow(targ),' variants remain.'))
  } else {
    log_add(log_file = log_file, message = 'Reported MAF column is not present, so discordance with reference cannot be determined.')
  }
  return(targ)
}

# Remove rows where N is X SD above or below the median N
filter_n <- function(targ, n_sd = 3, log_file = NULL){
  if(length(unique(targ$N)) > 1){
    thresh <- n_sd*sd(targ$N)
    targ <- targ[targ$N < median(targ$N) + thresh & targ$N > median(targ$N) - thresh,]
    
    log_add(log_file = log_file, message = paste0('After removal of SNPs with N > ',median(targ$N)+thresh,' or < ',median(targ$N)-thresh,', ',nrow(targ),' variants remain.'))
  } else {
    log_add(log_file = log_file, message = 'N column is not present or invariant.')
  }
}

# Convert OR or Z into BETA if not present already
insert_beta <- function(targ, log_file = NULL){
  # If OR is present but BETA is not, convert OR to logOR and name it BETA
  if('OR' %in% names(targ) & !('BETA' %in% names(targ))){
    targ$BETA<-log(targ$OR)
    log_add(log_file = log_file, message = 'BETA column inserted based on log(OR).')
  }

  # If Z is present but BETA is not, calculate it based on sample size, allele frequency and Z score
  if('Z' %in% names(targ) & !('BETA' %in% names(targ))){
    if('FREQ' %in% names(targ)){
        frq_tmp<-targ$FREQ
        log_add(log_file = log_file, message = 'BETA and SE column inserted based on Z, FREQ, and N.')
      } else {
        frq_tmp<-targ$REF.FREQ
        log_add(log_file = log_file, message = 'BETA and SE column inserted based on Z, REF.FREQ, and N')
    }

    targ$SE <- 1/sqrt((2*frq_tmp)*(1-(frq_tmp))*(targ$N + (targ$Z^2)))
    targ$BETA <- targ$Z * targ$SE
  }
  return(targ)
}

# Insert SE from BETA and P value if not already present
insert_se <- function(targ, log_file = NULL){
  if(sum(names(targ) == 'SE') == 0){
    if(any(!(c('BETA','P') %in% names(targ)))){
      stop('BETA and P columns must be present to compute SE.\n')
    }
    
    z_tmp<-abs(qnorm(targ$P/2))
    targ$SE<-abs(targ$BETA/z_tmp)

    log_add(log_file = log_file, message = 'SE column inserted based on BETA and P.')
  }
  return(targ)
}

# Detect whether P values have been adjusted using genomic control
# If so, recompute using SE if available
avoid_gc<-function(targ, log_file = NULL){
  if(sum(names(targ) == 'SE') == 1){  
    targ$Z<-targ$BETA/targ$SE
    targ$P_check<-2*pnorm(-abs(targ$Z))
    targ$Z<-NULL

    if(abs(mean(targ$P[!is.na(targ$P_check)]) - mean(targ$P_check[!is.na(targ$P_check)])) > 0.01){
      targ$P<-targ$P_check
      targ$P_check<-NULL
      log_add(log_file = log_file, message = 'Genomic control detected. P-value recomputed using BETA and SE.')
    } else {
      log_add(log_file = log_file, message = 'Genomic control was not detected.')
      targ$P_check<-NULL
    }
  } else {
    log_add(log_file = log_file, message = 'SE column is not present, genomic control cannot be detected.')
  }
  return(targ)
}
