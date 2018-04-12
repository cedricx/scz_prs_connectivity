## import data ##
wkdir <- '~/Google Drive/TDSlab/SCZ_gene_imaging/'
prs <- read.csv(paste0(wkdir,'genotype/pnc_EA_prs_20170501.csv'))
demographic <- read.csv(paste0(wkdir,'phenotype/EuropeanAmerican/demographics/n9498_demographics_go1_20161212_EA.csv'))
snp108 <- read.csv(paste0(wkdir,'genotype/illumina_EA_poly_scz108_converted.csv'))
cnb <- read.csv(paste0(wkdir,'phenotype/EuropeanAmerican/cnb/n9498_cnb_zscores_fr_20170202_EA.csv'))
cnbar <- read.csv(paste0(wkdir,'phenotype/EuropeanAmerican/cnb/n9498_cnb_zscores_frar_20170202_EA.csv'))
cnb_factor <- read.csv(paste0(wkdir,'phenotype/EuropeanAmerican/cnb/n9498_cnb_factor_scores_fr_20170202_EA.csv'))
cnbar_factor <- read.csv(paste0(wkdir,'phenotype/EuropeanAmerican/cnb/n9498_cnb_factor_scores_frar_20170202_EA.csv'))
bifactor <- read.csv(paste0(wkdir,'phenotype/EuropeanAmerican/clinical/n9498_goassess_itemwise_bifactor_scores_20161219_EA.csv'))    
envr <- read.csv(paste0(wkdir,'phenotype/EuropeanAmerican/environment/n9498_go1_environment_factor_scores_tymoore_20150909_EA.csv'))

## REML functions ##
prs_gam10<-function(vab_of_int, prs_df){
  out<-gam(vab_of_int ~ s(ageAtCnb1) +  sex + scz_prs +
             pc1 + pc2 + pc3 + pc4 + pc5 +
             pc6 + pc7 + pc8 + pc9 + pc10, data = prs_df, method="REML")
}


voi_reml <-function(voi) {
  prs_voi <- merge(prs,voi, by='bblid')
  range = 25:dim(prs_voi)[2]
  prs_voi_reml<-lapply(range,function(vab_of_int) prs_gam10(prs_voi[,vab_of_int],prs_voi))
  prs_voi_pval<-sapply(prs_voi_reml, function(reml_result) summary(reml_result)$p.pv)
  colnames(prs_voi_pval) <- colnames(prs_voi)[range]
  voi_of_sig <- prs_voi_pval['scz_prs',which(prs_voi_pval['scz_prs',] <0.05)]
  out <- list(data = prs_voi, pval = prs_voi_pval, sig_voi = voi_of_sig)
}

## combine PRS and demographics ##
prs <- merge(prs, demographic, by = 'bblid')
prs$sex <- as.ordered(as.factor(prs$sex))
## PRS vs CNB
cnb_out <-voi_reml(cnb)

## PRS vs CNB_ar
cnbar_out <-voi_reml(cnbar)

## PRS vs CNB_factor
cnbfactor_out <-voi_reml(cnb_factor)
cnbfactorar_out <-voi_reml(cnbar_factor)


## PRS vs bifactor
cnbbifactor_out <-voi_reml(bifactor)

## PRS vs env
cnbenv_out <-voi_reml(envr)

