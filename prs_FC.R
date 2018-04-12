## match subjects ##
## set directories ##
local_wkdir <- '~/Google Drive/TDSlab/SCZ_gene_imaging/'
remote_wkdir <- '~/Desktop/BBL/data/joy/BBL/studies/pnc/'

## import demographcis ##
prs_demo <- read.csv(paste0(local_wkdir,'phenotype/EuropeanAmerican/demographics/n9498_demographics_go1_20161212_EA.csv'))
img_demo <- read.csv(paste0(remote_wkdir,'n9498_dataFreeze/demographics/n9498_demographics_go1_20161212.csv'))


## import cnb ##
prs_cnb <- read.csv(paste0(local_wkdir,'phenotype/EuropeanAmerican/cnb/n9498_cnb_zscores_fr_20170202_EA.csv'))
img_cnb <- read.csv(paste0(remote_wkdir,'n9498_dataFreeze/cnb/n9498_cnb_zscores_fr_20170202.csv'))

## merge cnb and demo ##
prs_cnb_demo <- merge(prs_demo,prs_cnb, by = 'bblid')
img_cnb_demo <- merge(img_demo,img_cnb, by = 'bblid')


prs_img_matched <- merge(prs_cnb_demo,img_cnb_demo, by = c(colnames(prs_cnb)[2:29],colnames(prs_demo)[2:11]))
dim(prs_img_matched)

## test if merged matches ##
prs_cnb_factor <- read.csv(paste0(local_wkdir,'phenotype/EuropeanAmerican/cnb/n9498_cnb_factor_scores_fr_20170202_EA.csv'))
img_cnb_factor <- read.csv(paste0(remote_wkdir,'n9498_dataFreeze/cnb/n9498_cnb_factor_scores_fr_20170202.csv'))

prs_cnb_factor_match <- merge(prs_cnb_factor, prs_img_matched, by.x = 'bblid', by.y = 'bblid.x')
img_cnb_factor_match <- merge(img_cnb_factor, prs_img_matched, by.x = 'bblid', by.y = 'bblid.y')

prs_img_cnb_factor_match <- merge(prs_cnb_factor_match,img_cnb_factor_match, by.x = 'bblid', by.y = 'bblid.x')
plot(prs_img_cnb_factor_match$NAR_Overall_Accuracy.x,prs_img_cnb_factor_match$NAR_Overall_Accuracy.y)
###### !!!! match is successful !!!! ######

prs_img_match_key <- data.frame(org_bblid = prs_img_matched$bblid.y, new_bblid = prs_img_matched$bblid.x)
prs <- read.csv(paste0(local_wkdir,'genotype/pnc_EA_prs_20170501.csv'))


## Construct healthy Sample ##
  sample_demo <- read.csv(paste0(remote_wkdir,'n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv'))
  # apply subject-level exclusion
  hx_qa <- read.csv(paste0(remote_wkdir,"n1601_dataFreeze/health/n1601_health_20161214.csv"))
  sample_hx <- merge(sample_demo,hx_qa)
  sample_qa <- subset(sample_hx, healthExcludev2 == 0)
  
  # apply strc exclusion
  t1_qa <- read.csv(paste0(remote_wkdir,"n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv"))
  sample_t1 <- merge(sample_qa, t1_qa)
  sample_qa <- subset(sample_t1, t1Exclude == 0)

## Construct healthy Sample ##
  
  # load modality exclusion file from the data-freeze
  fc_qa <- read.csv(paste0(remote_wkdir,"n1601_dataFreeze/neuroimaging/rest/n1601_RestQAData_20170714.csv"))
  sample_fc <- merge(sample_qa,fc_qa)
  fc_sample_qa <- subset(sample_fc, restExclude ==0)  
  

## Construct SC Sample ##
  
  # load modality exclusion file from the data-freeze
  sc_qa <- read.csv(paste0(remote_wkdir,"n1601_dataFreeze/neuroimaging/dti/n1601_dti_qa_20170301.csv"))
  sample_sc <- merge(sample_qa,sc_qa)
  sc_sample_qa <- subset(sample_sc, dti64Exclude ==0)

## Construct FSC sample ##

fc_sc_sample_qa <- merge(fc_sample_qa,sc_sample_qa)

## add new bblid ##
fc_sc_sample_qa_newkey <- merge(fc_sc_sample_qa,prs_img_match_key, by.x = 'bblid', by.y='org_bblid' )
fc_sample_qa_newkey <- merge(fc_sample_qa,prs_img_match_key, by.x = 'bblid', by.y='org_bblid' )
sc_sample_qa_newkey <- merge(sc_sample_qa,prs_img_match_key, by.x = 'bblid', by.y='org_bblid' )

fc_sc_sample_qa_newkey <- merge(fc_sc_sample_qa_newkey,prs, by.x = 'new_bblid',by.y = 'bblid')
fc_sample_qa_newkey <- merge(fc_sample_qa_newkey,prs, by.x = 'new_bblid',by.y = 'bblid')
sc_sample_qa_newkey <- merge(sc_sample_qa_newkey,prs, by.x = 'new_bblid',by.y = 'bblid')

save(fc_sc_sample_qa_newkey,fc_sample_qa_newkey,sc_sample_qa_newkey, file = paste0(remote_wkdir,'../../projects/prsConnectivity/result/sample_qa_newkey.RData'))

### actual power analysis ##
parcellation = 'power'
sample = fc_sample_qa_newkey
resolution =  NA
modality = NA

sample_net <- get_net_from_sample(sample = sample,parcellation = parcellation,resolution = resolution, modality = modality)
community_info <-get_power_node_info(6)

## within and  between network power FC ##
sample_net_processed = lapply(sample_net, function(net) {diag(net)<-NA; net} )
power_prs_gams <- get_gams(community_info = community_info,sample_net = sample_net_processed)



has_default <- sapply(colnames(power_prs_gams$gam$pval),function(name) grepl('default',name))
has_default_pval <- power_prs_gams$gam$pval['scz_prs',has_default]

has_default_fdr <- p.adjust(power_prs_gams$gam$pval['scz_prs',has_default],method ='fdr')

within_coms <- diag(get_community_net_names(CommunityName =community_info$CommunityName ))
within_coms_pval <- power_prs_gams$gam$pval['scz_prs',within_coms]
within_coms_fdr <- p.adjust(power_prs_gams$gam$pval['scz_prs',within_coms],method ='fdr')

### actual scahefer fc 200 analysis ##
parcellation = 'schaefer'
sample = fc_sample_qa_newkey
resolution =  200
modality = 'fc'

fc_s200_sample_net <- get_net_from_sample(sample = sample,parcellation = parcellation,resolution = resolution, modality = modality)
community_info <-get_schafer_node_info(resolution = resolution,num_com = 17)

## within and  between network power FC ##
empty_net<-which(is.na(fc_s200_sample_net) == T)
non_empty_sample_net <- sc_sample_net
non_empty_sample_net[empty_net] <- NULL
non_empty_sample <- sc_sample_qa_newkey[-empty_net,]

fc_s200_prs_gams <- get_gams(community_info = community_info,sample_net = fc_s200_sample_net,sample = sample)
fc_prs_gams$gam$pval['scz_prs',which(sc_prs_gams$gam$pval['scz_prs',] <0.05)]

has_default <- sapply(colnames(fc_s200_prs_gams$gam$pval),function(name) grepl('default',name))
has_default_pval <- fc_s200_prs_gams$gam$pval['scz_prs',has_default]

has_default_fdr <- p.adjust(has_default_pval,method ='fdr')

within_coms <- diag(get_community_net_names(CommunityName =community_info$CommunityName ))
within_coms_pval <- fc_s200_prs_gams$gam$pval['scz_prs',within_coms]
within_coms_fdr <- p.adjust(within_coms_pval,method ='fdr')

### actual scahefer fc 400 analysis ##
parcellation = 'schaefer'
sample = fc_sample_qa_newkey
resolution =  400
modality = 'fc'

fc_s400_sample_net <- get_net_from_sample(sample = sample,parcellation = parcellation,resolution = resolution, modality = modality)
community_info <-get_schafer_node_info(resolution = resolution,num_com = 17)

## within and  between network power FC ##
empty_net<-which(is.na(fc_s200_sample_net) == T)
non_empty_sample_net <- sc_sample_net
non_empty_sample_net[empty_net] <- NULL
non_empty_sample <- sc_sample_qa_newkey[-empty_net,]

fc_s400_prs_gams <- get_gams(community_info = community_info,sample_net = fc_s400_sample_net,sample = sample)

has_default <- sapply(colnames(fc_s400_prs_gams$gam$pval),function(name) grepl('default',name))
has_default_pval <- fc_s400_prs_gams$gam$pval['scz_prs',has_default]

has_default_fdr <- p.adjust(has_default_pval,method ='fdr')

within_coms <- diag(get_community_net_names(CommunityName =community_info$CommunityName ))
within_coms_pval <- fc_s400_prs_gams$gam$pval['scz_prs',within_coms]
within_coms_fdr <- p.adjust(within_coms_pval,method ='fdr')




### actual structural analysis ##
parcellation = 'schaefer'
sample = sc_sample_qa_newkey
resolution =  200
modality = 'sc-d'

sc_sample_net <- get_net_from_sample(sample = sample,parcellation = parcellation,resolution = resolution, modality = modality)
community_info <-get_schafer_node_info(resolution = resolution,num_com = 17)

## within and  between network power FC ##
empty_net<-which(is.na(sc_sample_net) == T)
non_empty_sample_net <- sc_sample_net
non_empty_sample_net[empty_net] <- NULL
non_empty_sample <- sc_sample_qa_newkey[-empty_net,]

sc_sample_net_processed = lapply(sc_sample_net, function(net) {diag(net)<-NA; net} )
sc_prs_gams <- get_gams(community_info = community_info,sample_net = non_empty_sample_net,sample = non_empty_sample)
sc_prs_gams$gam$pval['scz_prs',which(sc_prs_gams$gam$pval['scz_prs',] <0.05)]


### gordon analysis ##
parcellation = 'gordon'
sample = fc_sample_qa_newkey
resolution =  NA
modality = NA

gordon_sample_net <- get_net_from_sample(sample = sample,parcellation = parcellation,resolution = resolution, modality = modality)
community_info <-get_schafer_node_info(resolution = resolution,num_com = 17)


