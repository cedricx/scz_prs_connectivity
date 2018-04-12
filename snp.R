## set directories ##
local_wkdir <- '~/Google Drive/TDSlab/SCZ_gene_imaging/'
remote_wkdir <- '~/Desktop/BBL/data/joy/BBL/studies/pnc/'


### power CCA analysis ##
parcellation = 'power'
sample = fc_sample_qa_newkey
resolution =  NA
modality = NA

snp_dosage <- read.csv(paste0(local_wkdir,'genotype/illumina_EA_poly_scz108_converted.csv'))
num_snp <- length(colnames(snp_dosage))-1
colnames(snp_dosage)[2:(num_snp+1)] <- sapply(colnames(snp_dosage)[2:(num_snp+1)],function(name) gsub("X", "ch", name))

snp_dosage_go1<- merge(snp_dosage,fc_sample_qa_newkey,by.x='bblid',by.y='new_bblid')
