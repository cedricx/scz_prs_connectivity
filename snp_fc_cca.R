## set directories ##
local_wkdir <- '~/Google Drive/TDSlab/SCZ_gene_imaging/'
remote_wkdir <- '~/Desktop/BBL/data/joy/BBL/studies/pnc/'


### power CCA analysis ##
parcellation- = 'power'
sample = fc_sample_qa_newkey
resolution =  NA
modality = NA

### load snp and merge with sample ##

snp_dosage <- read.csv(paste0(local_wkdir,'genotype/illumina_EA_poly_scz108_converted.csv'))
num_snp <- length(colnames(snp_dosage))-1
colnames(snp_dosage)[2:(num_snp+1)] <- sapply(colnames(snp_dosage)[2:(num_snp+1)],function(name) gsub("X", "chr", name))
snp_dosage_go1<- merge(snp_dosage,sample,by.x='bblid',by.y='new_bblid')

### pull sample net modual connectivity  ##
community_info <- get_power_node_info(6)
com_of_int <- c(5)

compute_com_edges<-function(net,coi){
    roi = which(community_info$CommunityAffiliation == coi)
    roi_net = net[roi,roi]
    roi_net_tri = roi_net[upper.tri(roi_net,diag = FALSE)]
}

sample_com_edges<-t(sapply(sample_net_processed,function(net) compute_com_edges(net,5)))


source(paste0(local_wkdir,'..//CEDRIC/sCCA/sCCA/code/final/cca_functions.R'))
require('PMA')
require('Matrix')
require('parallel')
require('emdbook')
require('caret')
require('R.matlab')
require('MASS')
require('permute')
require('matrixStats')
require('scales')
require('cowplot')
require('ggplot2')
require('ggrepel')
require('rasterVis')
x = snp_dosage_go1[,which(substr(colnames(snp_dosage_go1),start = 1,stop = 3) == 'chr')]
z = sample_com_edges

sampleid <- createDataPartition(sample$scz_prs, p = 0.667, list =T,times=10)

x_sample <- mclapply(sampleid, function(id) x[id,])
z_sample <- mclapply(sampleid, function(id) z[id,])






x_pen <- seq(0.1,1,length.out=10)
y_pen <- seq(0.1,1,length.out=10)

scca.gs<-ccaDWfoldgs(x_sample,z_sample,x_pen,y_pen)
gs.mat <- matrix(scca.gs$GS[,'COR_MEAN'], nrow = 10, ncol = 10)
rownames(gs.mat) <- x_pen
colnames(gs.mat) <- y_pen
image(gs.mat)


modenum <- dim(x)[2] #number of all possible canonical variates
scca.org <- ccaDW(x, z,1,0.2,modenum) #0.8 and 0.4 here are the best parameteres selected above in the grid search
x_std <- apply(x,2,scale) #make sure data are demeaned
z_std <- apply(z,2,scale)
covmat <- t(scca.org$u) %*% t(x_std) %*% z_std %*% scca.org$v #calculate covariance matrix
varE <- diag(covmat)^2 / sum(diag(covmat)^2) #calcualte covariance explained by each component
varE.df <- data.frame(modenum = as.factor(1:modenum), var = varE) #prepare dataframe for plotting
candnum = 2 #number selected based on the scree plot

p.var<-ggplot(varE.df,aes(modenum,var)) +
  geom_point(stat = 'identity',aes(color = var > varE[candnum+1], size = var)) +
  geom_hline(yintercept = 1/modenum,linetype="dashed") +
  scale_x_discrete(name ="Mode", limits=c(0:modenum),breaks =  c(1,seq(10,modenum,10))) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.075),labels = percent,name = "Variance Explained", breaks=seq(0,0.075,length=4)) +
  theme_classic(base_size = 20) +
  theme(legend.position = 'none') 

p.var



scca.cand <- ccaDW(x, z,1,0.2,candnum)
scca.cca<-mclapply(seq_along(sampleid),function(i) ccaDW(x_sample[[i]],z_sample[[i]],1,0.2,5)) #loop through split
scca.cca.ro <- sapply(scca.cca,function(x) reorderCCA(x,scca.cand,5)) #reorder the component so they match across splits
scca.cca.cor <- rowMeans(simplify2array(scca.cca.ro['cors',]),na.rm =T) #calculate average of cca correlations
scca.cca.cor.se <- rowSds(simplify2array(scca.cca.ro['cors',]),na.rm =T)/sqrt(dim(scca.cca.ro)[2]) #calculate standard error of correlations

cor.df <- data.frame(modenum = as.factor(1:candnum), cor = scca.cca.cor, se = scca.cca.cor.se)
cor.df.order <- cor.df[order(-cor.df$cor),]
cor.lim <- aes(ymax = cor.df.order$cor + cor.df$se, ymin = cor.df.order$cor - cor.df$se)
p.cor <- ggplot(cor.df.order,aes(1:length(modenum), cor, label = round(cor,2))) +
  geom_bar(width = 0.75, stat = 'identity',  fill = '#00BFA5') +
  geom_errorbar(cor.lim,  width=0.25) +
  geom_text(size = 4, position = position_dodge(width = 0.9), vjust= -1,color='grey')+
  scale_x_discrete(name ="Mode", limits=c(1:candnum) ) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.8),name = "CCA Correlation", breaks=seq(0,0.8,length=5)) +
  theme_classic(base_size = 20) +
  coord_cartesian(ylim=c(0.2,0.9)) +
  theme(legend.position = 'none') 
p.cor


num.perm <- 2000 #number of permutaitons to run
z.perm <- rlply(num.perm,z[sample(nrow(z)),]) #permute the clinical matrix by row
scca.perm.cca<-sapply(z.perm, function(z_perm){ out<-ccaDWpermorder(x,z_perm,1,0.2,candnum,scca.cand)} ) #run scca again with permuted clinical but with original connectivity
#load("~/Desktop/BBL/projects/xiaNetworkCca/sCCA/aim1/result/201701/pwr_perm_cca.RData")
perm.cor <- simplify2array(scca.perm.cca['cors',]) #extrac the correlations
perm.pval <- sapply(seq_along(cor.df$cor),function(x) (length(which(perm.cor[x,] >= cor.df$cor[x])) ) / length(which(is.na(perm.cor[x,]) == FALSE))) #calcualte the empirical p-val

perm.cor.df <- as.data.frame(t(perm.cor))
perm.pass <- which(perm.pval < 0.05)
permplots <-lapply(perm.pass,function(x) perm.plot(perm.cor.df,cor.df,perm.pval,x))


bootnum = 1000 #number of bootstraps to run
bootid<-createResample(sample$scz_prs, list = T, times = bootnum) #create lists of subjs for bootstrap
x_boot <- lapply(bootid, function(id) x[id,]) #creat bootstrap samples for connectivity features
z_boot <- lapply(bootid, function(id) z[id,]) #creat bootstrap samples for clinical features

scca.boot<- mclapply(seq_along(bootid),function(i) ccaDW(x_boot[[i]],z_boot[[i]],1,0.2,5),mc.cores = 5) #run scca on these bootstrap sample

scca.boot.ro<- lapply(1:bootnum,function(i) reorderCCA(scca.boot[[i]],scca.cand,5)) #reorder to match components across samples
scca.boot.u <- lapply(perm.pass, function(x) sapply(1:bootnum, function(i) scca.boot.ro[[i]]$u[,x])) #extract loadings on connectivity features
scca.boot.v <- lapply(perm.pass, function(x) sapply(1:bootnum, function(i) scca.boot.ro[[i]]$v[,x])) #extract loadings on clinical features
#scca.boot.cor <-  sapply(1:1000, function(i) scca.boot.ro[[i]]$cor)
u.boot.plot <- lapply(seq_along(perm.pass), function(x) bootplot_u(scca.cand$u[,perm.pass[x]], scca.boot.u[[x]] ))  #apply 99.5% confidence interval
v.boot.plot <- lapply(seq_along(perm.pass), function(x) bootplot(scca.cand$v[,perm.pass[x]], scca.boot.v[[x]] )) #apply 95% confidence interval

snp_info_pnc <- read.table(paste0(local_wkdir,'genotype/scz2_103_rsids_ref_alt_info.txt'),header = T)
snp_info108 <- read.table(paste0(local_wkdir,'genotype/UNC_PGC/scz2.regions/scz2.anneal.108.txt'),header = T)
snp_info128 <- read.table(paste0(local_wkdir,'genotype/UNC_PGC/scz2.regions/scz2.rep.128.txt'),header = T)
snp_info_merge <- merge(snp_info108, snp_info128, by.x ='bestsnp',by.y='snpid')

dim1_snp <- colnames(x)[u.boot.plot[[1]]$fea]
dim1_snp_locus <- sapply(dim1_snp, function(snp) unlist(strsplit(snp,'[.]')))

match_pnc_snp <- function(pnc_snp) {
  pnc_snp_chr <- as.numeric(substr(pnc_snp[1],start = 3, stop = nchar(pnc_snp[1])))
  pnc_snp_locus <- as.numeric(pnc_snp[2])
  snp_info_merge_chr <- snp_info_merge[which(snp_info_merge$chr == pnc_snp_chr),]
  start_locus <-which(snp_info_merge_chr$anneal1 < pnc_snp_locus)
  stop_locus <-which(snp_info_merge_chr$anneal2 > pnc_snp_locus)
  out <- snp_info_merge_chr[start_locus[start_locus %in% stop_locus],]
  out
}

dim1_snp_locus_pgc  <- apply(dim1_snp_locus,2,match_pnc_snp)


