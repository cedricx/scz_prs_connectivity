require('stringr')
require('R.matlab')
require('mgcv')
require('visreg')
##############################
####                      ####
#### Schaefer Specifics   ####
####                      ####
##############################
get_schaefer_netpath <- function(resolution, modality, bblid, scanid){
  if (modality == 'fc' & (resolution == 200 | resolution == 100)) {
    netpath <- paste0(remote_wkdir,"processedData/restbold/restbold_201607151621/",bblid,"/*",scanid,"/net/Schaefer",resolution,"PNC/*_network.txt")
  } else if (modality == 'fc' & resolution == 400) {
    netpath <- paste0(remote_wkdir,"processedData/restbold/restbold_201607151621/",bblid,"/*",scanid,"/net/SchaeferPNC/*_network.txt")
  } else if (modality == 'sc-d') {
    netpath <- paste0(remote_wkdir,"processedData/diffusion/deterministic_20171118/",bblid,"/*",scanid,"/tractography/connectivity/*_",resolution,"_*_fa_connectivity.mat")
  } else if (modality == 'sc-p' & resolution == 400) {
    netpath <- paste0(remote_wkdir,"processedData/diffusion/probabilistic_20171118/",bblid,"/*",scanid,"/output/connectivity/*Schaefer",resolution,"_17net.mat")
  } else if (modality == 'sc-p' & resolution == 200) {
    netpath <- paste0(remote_wkdir,"processedData/diffusion/probabilistic_20171118/",bblid,"/*",scanid,"/output/schaefer",resolution,"/connectivity/*Schaefer",resolution,"_17net.mat")
  }
  netpath
}

get_schafer_node_info <- function(resolution,num_com){
  CommunityAffiliation <- read.csv2(paste0(local_wkdir,'imaging/schaefer',resolution,'/schaefer',resolution,'x',num_com,'CommunityAffiliation.1D'),header = FALSE)$V1
  CommunityName <- read.csv2(paste0(local_wkdir,'imaging/schaefer',resolution,'/schaefer',resolution,'x',num_com,'CommunityNames.txt'),header = FALSE)$V1
  NodeName <- read.csv2(paste0(local_wkdir,'imaging/schaefer',resolution,'/schaefer',resolution,'NodeNames.txt'),header = FALSE)
  out <- list(CommunityAffiliation = CommunityAffiliation, CommunityName = CommunityName, NodeName= NodeName)
}

##############################
####                      ####
####   Power Specific     ####
####                      ####
##############################

get_power_netpath <- function(scanid) {
  netpath <- paste0(remote_wkdir, 'n1601_dataFreeze/neuroimaging/rest/restNetwork_264PowerPNC/264PowerPNCNetworks/',scanid,'*.txt')
}

get_power_node_info <- function(num_com){
  CommunityAffiliation <- read.csv2(paste0(local_wkdir,'imaging/power264/power264CommunityAffiliation.1D'),header = FALSE)$V1
  CommunityName <- read.csv2(paste0(local_wkdir,'imaging/power264/power264CommunityNames.txt'),header = FALSE)$V1
  NodeName <- read.csv2(paste0(local_wkdir,'imaging/power264/power264NodeNames.txt'),header = FALSE)$V1
  NodeCoord <- read.delim(paste0(local_wkdir,'imaging/power264/power264CoorMNI.sclib'),header = F,dec = ',',skip = 2)$V1
  NodeCoord_df <- t(sapply(NodeCoord,function(coord) as.numeric(unlist(strsplit(as.character(coord),"\\,|\\#"))[3:5])))
  colnames(NodeCoord_df) <-c('X','Y','Z')
  if (num_com == 14){
  out <- list(CommunityAffiliation = CommunityAffiliation, CommunityName = CommunityName, NodeName= NodeName, NodeCoord = NodeCoord_df)
  }
  else if (num_com == 6){
    CommunityName2 <- CommunityName
    CommunityAffiliation2 <- array(NA,length(CommunityAffiliation))
    CommunityName2 <- c('somatomotor_aud','COP_SAL_vATT','FrontoParietal','vis','default','dATT')
    CommunityAffiliation2[which(CommunityAffiliation == 1 | CommunityAffiliation == 2 | CommunityAffiliation == 4)] <- 1  # merge somato
    CommunityAffiliation2[which(CommunityAffiliation == 3 | CommunityAffiliation == 9 | CommunityAffiliation == 11)] <- 2  # merge COP_salience_vatt
    CommunityAffiliation2[which(CommunityAffiliation == 8)] <- 3  # fronto
    CommunityAffiliation2[which(CommunityAffiliation  == 7 )] <- 4  # vis
    CommunityAffiliation2[which(CommunityAffiliation == 5 | CommunityAffiliation == 6)] <- 5  # default
    CommunityAffiliation2[which(CommunityAffiliation == 12)] <- 6  # datt
    out <- list(CommunityAffiliation = CommunityAffiliation2, CommunityName = CommunityName2, NodeName= NodeName, NodeCoord = NodeCoord_df)
  }
  else if (num_com == 3){
    CommunityName2 <- CommunityName
    CommunityAffiliation2 <- array(NA,length(CommunityAffiliation))
    CommunityName2 <- c('somatomotor','cogcontrol','default')
    CommunityAffiliation2[which(CommunityAffiliation == 1 | CommunityAffiliation == 2 | CommunityAffiliation == 4  | CommunityAffiliation == 7)] <- 1  # somatomotor
    CommunityAffiliation2[which(CommunityAffiliation == 3 | CommunityAffiliation == 8 | CommunityAffiliation == 9 | CommunityAffiliation == 11  | CommunityAffiliation == 12)] <- 2  #cogcontrol 
    CommunityAffiliation2[which(CommunityAffiliation == 5 | CommunityAffiliation == 6)] <- 3  # default
    out <- list(CommunityAffiliation = CommunityAffiliation2, CommunityName = CommunityName2, NodeName= NodeName, NodeCoord = NodeCoord_df)
  }
  out
}


##############################
####                      ####
####  Gordon Specific     ####
####                      ####
##############################
get_gordon_netpath <- function(scanid) {
  netpath <- paste0(remote_wkdir, 'n1601_dataFreeze/neuroimaging/rest/restNetwork_gordon/GordonPNCNetworks/',scanid,'*.txt')
}

get_gordon_node_info <- function(resolution,num_com){
  CommunityAffiliation <- read.csv2(paste0(local_wkdir,'imaging/gordon333/gordon333CommunityAffiliation.1D'),header = FALSE)$V1
  CommunityName <- read.csv2(paste0(local_wkdir,'imaging/gordon333/gordon333CommunityNames.txt'),header = FALSE)$V1
  NodeName <- read.csv2(paste0(local_wkdir,'imaging/gordon333/gordon333NodeNames.txt'),header = FALSE)
  out <- list(CommunityAffiliation = CommunityAffiliation, CommunityName = CommunityName, NodeName= NodeName)
}

##############################
####                      ####
#### grab network matrix  ####
####                      ####
##############################

grab_net_from_path <- function(netpath){
  net_file_type = str_sub(netpath,-3)
  if (identical(Sys.glob(netpath), character(0))) {
    temp_net <- NA
  } else  {
    if (net_file_type == 'txt') {  
      temp_net <- as.matrix(read.table(Sys.glob(netpath)))
    } else if (net_file_type == 'mat') {
      temp_net <- readMat(Sys.glob(netpath))
    }
  }
}

get_net_from_sample <- function(sample,parcellation,resolution,modality) {
  n_sample <- dim(sample)[1]
  sample_net<-list()
  for (i in 1:n_sample) {
    bblid = sample[i,'bblid']
    scanid = sample[i,'scanid']
    # set up the correct path by modality, and resolution#
    
    
    # import the network data #
    
    if (parcellation == 'power') {
      print(paste0(i,"/",n_sample,": copying ",bblid,'_',scanid, ' of ',parcellation,' atlas'))
      netpath <- get_power_netpath(scanid)
      sample_net[[i]] <- grab_net_from_path(netpath)
    } else if (parcellation == 'gordon') {
      print(paste0(i,"/",n_sample,": copying ",bblid,'_',scanid, ' of ',parcellation,' atlas'))
      netpath <- get_gordon_netpath(scanid)
      sample_net[[i]] <- grab_net_from_path(netpath)
    } else if (parcellation == 'schaefer') {
      netpath <- get_schaefer_netpath(resolution,modality,bblid,scanid)
      print(paste0(i,"/",n_sample,": copying ",bblid,'_',scanid, ' of resolution ', resolution,' in ', modality,' of ', parcellation,' atlas'))
      if (modality == 'sc-d') {
        temp_net <- grab_net_from_path(netpath)
        if (is.na(temp_net)) {
        sample_net[[i]] <- NA
        } else {
        sample_net[[i]] <- temp_net$connectivity
        }
      } else if  (modality == 'sc-p') {
        temp_net <- grab_net_from_path(netpath)
        if (is.na(temp_net)) {
          sample_net[[i]] <- NA
        } else {
          sample_net[[i]] <- temp_net$streamlineCount.mat
        }
      } else if (modality == 'fc') {
        sample_net[[i]] <- grab_net_from_path(netpath)
      }
    }
    
  }
  
  if (parcellation == 'power'){
  save_file_path <- paste0(remote_wkdir,'../../projects/prsConnectivity/result/',deparse(substitute(sample)),'_power_network.RData')
  }
  else if (parcellation == 'schaefer'){
  save_file_path <- paste0(remote_wkdir,'../../projects/prsConnectivity/result/',deparse(substitute(sample)),'_',modality,'_',resolution,'_schaefer_network.RData')
  } else if (parcellation == 'gordon'){
    save_file_path <- paste0(remote_wkdir,'../../projects/prsConnectivity/result/',deparse(substitute(sample)),'_gordon_network.RData')
  }
  save(sample_net,sample,file = save_file_path)
  sample_net
}


##############################
####                      ####
## calculate bet/within net ##
####                      ####
##############################


get_community_net_names <- function(CommunityName){
  community_net_names <-  array(NA, dim=c(length(CommunityName),length(CommunityName)))
  for (community_i in 1:length(CommunityName)){
    for (community_j in 1:length(CommunityName)) {
      community_net_names[community_i,community_j] = paste(CommunityName[community_i],CommunityName[community_j],sep = '_')
    }
  }
  community_net_names
}

get_community_net<-function(netmat,com_name_vector,node_aff_vector){
  community_net =  array(NA, dim=c(length(com_name_vector),length(com_name_vector)))
  for (community_i in 1:length(com_name_vector)){
    for (community_j in 1:length(com_name_vector)) {
      nodes_in_community_i = which(node_aff_vector == community_i)
      nodes_in_community_j = which(node_aff_vector == community_j)
      community_net[community_i,community_j] = 
        mean(netmat[nodes_in_community_i,nodes_in_community_j],na.rm=TRUE)
    }
  }
  colnames(community_net) <- com_name_vector
  rownames(community_net) <- com_name_vector
  community_net
}


##############################
####                      ####
####       PRS GAMs       ####
####                      ####
##############################
prs_img_gam10<-function(vab_of_int, prs_df){
  out<-gam(vab_of_int ~ s(ageAtScan1) +  sex + scz_prs + restRelMeanRMSMotion +
             pc1 + pc2 + pc3 + pc4 + pc5 +
             pc6 + pc7 + pc8 + pc9 + pc10, data = prs_df, method="REML")
}

prs_gam<-function(vab_of_int, prs_df){
  out<-gam(vab_of_int ~ s(ageAtScan1, k =4) + sex +  scz_prs , data = prs_df, method="REML")
}

apply_prs_gam <-function(prs_voi_df,voi_names) {
  prs_voi_reml<-lapply(voi_names,function(vab_of_int) prs_img_gam10(prs_voi_df[,vab_of_int],prs_voi_df))
  prs_voi_pval<-sapply(prs_voi_reml, function(reml_result) summary(reml_result)$p.pv)
  colnames(prs_voi_pval) <- voi_names
  voi_of_sig <- prs_voi_pval[,which(prs_voi_pval['scz_prs',] <0.05)]
  out <- list(data = prs_voi_df, reml_result = prs_voi_reml,pval = prs_voi_pval, sig_voi = voi_of_sig)
}

get_gams <- function(community_info,sample_net, sample){
  com_net_names <- get_community_net_names(community_info$CommunityName)
  com_net<-lapply(sample_net,function(netmat) get_community_net(netmat,community_info$CommunityName,community_info$CommunityAffiliation))
  com_net_flat<-t(sapply(com_net,function(net) net[upper.tri(net,diag = TRUE)]))
  colnames(com_net_flat) <- com_net_names[upper.tri(com_net_names,diag = TRUE)]
  com_net_names_flat <- com_net_names[upper.tri(com_net_names,diag = TRUE)]
  com_net_flat <- cbind(com_net_flat,sample)
  com_net_flat$sex <- as.ordered(as.factor(com_net_flat$sex))
  prs_gams <- apply_prs_gam(com_net_flat,com_net_names_flat)
  out <- list(gam = prs_gams, data = com_net_flat)
}

##############################
####                      ####
#### grab community edges ####
####                      ####
##############################
