## set directories ##
local_wkdir <- '~/Google Drive/TDSlab/SCZ_gene_imaging/'
remote_wkdir <- '~/Desktop/BBL/data/joy/BBL/studies/pnc/'

load(paste0(remote_wkdir,'../../projects/prsConnectivity/result/sample_qa_newkey.RData'))

parcellation = 'power'
sample = fc_sample_qa_newkey
resolution =  NA
modality = NA

load(paste0(remote_wkdir,'../../projects/prsConnectivity/result/fc_sample_qa_newkey_power_network.RData'))

#compute bet/within gams#
community_info <-get_power_node_info(6)
sample_net_processed = lapply(sample_net, function(net) {diag(net)<-NA; net} )
power_prs_gams <- get_gams(community_info = community_info,sample_net = sample_net_processed,sample = sample)

#p-val correction#
has_default <- sapply(colnames(power_prs_gams$gam$pval),function(name) grepl('default',name))
has_default_pval <- power_prs_gams$gam$pval['scz_prs',has_default]
has_default_fdr <- p.adjust(has_default_pval,method ='fdr')

within_coms <- diag(get_community_net_names(CommunityName =community_info$CommunityName ))
within_coms_pval <- power_prs_gams$gam$pval['scz_prs',within_coms]
within_coms_fdr <- p.adjust(has_default_pval,method ='fdr')

#nodes and edges within DMN#
com_of_interest  = 'default'
get_coi_stats <-  function(net,com_of_interest){
  com_num_of_int <- which(community_info$CommunityName == com_of_interest)
  coi_nodes <- which(community_info$CommunityAffiliation == com_num_of_int)
  coi_net <- net[coi_nodes,coi_nodes]
  coi_net_graph <- graph_from_adjacency_matrix(coi_net,mode = 'undirected',diag = FALSE,weighted = TRUE)
  coi_net_strength <- graph.strength(coi_net_graph,loops = FALSE)
  coi_net_cent <- betweenness(coi_net_graph,weights = abs(E(coi_net_graph)$weight))
  coi_edge <- coi_net[upper.tri(coi_net)]
  edge_names_net <- get_community_net_names(colnames(coi_net))
  edge_names <- edge_names_net[upper.tri(edge_names_net)]
  out <- c(coi_net_strength,coi_net_cent,coi_edge)
  names(out) <-  c(paste0('stg',names(coi_net_strength)),paste0('stg',names(coi_net_cent)),edge_names)
  out
}

power_coi_stats <-  t(sapply(sample_net_processed,function(net) get_coi_stats(net, 'default')))
power_coi_stats_prs <- cbind(fc_sample_qa_newkey,power_coi_stats)
power_coi_stats_prs$sex <- as.ordered(as.factor(power_coi_stats_prs$sex))
power_coi_gams <- apply_prs_gam(power_coi_stats_prs,colnames(power_coi_stats_prs)[61:length(colnames(power_coi_stats_prs))])
