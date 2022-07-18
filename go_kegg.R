
library(clusterProfiler)
library("org.Hs.eg.db")

library(biomaRt)
# Retrieve the ensembl info
ensembl <- useMart("ENSEMBL_MART_ENSEMBL") # host="asia.ensembl.org"
dsets = listDatasets(ensembl)
bm <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

pma_type = "minus"
method = "DESeq2"
regr_type = 'poly'

#################### conc = const ####################

conds = c('0', '001', '01', '1', '10')
for (cond in conds)
{
  time_0 = read.table(paste0(pma_type,"PMA_",method,"_time_",cond,"_clusters",".csv"), sep = '\t', header= T)
  
  clusters = unique(time_0['cluster'])$cluster
  for (cluster in clusters)
  {
    genes <- time_0[time_0['cluster'] == cluster,]$geneID
    
    ego_MF <- enrichGO(gene         = genes,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05)
    
    write.table(data.frame(ego_MF),
                paste0(pma_type,"PMA_",method,"_time_",cond,"_cluster",cluster,"_GO_MF_",regr_type,".csv"),
                sep = '\t',  quote = F,  row.names = F)
    
    ego_BP <- enrichGO(gene         = genes,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05)
    
    write.table(data.frame(ego_BP),
                paste0(pma_type,"PMA_",method,"_time_",cond,"_cluster",cluster,"_GO_BP_",regr_type,".csv"),
                sep = '\t',  quote = F,  row.names = F)
    
    ego_CC <- enrichGO(gene         = genes,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05)
    
    write.table(data.frame(ego_CC),
                paste0(pma_type,"PMA_",method,"_time_",cond,"_cluster",cluster,"_GO_CC_",regr_type,".csv"),
                sep = '\t',  quote = F,  row.names = F)
    
    # Filter the ensembl ID
    ids <- getBM(filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id", 'entrezgene_id', 'entrezgene_accession'), 
                 values = genes, mart = bm)
    
    kk <- enrichKEGG(gene         = ids$entrezgene_id,
                     organism     = 'hsa',
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05)
    
    write.table(data.frame(kk),
                paste0(pma_type,"PMA_",method,"_time_",cond,"_cluster",cluster,"_KEGG_",regr_type,".csv"),
                sep = '\t',  quote = F,  row.names = F)
  }
}

#################### time = const ####################

conds = c('2', '4', '8', '24', '48')
for (cond in conds)
{
  time_0 = read.table(paste0(pma_type,"PMA_",method,"_conc_",cond,"_clusters_",regr_type,".csv"), sep = '\t', header= T)
  
  clusters = unique(time_0['cluster'])$cluster
  for (cluster in clusters)
  {
    genes <- time_0[time_0['cluster'] == cluster,]$geneID
    
    ego_MF <- enrichGO(gene         = genes,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05)
    
    write.table(data.frame(ego_MF),
                paste0(pma_type,"PMA_",method,"_conc_",cond,"_cluster",cluster,"_GO_MF_",regr_type,".csv"),
                sep = '\t',  quote = F,  row.names = F)
    
    ego_BP <- enrichGO(gene         = genes,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05)
    
    write.table(data.frame(ego_BP),
                paste0(pma_type,"PMA_",method,"_conc_",cond,"_cluster",cluster,"_GO_BP_",regr_type,".csv"),
                sep = '\t',  quote = F,  row.names = F)
    
    ego_CC <- enrichGO(gene         = genes,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05)
    
    write.table(data.frame(ego_CC),
                paste0(pma_type,"PMA_",method,"_conc_",cond,"_cluster",cluster,"_GO_CC_",regr_type,".csv"),
                sep = '\t',  quote = F,  row.names = F)
    
    # Filter the ensembl ID
    ids <- getBM(filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id", 'entrezgene_id', 'entrezgene_accession'), 
                 values = genes, mart = bm)
    
    kk <- enrichKEGG(gene         = ids$entrezgene_id,
                     organism     = 'hsa',
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05)
    
    write.table(data.frame(kk),
                paste0(pma_type,"PMA_",method,"_conc_",cond,"_cluster",cluster,"_KEGG_",regr_type,".csv"),
                sep = '\t',  quote = F,  row.names = F)
  }
}