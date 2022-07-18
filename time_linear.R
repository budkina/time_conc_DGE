library(dplyr)
library(edgeR)
library(DESeq2)
library(DEGreport)
library(BiocParallel)
library(ggplot2)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 16, progressbar = TRUE))

type = 'minus'

#################### read gene id table ####################

id_to_gene <- read.table('enseble_annotation_genes_result.csv', sep = '\t')
rownames(id_to_gene) <- id_to_gene$V1
na_rownames <- rownames(id_to_gene[id_to_gene$V2=="",])
id_to_gene[na_rownames,'V2']<- na_rownames

#################### read count matrix ####################

file_table <- read.table(paste0("PMA_",type), header = T)
time_conc_matrix <-c()
for (i in file_table$filename)
{
		r <- read.table(i, header = T)
		time_conc_matrix <- c(time_conc_matrix, r[7])
}

time_conc <- data.frame(time_conc_matrix)
rownames(time_conc) <- r[,1]

#################### Prepare covariates and count matrix ####################

group <- paste(file_table$time, file_table$conc, sep=".")
keep <- rowSums(time_conc) >= 50
time_conc_filtered <- time_conc[keep,]

#################### Clusterization function ####################

clustering <- function(ma, metadata, id_to_gene, is_top, type, method, conc)
{
	tma <-t(ma)
	tma <- tma[, ! apply(tma, 2 , function(x) sd(x, na.rm = TRUE)==0 ) ]
	ma <- t(tma)
	suffix <- ""
	if(is_top)
	{
		ma <- head(ma, 2000)
		suffix <- "_2000"
	}

	res <- degPatterns(ma, metadata, time = "condition")
	ggsave(paste0(type,'PMA_', method, '_time_', conc, '_linear', suffix, '.png'),  plot = res$plot,  device = "png", width = 30,  height = 20,  units = "cm")
	save(res, file = paste0(type, 'PMA_', method, '_time_', conc, '_linear', suffix, '.RData'))
	degPlotCluster(res$normalized, time = "condition")
	dev.off()
	file.rename('Rplots.pdf', paste0(type,'PMA_', method, '_time_', conc, '_linear', suffix, '.pdf'))
	clusters <- res$df 
	rownames(clusters) <- clusters$genes
	clusters$geneID <- clusters$genes
	clusters$geneName <- id_to_gene[rownames(clusters), ]$V2
	clusters$genes <- NULL
	clusters <- clusters[,c(3,2,1)]
	write.table(clusters, paste0(type, 'PMA_', method, '_time_', conc, '_clusters_linear', suffix, '.csv'), sep = '\t',  quote = F,  row.names = F)

	#################### Save plots to separate pictures ####################

	cluster_ids <- unique(clusters$cluster)
	for (cluster in cluster_ids)
	{
		subset_genes <- clusters[clusters$cluster == cluster,]$geneID
		ma_subset <- ma[subset_genes,]
		res <- degPatterns(ma_subset, metadata, time = "condition", minc = 2)
		ggsave(paste0("const_conc_sep_cluster_pictures/", type, 'PMA_', method, '_time_', conc, '_linear_', cluster, '.png'),
			plot = res$plot + ggtitle(paste0(type,'PMA, ', method, ', linear, conc = ', conc, ' cluster ', cluster)),
			device = "png", width = 30, height = 20, units = "cm")
	}
}

#################### Linear model ####################

concs <- c('0','001','01','1','10')
concs_num <- c('0','10','100','1000','10000')

de_num <- data.frame()
x <- c("method", "conc", "de_num")
colnames(de_num) <- x

for (i in 1:5)
{
	conc = concs[i]
	conc_num = concs_num[i]
	time_conc_filtered_case <- select(time_conc_filtered, starts_with(paste0('X',conc,'_')))
	group <- file_table[file_table$conc == conc_num,]$time

	m<- poly(group, df=1)
	design <- model.matrix(~group)

	#################### DESeq2 LRT ####################

	colnames(m)<-c('m1')
	condition <- factor(group)
	colData <- data.frame(row.names=colnames(time_conc_filtered_case), m)
	colData <- cbind(colData, condition)
	dds1 <- DESeqDataSetFromMatrix(countData = time_conc_filtered_case, colData = colData , design = ~m1) 
	dds <-DESeq(dds1,  test="LRT",  reduced= ~1, parallel = T)
	top_table <- results(dds)
	top_table <- top_table[order(top_table$padj),]
	newdata <- top_table[which(top_table$padj < 0.05),]
	newdata$geneID <- rownames(newdata)
	newdata$geneName <- id_to_gene[rownames(newdata), ]$V2
	newdata <- newdata[,c(8,7,1,2,3,4,5,6)]
	write.table(newdata, paste0(type,'PMA_DESeq2_time_',conc,'_linear.csv'), sep = '\t',  quote = F,  row.names = F)
	print(dim(newdata))
	de_num <- rbind(de_num, c('DESeq2', conc, nrow(newdata)))

	#################### clustering ####################

	ma = assay(rlog(dds))[row.names(newdata),]
	ma <- ma[ rowSums(ma)!=0, ] 
	metadata <- as.data.frame(colData(dds))
	clustering(ma, metadata, id_to_gene, TRUE, type, 'DESeq2', conc)
	clustering(ma, metadata, id_to_gene, FALSE, type, 'DESeq2', conc)

	#################### edgeR ####################

	dge_list <- DGEList(counts=time_conc_filtered_case,  group = group)
	dge_list_norm <- calcNormFactors(dge_list)
	y <- estimateDisp(dge_list_norm, design)
	fit <- glmQLFit(y, design, robust=TRUE)
	fit <- glmQLFTest(fit, coef=2)
	top_tags <- topTags(fit, n = nrow( fit$table ), sort.by = "PValue", p.value = 0.05)$table
	result_edger <- top_tags[which(top_tags$FDR < 0.05),]
	newdata <- result_edger[order(-abs(result_edger$logFC)),]
	newdata$geneID <- rownames(newdata)
	newdata$geneName <- id_to_gene[rownames(newdata), ]$V2
	newdata <- newdata[,c(7,6,1,2,3,4,5)]
	write.table(newdata, paste0(type,'PMA_edgeR_time_',conc,'_linear.csv'), sep = '\t',  quote = F,  row.names = F)
	print(dim(newdata))
	de_num <- rbind(de_num, c('edgeR', conc, nrow(newdata)))

	#################### clustering ####################

	ma = cpm(y, log=TRUE)[row.names(newdata),]
	ma <- ma[ rowSums(ma)!=0, ]
	metadata <- as.data.frame(colData(dds))
	clustering(ma, metadata, id_to_gene, TRUE, type, 'edgeR', conc)
	clustering(ma, metadata, id_to_gene, FALSE, type, 'edgeR', conc)

	#################### limma + voom ####################

	dge_voom <- voom(dge_list_norm, design)
	fit <- lmFit(dge_voom, design)
	fit <- eBayes(fit)
	top_table <- topTable(fit, coef=2, n = Inf)
	result_voom <-top_table[which(top_table$adj.P.Val < 0.05),]
	newdata <- result_voom[order(-abs(result_voom$logFC)),]
	newdata$geneID <- rownames(newdata)
	newdata$geneName <- id_to_gene[rownames(newdata), ]$V2
	newdata <- newdata[,c(8,7,1,2,3,4,5,6)]
	write.table(newdata, paste0(type,'PMA_voom_time_',conc,'_linear.csv'), sep = '\t',  quote = F,  row.names = F)
	print(dim(newdata))
	de_num <- rbind(de_num, c('voom', conc, nrow(newdata)))

	#################### clustering ####################

	ma = cpm(y, log=TRUE)[row.names(newdata),]
	ma <- ma[ rowSums(ma)!=0, ] 
	metadata <- as.data.frame(colData(dds))
	clustering(ma, metadata, id_to_gene, TRUE, type, 'voom', conc)
	clustering(ma, metadata, id_to_gene, FALSE, type, 'voom', conc)
}

de_num_filename <- paste0(type,'PMA_time_linear_linear_de_num.csv')
colnames(de_num) <- c("method", "case", "de_num")
write.table(de_num, de_num_filename, sep = '\t',  quote = F,  row.names = F)
