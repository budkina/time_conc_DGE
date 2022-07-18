library(DESeq2)
library(edgeR)
library(ggplot2)

#################### WRITE DESIGN TABLE MINUS ####################

time<-c(2,4,8,24,48)
conc<-c(0, 10, 100, 1000, 10000)
conc_str<-c("0_", "001_", "01_", "1_", "10_")
replicates <- c(1,2,3,4,5,6,7)
df<-expand.grid(replicates, time, conc_str)
zero_point <- expand.grid(replicates, c(0), c("0_"))
df <- rbind(df,zero_point)
colnames(df) <- c("replicate", "time", "conc")
filename <- paste0(df$conc, df$time,"h/", df$conc, df$time, "h_", df$replicate, "_featurecounts.txt")
file_table <- expand.grid(replicates, time, conc)
zero_point <- expand.grid(replicates, c(0), c(0))
file_table <- rbind(file_table,zero_point)
file_table$filename <- filename
colnames(file_table) <- c("replicate", "time", "conc", "filename")
write.table(file_table, "PMA_minus", row.names = F, quote = F)
write.table(file_table, "PMA_plus", row.names = F, quote = F)

#################### WRITE DESIGN TABLE PLUS ####################

time<-c(2,4,8,24,48)
conc<-c(0, 10, 100, 1000, 10000)
conc_str<-c("0_", "001_", "01_", "1_", "10_")
replicates <- c(1,2,3,4,5,6,7)
df<-expand.grid(replicates, conc_str, time)
zero_point <- expand.grid(replicates, c("0_"), c(0))
df <- rbind(df,zero_point)
colnames(df) <- c("replicate", "conc", "time")
filename <- paste0(df$conc, df$time,"h/", df$conc, df$time, "h_", df$replicate, "_featurecounts.txt")
file_table <- expand.grid(replicates, conc, time)
zero_point <- expand.grid(replicates, c(0), c(0))
file_table <- rbind(file_table,zero_point)
file_table$filename <- filename
colnames(file_table) <- c("replicate", "conc", "time",  "filename")
write.table(file_table, "PMA_minus_rev", row.names = F, quote = F)
write.table(file_table, "PMA_plus_rev", row.names = F, quote = F)

#################### set type ####################

type <- 'minus'

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

#################### Prepage covariates and count matrix ####################

group <- paste(file_table$time, file_table$conc, sep=".")
keep <- rowSums(time_conc) >= 50
time_conc_filtered <- time_conc[keep,]

#################### PCA plot ####################

condition <- factor(group, levels = c('0.0', '2.0', '2.10', '2.100', '2.1000', '2.10000', '4.0', '4.10', '4.100', '4.1000', '4.10000', '8.0', '8.10', '8.100', '8.1000', '8.10000', '24.0', '24.10', '24.100', '24.1000', '24.10000', '48.0', '48.10', '48.100', '48.1000', '48.10000'))
colData <- data.frame(row.names=colnames(time_conc), condition)
dds1 <- DESeqDataSetFromMatrix(countData = time_conc_filtered, colData = colData ,  design=~condition) 
dds <-DESeq(dds1, parallel = T)
rld <- vst(dds, blind=FALSE)
g <- plotPCA(rld, intgroup="condition")
g <- g + ggtitle(paste0(type, 'PMA PCA plot'))
ggsave(paste0(type,'PMA_PCA.png'),  plot = g,  device = "png", width = 30,  height = 20,  units = "cm")