library(docopt)
library(DESeq2)

'Usage: run_deseq2.R <count_matrix>' -> doc
opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

#initialize matrix of counts 
counts <- read.table(opts$count_matrix, header=TRUE, row.names='gene_id')
counts <- as.matrix(counts)

#make coldata df
samplenames <- colnames(counts)
conditions <- c(rep('ctl', 4), rep('kd', 4))
conditions <- as.factor(conditions)
coldata <- data.frame(samplenames, conditions)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ conditions)
dds$conditions <- relevel(dds$conditions, ref='ctl')

#construct file names
day <- unlist(strsplit(basename(opts$count_matrix), '_'))[4]
deseqout <- paste('verse', day, 'deseq2', 'lfcthreshold', 'de', 'tsv', sep = '.')
rldname <- paste('verse', day, 'rld', 'RDS', sep = '.')

#DESeq analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c('conditions', 'kd', 'ctl'), lfcThreshold=0.5, altHypothesis='greaterAbs')
resOrdered <- res[order(res$padj),]
write.table(as.data.frame(resOrdered), file=file.path('/projectnb/perissilab/hmads_gps2_shrna', deseqout), sep='\t', quote=FALSE)  

#write out regularized log transform
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file=file.path('/projectnb/perissilab/hmads_gps2_shrna', rldname))

#extracting normalized counts
dds <- estimateSizeFactors(dds)
norm.counts <- counts(dds, normalized=TRUE)
name <- paste('norm', day, 'verse_counts', 'tsv', sep='.')
colnames(norm.counts) <- samplenames
write.table(norm.counts, file=file.path('/projectnb/perissilab/hmads_gps2_shrna', name), sep='\t', col.names=TRUE, quote=FALSE)
