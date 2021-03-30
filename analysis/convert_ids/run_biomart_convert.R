library(docopt)
library(biomaRt)

'Usage: run_biomart_convert.R <gene_list>' -> doc
opts <- docopt(doc, commandArgs(trailingOnly = TRUE))

convertgenes <- function(x){
	require('biomaRt')
	ensembl97mus <- useMart(host = 'jul2019.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
	ensembl97hum <- useMart(host = 'jul2019.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
	genes = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = x , mart = ensembl97hum, attributesL = c("ensembl_gene_id"), martL = ensembl97mus, uniqueRows=T) 
	return(genes)
}

convertsymbols <- function(y){
	require('biomaRt')
   	ensembl97mus <- useMart(host = 'jul2019.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl') 
	symbols = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = y , mart = ensembl97mus, attributesL = c("mgi_symbol"), martL = ensembl97mus, uniqueRows = T)	
	return(symbols)
}

gene_list <- readLines(opts$gene_list)
idmat <- convertgenes(gene_list)
symbolmat <- convertsymbols(idmat[,2])

method <- unlist(strsplit(basename(opts$gene_list), '.', fixed=TRUE))[[1]]
day <- unlist(strsplit(basename(opts$gene_list), '.', fixed=TRUE))[[2]]

idname <- paste(method, day, 'hum-id_to_mus-id', 'tsv', sep='.')
symbolname <- paste(method, day, 'mus-id_to_mgi', 'tsv', sep='.')

write.table(idmat, file = file.path('samples', idname), sep='\t', row.names=FALSE, quote=FALSE, col.names = c('hum_id', 'mus_id'))
write.table(symbolmat, file = file.path('samples', symbolname), sep='\t', row.names=FALSE, quote=FALSE, col.names= c('mus_id', 'mgi')) 
