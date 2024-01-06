setwd("C:\\Users\\Dell\\Desktop\\")
library(dplyr)
library(data.table)
library(ggplot2)
library(limma)
library(DESeq2)

#UCSC Xena dealing for LUAD
###mRNA dealing
LUAD_mrna<-data.table::fread(".\\TCGA-LUAD.htseq_counts.tsv.gz")
LUAD_gene_map<-data.table::fread(".\\gencode.v22.annotation.gene.probeMap")

TCGA_rawdata<-inner_join(LUAD_mrna,LUAD_gene_map,by=c("Ensembl_ID" = "id")) %>% select(.,gene,starts_with('TCGA'))

TCGA_rawdata<-data.table::setDF(TCGA_rawdata)
TCGA_rawdata[,2:ncol(TCGA_rawdata)]<-apply(TCGA_rawdata[,2:ncol(TCGA_rawdata)],2,function(x){2^x-1})

TCGA_gene<-avereps(TCGA_rawdata[,-1],ID = TCGA_rawdata$gene) %>% as.data.frame()
table(duplicated(colnames(TCGA_gene)))
table(duplicated(rownames(TCGA_gene)))
TCGA_gene$gene_name<-rownames(TCGA_gene)
TCGA_gene<-TCGA_gene[,c(586,1:585)]
#reference gene sets
ref_set<-rtracklayer::import("Homo_sapiens.GRCh38.107.gtf.gz") %>% as.data.frame()
mRNA_exprSet <- ref_set %>%
  dplyr::filter(type=='gene',gene_biotype=='protein_coding') %>% 
  dplyr::select(c(gene_name,gene_id,gene_biotype)) %>%
  dplyr::inner_join(TCGA_gene,by ='gene_name') 


mRNAexp<-mRNA_exprSet[,-c(2,3)]
metadata <- data.frame(TCGA_id =colnames(mRNAexp)[-1])
table(substring(metadata$TCGA_id,14,15))
mRNAexp<-mRNAexp[,-which(substring(colnames(mRNAexp),14,15)=="02")]

metadata1 <- data.frame(TCGA_id =colnames(mRNAexp)[-1])
table(substring(metadata1$TCGA_id,14,15))

metadata1$sample <- ifelse(substring(metadata1$TCGA_id,14,15)=="11","control","cancer")
metadata1$sample <- factor(metadata1$sample,levels = c('control','cancer'))

expr2 = avereps(mRNAexp[,-1],ID = mRNAexp$gene_name) %>% as.data.frame()
expr3 <- expr2[rowMeans(expr2)>1,] 
expr4<-apply(expr3,MARGIN = 1,round) %>%  t() %>% as.data.frame()
dds <-DESeqDataSetFromMatrix(countData=expr4, colData=metadata1,design=~sample)

dds$sample<- relevel(dds$sample, ref = "control")
dds <- DESeq(dds)
LUAD_allDEG2 <- as.data.frame(results(dds))
#save(LUAD_allDEG2,file = "LUAD_allDEG2.Rdata")

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, "sample")
LUAD_mRNA_exprSet_vst <- as.data.frame(assay(vsd))
#save(LUAD_mRNA_exprSet_vst,file='LUAD_mRNA_exprSet_vst.Rdata')


padj = 0.05
foldChange= 2
diff_signif_DESeq2 = LUAD_allDEG2[(LUAD_allDEG2$padj < padj&(LUAD_allDEG2$log2FoldChange>foldChange 
                                                             |LUAD_allDEG2$log2FoldChange<(-foldChange))),]
mRNA_expr<-LUAD_mRNA_exprSet_vst[rownames(diff_signif_DESeq2),] %>% t() %>% as.data.frame()


