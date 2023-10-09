library(tidyverse)
library(DESeq2)
library(tximport)

#library(GenomicFeatures)
#txdb <- makeTxDbFromGFF(file="reference/gencode.v42.primary_assembly.annotation.gtf")
#saveDb(x=txdb, file = "reference/gencode.v42.primary_assembly.annotation.TxDb")
#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME") 

#write_tsv(tx2gene,"reference/tx2gene.tsv")

tx2gene <- read_tsv("reference/tx2gene.tsv")

files <- c("results/KGN-A/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           "results/KGN-B/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           "results/KGN-C/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           "results/C1B6-A/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           "results/C1B6-B/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           "results/C1B6-C/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results")
           #"results/C1D10-A/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           #"results/C1D10-B/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           #"results/C1D10-C/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           #"results/B2D8-A/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           #"results/B2D8-B/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results",
           #"results/B2D8-C/rsem_quant/rep1/rep1PE_stranded_anno_rsem.isoforms.ENSG.results")

names(files) <- c("KGN.A","KGN.B","KGN.C","C1B6.A","C1B6.B","C1B6.C")

txi <- tximport(files, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)

q <- data.frame(txi$counts)
q$gene <- rownames(q)

sampleTable <- data.frame(condition = factor(c("WT","WT","WT","DKO","DKO","DKO")))

rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds$condition <- relevel(dds$condition,ref="WT")

keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]
dseq_obj <- DESeq(dds)

DKOvSKO <- results(dseq_obj, contrast=c("condition", "WT", "DKO")) #,

res <- lfcShrink(dseq_obj,
                 contrast=c("condition", "DKO", "WT"), res=DKOvSKO, type = 'normal')

t <- data.frame(res)
t$gene <- rownames(t)
t <- t %>% filter(!grepl('PAR_Y', gene))

t$gene <- gsub("\\..*","",rownames(t))

write_tsv(t,"DEseq2_lfc_DKOvsWT_030223.tsv", col_names = TRUE)

library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
ens2symbol <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id","gene_biotype"),
      filters=c('ensembl_gene_id'),
      values= t$gene,
      mart=ensembl)
colnames(ens2symbol) <- c("SYMBOL","ENSEMBL","GENE_BIOTYPE")

ens2symbol <- as_tibble(ens2symbol)
t <- inner_join(t, ens2symbol, by=c("gene"="ENSEMBL"))

#write_tsv(t,"C1B6vC1D10_diffexp_lfcSHrink_010823.tsv")
write_tsv(t,"C1B6vALLOTHER_diffexp_lfcSHrink_012723_Symbol.tsv", col_names = TRUE)

t <- read_tsv("C1B6vALLOTHER_diffexp_lfcSHrink_012723_Symbol.tsv", col_names = TRUE)

t.protein_coding <- t %>%
  filter(GENE_BIOTYPE == "protein_coding")

emt <- c("ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","TWIST2","CDH1","CDH2","FN1","VIM")
mAGCT <- c(unlist(mAGCTgenes$mGeneID))


mAGCTgenes <- read_tsv("Llano_mAGCTs_vs_WT_UP.txt", col_names = TRUE)

KEGG_TGFbgenes <- read_tsv("KEGG_TGF_BETA_SIGNALING_PATHWAY.v2022.1.Hs.grp", col_names = TRUE)
KEGG_TGFbgenes <- KEGG_TGFb[-1,]
colnames(KEGG_TGFbgenes) <- c("geneID")
KEGG_TGFb <- c(unlist(KEGG_TGFbgenes$geneID))

t.protein_coding <- t %>%
  filter(SYMBOL %in% emt)

library(EnhancedVolcano)
EnhancedVolcano(t.protein_coding,
                lab = t.protein_coding$SYMBOL,
                selectLab = emt,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 1e-1,)

ggsave(file="/Users/rthillman/Desktop/KGN_data_hub/encode_rna_seq/figures/EMT_genes_volcano_013123.png", width=8, height=8, dpi=300)


t <- data.frame(DKOvSKO)

norm_counts <- counts(dseq_obj, normalized = T)
fid <- "norm_counts_012723.gct" 
writeLines(c("#1.2", paste(nrow(norm_counts), ncol(norm_counts) - 2, collapse="\t")), fid, sep="\n")
write.table(norm_counts, file=fid, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", append = TRUE)

library(fgsea)

filt <- c("HERC3","TBCE")

res2 <- t.protein_coding %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>%
  dplyr::filter(!SYMBOL %in% filt)
#%>% 
  #distinct() #%>% 
  #group_by(SYMBOL) %>% 
  #summarize(stat=mean(stat))
    
 #   res2 %>% group_by(SYMBOL) %>% filter(n()>1)

ranks <- deframe(res2)
pathways.hallmark <- gmtPathways("msigdb/h.all.v2022.1.Hs.symbols.gmt")
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  dplyr::filter(padj < 0.1)

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

ggsave(file="/Users/rthillman/Desktop/KGN_data_hub/encode_rna_seq/figures/GSEA_proteincoding_013123.png", width=8, height=8, dpi=300)


library(pheatmap)
library(RColorBrewer)
vsd <- vst(dseq_obj, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("condition")) + theme_classic()



