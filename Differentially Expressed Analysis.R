library(DESeq2)
library(tidyverse)
### raed the count table 
raw_count <- read.csv("gene_count.csv", header = T, row.names = "gene_id")
#### extract the desired samples
WT_raw <- cbind(raw_count[,1:7], raw_count[,13:19])
#### make the relevant meta data 
WT_meta <- as.data.frame(rep(c("WT mock", "WT day 7"),c(7,7)))
colnames(WT_meta) <- "condition"
row.names(WT_meta) <- colnames(WT_raw)
WT_meta$condition <- as.factor(WT_meta$condition)
### generate the dds object
dds <- DESeqDataSetFromMatrix(countData = WT_raw,
                              colData = WT_meta,
                              design = ~ condition)
dds <- DESeq(dds, quiet = T)

## function for DEG 
generate_DE_results <- function (dds, comparisons, padjcutoff = 1, log2cutoff = 0, cpmcutoff = 0) {
  # generate average counts per million metric from raw count data 
  raw_counts <- counts(dds, normalized = F)
  cpms <- enframe(rowMeans(edgeR::cpm(raw_counts)))
  colnames(cpms) <- c("gene_id", "avg_cpm")
  
  # extract DESeq results between the comparisons indicated
  res <- results(dds, contrast = c("condition", comparisons[1], comparisons[2]))[,-c(3,4)]
  
  # annotate the data with gene name and average counts per million value
  res <- as_tibble(res, rownames = "gene_id")
  # read in the annotation and append it to the data
  my_annotation <-  read.csv("gene_count.csv", header = T, row.names = NULL)
  my_annotation <- my_annotation [,c(1,28,34)] ## numbers of cols should be for: gene_id,gene_name,gene_biotype
  res <- left_join(res, my_annotation, by = c("gene_id" = "gene_id"))
  # append the average cpm value to the results data
  res <- left_join(res, cpms, by = c("gene_id" = "gene_id"))
  
  # combine normalized counts with entire DE list
  normalized_counts <- round(counts(dds, normalized = TRUE),3)
  pattern <- str_c(comparisons[1], "|", comparisons[2])
  combined_data <- as_tibble(cbind(res, normalized_counts[,grep(pattern, colnames(normalized_counts))] ))
  combined_data <- combined_data[order(combined_data$log2FoldChange, decreasing = T),]
  
  # make ordered rank file for GSEA, selecting only protein coding genes
  res_prot <- res[which(res$gene_biotype == "protein_coding"),]
  res_prot_ranked <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("gene_name", "log2FoldChange")]
  res_prot_ranked <- na.omit(res_prot_ranked)
  res_prot_ranked$gene_name <- str_to_upper(res_prot_ranked$gene_name)
  
  # generate sorted lists with the indicated cutoff values
  res <- res[order(res$padj, decreasing=FALSE ),]
  de_genes_padj <- res[which(res$padj < padjcutoff),]
  de_genes_log2f <- res[which(abs(res$log2FoldChange) > log2cutoff & res$padj < padjcutoff),]
  de_genes_cpm <- res[which(res$avg_cpm > cpmcutoff & res$padj < padjcutoff),]
  
  # write output to files
  write.csv (de_genes_padj, file = paste0(comparisons[1], "_vs_", comparisons[2], "_padj_cutoff.csv"), row.names =F)
  write.csv (de_genes_log2f, file = paste0(comparisons[1], "_vs_", comparisons[2], "_log2f_cutoff.csv"), row.names =F)
  write.csv (de_genes_cpm, file = paste0(comparisons[1], "_vs_", comparisons[2], "_cpm_cutoff.csv"), row.names =F)
  write.csv (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_allgenes.csv"), row.names =F)
  write.table (res_prot_ranked, file = paste0(comparisons[1], "_vs_", comparisons[2], "_rank.rnk"), sep = "\t", row.names = F, quote = F)
  
  writeLines( paste0("For the comparison: ", comparisons[1], "_vs_", comparisons[2], ", out of ", nrow(combined_data), " genes, there were: \n", 
                     nrow(de_genes_padj), " genes below padj ", padjcutoff, "\n",
                     nrow(de_genes_log2f), " genes below padj ", padjcutoff, " and above a log2FoldChange of ", log2cutoff, "\n",
                     nrow(de_genes_cpm), " genes below padj ", padjcutoff, " and above an avg cpm of ", cpmcutoff, "\n",
                     "Gene lists ordered by log2fchange with the cutoffs above have been generated.") )
  gene_count <- tibble (cutoff_parameter = c("padj", "log2fc", "avg_cpm" ), 
                        cutoff_value = c(padjcutoff, log2cutoff, cpmcutoff), 
                        signif_genes = c(nrow(de_genes_padj), nrow(de_genes_log2f), nrow(de_genes_cpm)))
  invisible(gene_count)
}

generate_DE_results(dds,comparisons = c("WT day 7", "WT mock"),padjcutoff = 0.05, log2cutoff = 1, cpmcutoff = 0)

