{
library(Biobase)
library(BisqueRNA)
library(dplyr)
library(tidyverse)
library(readxl)
library(msigdbr)
library(GSEABase)
library("AnnotationDbi")
library("org.Mm.eg.db")
library(sp)
library(Seurat)
library(reshape)
library(ggpubr)
}
####
msig_cell_markers = msigdbr(species = "Mus musculus", category = "C8")

######bulk data with FPMK
bulk.matrix <- raw_fpkm %>% 
  dplyr::select(1:27) %>% 
  column_to_rownames(var="gene_id") %>% 
  filter(rowSums(.) >= 10) %>%
  as.matrix()

dim(bulk.matrix)
# bulk data with counts 
bulk.matrix <- raw_count %>% 
  filter(gene_biotype == "protein_coding") %>% 
  distinct(gene_name, .keep_all = TRUE) %>% 
  rownames_to_column() %>% 
  column_to_rownames(var = "gene_name") %>% 
  dplyr::select(2:27) %>% 
  as.matrix()


bulk <- Biobase::ExpressionSet(assayData = bulk.matrix)

####marker based 
{
cell_marker <- read_excel("cell markers_mice heart.cell_markers.db.xlsx")

cell_marker$ENSEMBL <- as.character(mapIds(org.Mm.eg.db, # convert IDs
                                          keys = cell_marker$`Cell marker`, 
                                          keytype="SYMBOL", 
                                          column = "ENSEMBL")) 

markers <- cell_marker %>%
  dplyr::select(`Cell marker`,`Cell name`) %>% 
  distinct(`Cell marker`, .keep_all = TRUE)

colnames(markers) <- c("gene", "cluster")



##
deconvo_res <- BisqueRNA::MarkerBasedDecomposition(bulk.eset= bulk, 
                                                   markers = markers, 
                                                   weighted=F,
                                                   ct_col = "cluster", 
                                                   gene_col = "gene", 
                                                   unique_markers = FALSE, 
                                                   verbose = TRUE)

length(intersect(rownames(bulk), markers$gene))

sampleNames(bulk.eset)
featureNames(bulk.eset)

length(intersect(unique(markers$gene), rownames(bulk)))
}
################ R object saurat: 
load("facs_Heart_seurat_tiss.Robj")
rm(tiss)

### creating a seurat object with tables

sc_count <- read.csv("Heart-counts-scRNA.csv")


sc_count <- sc_count %>% 
  column_to_rownames("X")

sc_meta <- read.csv("annotations_facs-scRNA.csv")

heart_cells <- intersect(colnames(sc_count),sc_meta$cell)


sc_meta <- sc_meta %>% 
  filter(cell %in% heart_cells) %>% 
  filter(cell_ontology_class != "")

sc_count <- sc_count[,sc_meta$cell]

colnames(sc_count)[1]
sub("\\.", "-",sub("\\.", "_",colnames(sc_count)[1]))


colnames(sc_count) <- sub("\\.", "-",sub("\\.", "_",colnames(sc_count)))

sc_meta <- sc_meta %>% 
  mutate(cell.new = sub("\\.", "-",sub("\\.", "_",cell))) %>% 
  column_to_rownames("cell.new") %>% 
  dplyr::select(cell_ontology_class)

identical(sort(colnames(sc_count)), sort(row.names(sc_meta)))

sc.seu <- CreateSeuratObject(counts = sc_count, 
                             meta.data = sc_meta)

sc.seu@meta.data <- sc.seu@meta.data %>% 
  dplyr::mutate(cell.ident = cell_ontology_class) %>% 
  dplyr::select(cell.ident,nCount_RNA, nFeature_RNA)

colnames(sc.seu)
Cells(sc.seu)
Features(sc.seu)
Layers(sc.seu)
Assays(sc.seu)
Idents(sc.seu) <- "cell.ident"

sc.eset <- BisqueRNA::SeuratToExpressionSet(sc.seu, delimiter="-", position=2, version="v3")

sc.eset$SubjectName
sc.eset$cellType

table(sc.eset[["SubjectName"]])
table(sc.eset[["cellType"]])
sum(is.na(sc.eset[["SubjectName"]]))
sum(is.na(sc.eset[["cellType"]]))
####

res <- BisqueRNA::ReferenceBasedDecomposition(bulk, 
                                              sc.eset, 
                                              use.overlap=FALSE)
res$bulk.props

estimated_cell_prop <- as.data.frame(t(res$bulk.props))
estimated_cell_prop$group <- rep(c("WT Mock", "dbdb Mock", "WT Day 7", "dbdb Day 7"), c(7,5,7,7))

write.csv(estimated_cell_prop, "Deconvo scores.csv")

melted_est <- melt(estimated_cell_prop)
melted_est$group <- factor(melted_est$group, levels = c("WT Mock","WT Day 7", "dbdb Mock", "dbdb Day 7"))


### box plot
ggplot(melted_est, aes(x= variable, y=value, colour= group))+
  geom_boxplot(outlier.shape = NA, width = 0.8)+
  geom_point(position=position_jitterdodge(dodge.width = 0.8, jitter.width = 0.5),
             size = 2,
             alpha = 0.3)+
  theme_classic()+
  scale_color_manual(values=c("#3D3D3D","#e01327", "#8913ff", "#15a354"))+
  ylab("Deconv score") +
  xlab("")+
  ylim(0, 1.2) +
  theme(text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1))+
  ggtitle('')+
  geom_pwc(aes(group=group),
           method = "t_test", label ="p.adj.signif",
           p.adjust.method = 'none', tip.length = 0, 
           group.by = 'x.var', hide.ns = TRUE, label.size = 6)


















