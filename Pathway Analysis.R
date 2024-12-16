{
library(msigdbr)
library(dplyr)
library(clusterProfiler)
library("AnnotationDbi")
library("org.Mm.eg.db")
library(ggplot2)
library(orthogene)
library(enrichplot)
library("pathview")
}

## retrieve the HALLMARK and BioPlanet databases:
hall_gene_sets = msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
hall_gene_sets$gs_name <- gsub("HALLMARK_", "", hall_gene_sets[, 1])


Bioplanet_database <- read.csv("Bioplanet database.csv")
Bioplanet_database <- Bioplanet_database %>% 
  dplyr::select(PATHWAY_NAME, GENE_SYMBOL)
################### Prepare genelist data frames:
db7vsdbmock_deg

genelist_up <- db7vsdbmock_deg %>%
  filter(pvalue < 0.05 & log2FoldChange > 0) %>% 
  dplyr::select(gene_id, gene_name, log2FoldChange, gene_ENTREZ)


genelist_down <- db7vsdbmock_deg %>%
  filter(pvalue < 0.05 & log2FoldChange < 0) %>% 
  dplyr::select(gene_id, gene_name,log2FoldChange, gene_ENTREZ)
#genelist for FCvsFC genes: 
 FC_WT_db
genelist_up <- FC_WT_db %>% 
   filter(Category == "Sig in both (3665)" &
            `log2foldchange(WT)` > 0 &
            `log2foldchange(db)` > 0) %>% 
  dplyr::select(`gene_id (Wt)`, `gene_name(WT)`, gene_ENTREZ)

genelist_down <- FC_WT_db %>% 
  filter(Category == "Sig in both (3665)" &
           `log2foldchange(WT)` < 0 &
           `log2foldchange(db)` < 0) %>% 
  dplyr::select(`gene_id (Wt)`, `gene_name(WT)`, gene_ENTREZ)

############## ORA analysis
## GO database
{
### Up: 
GO_ORA_up <- enrichGO(gene = genelist_up$gene_id,
                           universe = db7vsWT7_deg$gene_id,
                           keyType = "ENSEMBL",
                           OrgDb = org.Mm.eg.db,
                           ont = "ALL",
                           pAdjustMethod = "fdr",
                           pvalueCutoff = 1,
                           qvalueCutoff = 1,
                           readable = TRUE)

dotplot(GO_ORA_up, 
        showCategory=10,
        split = "ONTOLOGY",
        orderBy = "ONTOLOGY",
        label_format = 100)+
  theme(text = element_text(size = 16, face = "bold"))

# cnetplot
## create a named vector:
genelist_up.vec <- genelist_up$log2FoldChange
names(genelist_up.vec) <- genelist_up$gene_name

cnetplot(GO_ORA_up,
         showCategory = 6,
         color.params = list(foldChange = genelist_up.vec))
#emapplot
GO_ORA_up <- pairwise_termsim(GO_ORA_up)
emapplot(GO_ORA_up,
         showCategory= 25)

write.csv(GO_ORA_up@result, 
          file= "Pathway Analysis/FCvsFC_GO_up (up in both).csv")

### DoWn:
GO_ORA_down <- enrichGO(gene = genelist_down$gene_id,
                      universe = db7vsWT7_deg$gene_id,
                      keyType = "ENSEMBL",
                      OrgDb = org.Mm.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "fdr",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      readable = TRUE)

dotplot(GO_ORA_down, 
        showCategory=10,
        split = "ONTOLOGY",
        orderBy = "ONTOLOGY",
        label_format = 100)+
  theme(text = element_text(size = 16, face = "bold"))

write.csv(GO_ORA_down@result, 
          file= "Pathway Analysis/FCvsFC_GO_down (down in both).csv")
}
## KEGG database
{
  
  WT7vsWTmock_deg$gene_ENTREZ <- as.character(mapIds(org.Mm.eg.db, # convert IDs
                                                        keys = WT7vsWTmock_deg$gene_id, 
                                                        keytype="ENSEMBL", 
                                                        column = "ENTREZID"))  

### Up:
  
  genelist_up$gene_ENTREZ <- as.character(mapIds(org.Mm.eg.db, # convert IDs
                                                 keys = genelist_up$gene_id, 
                                                 keytype="ENSEMBL", 
                                                 column = "ENTREZID"))

KEGG_ORA_up <- enrichKEGG(gene = genelist_up$gene_ENTREZ,
                               universe = WT7vsWTmock_deg$gene_ENTREZ,
                               keyType = "kegg",
                               organism  = "mmu",
                               pAdjustMethod = "fdr",
                               pvalueCutoff = 1,
                               qvalueCutoff = 1
                               )

KEGG_ORA_up <- KEGG_ORA_up %>% 
  mutate(Description = sub(" - [^-]*$", "",KEGG_ORA_up@result$Description))

dotplot(KEGG_ORA_up, 
        showCategory=20,
        label_format = 100)+
  theme(text = element_text(size = 16, face = "bold"))

##creat vector for foldchange:
genelist_up.vec <- genelist_up$log2FoldChange
names(genelist_up.vec) <- genelist_up$gene_ENTREZ

pathview(gene.data  = genelist_up.vec,
         pathway.id = "mmu05171",
         species    = "mmu",
         limit      = list(gene=max(abs(genelist_up.vec)), cpd=1))


write.csv(KEGG_ORA_up@result, 
          file= "Pathway Analysis/FCvsFC_KEGG_up (up in both).csv")

###Down: 

genelist_down$gene_ENTREZ <- as.character(mapIds(org.Mm.eg.db, # convert IDs
                                                 keys = genelist_down$`gene_id (Wt)`, 
                                                 keytype="ENSEMBL", 
                                                 column = "ENTREZID"))

KEGG_ORA_down <- enrichKEGG(gene = genelist_down$gene_ENTREZ,
                          universe = FC_WT_db$gene_ENTREZ,
                          keyType = "kegg",
                          organism  = "mmu",
                          pAdjustMethod = "fdr",
                          pvalueCutoff = 1,
                          qvalueCutoff = 1
)

KEGG_ORA_down <- KEGG_ORA_down %>% 
  mutate(Description = sub(" - [^-]*$", "",KEGG_ORA_down@result$Description))

dotplot(KEGG_ORA_down, 
        showCategory=20,
        label_format = 100)+
  theme(text = element_text(size = 16, face = "bold"))

write.csv(KEGG_ORA_down@result, 
          file= "Pathway Analysis/FCvsFC_KEGG_down (down in both).csv")
}
## Hallmark database
{
  ###Up: 
  HALL_ORA_up <- enricher(gene= genelist_up$gene_ENTREZ, 
                               universe = FC_WT_db$gene_ENTREZ,
                               TERM2GENE = hall_gene_sets,
                               pvalueCutoff = 1,
                               pAdjustMethod = 'fdr')
  HALL_ORA_up@result$Description <- gsub('HALLMARK_', '', HALL_ORA_up@result$Description)
  
  dotplot(HALL_ORA_up, 
          showCategory=20,
          label_format = 100)+
    theme(text = element_text(size = 16, face = "bold"))
  
  write.csv(HALL_ORA_up@result, 
            file= "Pathway Analysis/HALLMARK correct labels/FC_WT_db_HALLMARK_up (up in both).csv")
  
  ###Down: 
  HALL_ORA_down <- enricher(gene= genelist_down$gene_ENTREZ, 
                          universe = FC_WT_db$gene_ENTREZ,
                          TERM2GENE = hall_gene_sets,
                          pvalueCutoff = 1,
                          pAdjustMethod = 'fdr')
  
  HALL_ORA_down@result$Description <- gsub('HALLMARK_', '', HALL_ORA_down@result$Description)
  
  dotplot(HALL_ORA_down, 
          showCategory=20,
          label_format = 100)+
    theme(text = element_text(size = 16, face = "bold"))
  
  write.csv(HALL_ORA_down@result, 
            file= "Pathway Analysis/HALLMARK correct labels/FC_WT_db_HALLMARK_down (down in both).csv")
}
## Bioplanet database
{
  
  ortho_for_deg <- map_genes(genes = FC_WT_db$`gene_name(WT)`,
                             species = 'human',
                             drop_na = FALSE)
  
  ortho_for_deg <- ortho_for_deg %>% 
    distinct(input_number, .keep_all = TRUE) %>% 
    dplyr::select(input, name)
  
  FC_WT_db$gene_ortho <- ortho_for_deg$name
##Up:
##preparing datasets:
orthologe_for_Bio <- map_genes(genes = genelist_up$`gene_name(WT)`,
                               species = "human",
                               drop_na = FALSE)

orthologe_for_Bio <- orthologe_for_Bio %>% 
  distinct(input_number, .keep_all = TRUE) %>% 
  dplyr::select(input, name)

genelist_up$gene_ortho <- orthologe_for_Bio$name

## Run ORA analysis 

BP_ORA_up <- enricher(gene = genelist_up$gene_ortho,
                   universe = FC_WT_db$gene_ortho,
                   TERM2GENE = Bioplanet_database,
                   pvalueCutoff = 1,
                   pAdjustMethod = 'fdr')

dotplot(BP_ORA_up, 
        showCategory=20,
        label_format = 80)+
  theme(text = element_text(size = 16, face = "bold"))

write.csv(BP_ORA_up@result, 
          file= "Pathway Analysis/FCvsFC_Biopla_up (up in both).csv")

## Down: 
##preparing datasets:
orthologe_for_Bio <- map_genes(genes = genelist_down$`gene_name(WT)`,
                               species = "human",
                               drop_na = FALSE)

orthologe_for_Bio <- orthologe_for_Bio %>% 
  distinct(input_number, .keep_all = TRUE) %>% 
  dplyr::select(input, name)

genelist_down$gene_ortho <- orthologe_for_Bio$name

## Run ORA analysis 

BP_ORA_down <- enricher(gene = genelist_down$gene_ortho,
                      universe = FC_WT_db$gene_ortho,
                      TERM2GENE = Bioplanet_database,
                      pvalueCutoff = 1,
                      pAdjustMethod = 'fdr')

dotplot(BP_ORA_down, 
        showCategory=20,
        label_format = 80)+
  theme(text = element_text(size = 16, face = "bold"))

write.csv(BP_ORA_down@result, 
          file= "Pathway Analysis/FCvsFC_Biopla_down (down in both).csv")
}

############### GSEA analysis 
{
# making the required gene list 
## with FC in DEG:
genelist_WT7vsWTm <- WT7_WTmock %>% 
  select(log2FoldChange) %>% 
  pull()

names(genelist_WT7vsWTm) <- as.character(mapIds(org.Mm.eg.db, # convert IDs
                                                keys = WT7_WTmock$gene_id, 
                                                keytype="ENSEMBL", 
                                                column = "ENTREZID"))

genelist_WT7vsWTm <- sort(genelist_WT7vsWTm, decreasing = TRUE)

head(genelist_WT7vsWTm)
###### with the gene list used in NOVOgene
new_genelist_WT7vsWTm <- readxl::read_xlsx("ranked_gene_list_WT_Day7_versus_WT_Mock.tsv.xlsx")

new_genelist_WT7vsWTm$ENTREZ <- as.character(mapIds(org.Mm.eg.db,
                                                    keys = new_genelist_WT7vsWTm$TITLE, 
                                                    keytype="ENSEMBL", 
                                                    column = "ENTREZID"))
sum(is.na(new_genelist_WT7vsWTm$ENTREZ))

new_genelist_WT7vsWTm_vec <- as.vector(new_genelist_WT7vsWTm$SCORE) 

names(new_genelist_WT7vsWTm_vec) <- new_genelist_WT7vsWTm$TITLE

new_genelist_WT7vsWTm_vec <- sort(new_genelist_WT7vsWTm_vec, decreasing = TRUE)


head(new_genelist_WT7vsWTm_vec)
sum(is.na(names(new_genelist_WT7vsWTm_vec)))

### run the GSEA with gseGO function:
GO_GSEA_WT7_WTm <- gseGO(geneList = new_genelist_WT7vsWTm_vec,
                         OrgDb = org.Mm.eg.db,
                         keyType = "ENSEMBL",
                         minGSSize = 15,
                         maxGSSize = 5000,
                         pvalueCutoff = 0.01,
                         ont = "ALL",
                         pAdjustMethod = "fdr"
                         )


GO_GSEA_WT7_WTm@result %>% 
  filter(NES > 0 & p.adjust < 0.25) %>% 
  dplyr::select(ID) %>% 
  count()

}


