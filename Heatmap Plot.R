{
library(dplyr)
library(ComplexHeatmap)
library(DESeq2)
require(circlize)
library(tibble)
}
### Heatmap for shared genes in the two comparisons:
{
gene_heatmap_shared_incomparisons <- FC_WT_db %>% 
  filter(Category == 'Sig in both (3665)') %>% 
  select(`gene_id (Wt)`) %>% 
  pull()

heatmap_data <- raw_fpkm %>% 
  filter(gene_id %in% gene_heatmap_shared_incomparisons) %>% 
  select(c(1:28)) 

heatmap_data <- column_to_rownames(heatmap_data, var = 'gene_name')

heatmap_data <- heatmap_data %>% 
  select(WT_M3, WT_M4,WT_M5,WT_M6,WT_M7,WT_M8,WT_M9, 
         WT_D7_1,WT_D7_2,WT_D7_3,WT_D7_4,WT_D7_5,WT_D7_6,WT_D7_7,
         db_M1,db_M2,db_M3,db_M4,db_M5,
         db_D7_1,db_D7_2,db_D7_3,db_D7_4,db_D7_5,db_D7_6,db_D7_7)

meta_heat <- data.frame(
  row.names = colnames(heatmap_data), 
  groups = rep(c('WT Mock', 'WT Day 7','db/db Mock', 'db/db Day 7' ), c(7,7,5,7)),
  status = factor(rep(c('WT', 'db/db'), c(14,12)), levels= c('WT', 'db/db'))
)

heat <- t(scale(t(heatmap_data)))

myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

col <- list(status = c('WT' = "pink", 'db/db' = '#abe3f3'),
            groups = c('WT Mock' = '#3D3D3D','WT Day 7' = '#e01327', 'db/db Mock'= '#8913ff', 'db/db Day 7' = '#15a354'))

colan <- HeatmapAnnotation(
  df= meta_heat,
  which = 'col',
  col = col,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm')
)


both_sig <- FC_WT_db %>%
  filter(Category =="Sig in both (3665)") %>% 
  select(`gene_id (Wt)`, `gene_name(WT)`, `log2foldchange(WT)`, `log2foldchange(db)`, Category)

both_sig <- both_sig %>% 
  mutate(direction = ifelse(`log2foldchange(WT)` > 0 & `log2foldchange(db)` > 0,
                            "up in both", 
                            ifelse(`log2foldchange(WT)` < 0 & `log2foldchange(db)` > 0,
                                   "down in WT \n up in db/db", 
                                   ifelse(`log2foldchange(WT)` > 0 & `log2foldchange(db)` < 0,
                                          "up in WT \n down in db/db",
                                          ifelse(`log2foldchange(WT)`< 0 & `log2foldchange(db)` < 0,
                                                 "Down in both", NA)))))
clusters <- both_sig %>% 
  select(`gene_name(WT)`, direction) %>% 
  column_to_rownames(var = "gene_name(WT)") %>% 
  arrange(direction)


heat <- heat[rownames(clusters),]

hmap <- Heatmap(heat,
                split = clusters$direction,
                show_row_dend = FALSE,
                column_split = meta_heat$status,
                cluster_row_slices = FALSE,
                col = colorRamp2(myBreaks, myCol),
                top_annotation = colan,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                row_names_gp = gpar(fontsize = 7),
                heatmap_legend_param = list(
                  title='z-score'
                ))


hmap
}

### Heatmap for top changed genes in comparisons: 
{
top_genes <- db7vsdbmock_deg %>% 
  arrange(pvalue) %>% 
  dplyr::slice(1:50) %>% 
  select(gene_id) %>% 
  pull()

heatmap_data <- raw_fpkm %>% 
  filter(gene_id %in% top_genes) %>% 
  select(c(1:28)) 
  
heatmap_data <- column_to_rownames(heatmap_data, var = 'gene_name')

heatmap_data <- heatmap_data %>% 
  select(db_M1,db_M2,db_M3,db_M4,db_M5,
         db_D7_1,db_D7_2,db_D7_3,db_D7_4,db_D7_5,db_D7_6,db_D7_7)

meta_heat <- data.frame(
  row.names = colnames(heatmap_data), 
  groups = rep(c('db/db Mock','db/db Day 7'), c(5,7))
  )

heat <- t(scale(t(heatmap_data)))

myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

col <- list(groups = c('db/db Mock'= '#8913ff', 'db/db Day 7' = '#15a354'))

colan <- HeatmapAnnotation(
  df= meta_heat,
  which = 'col',
  col = col,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm')
)

hmap <- Heatmap(heat,
                show_row_dend = TRUE,
                cluster_row_slices = TRUE,
                col = colorRamp2(myBreaks, myCol),
                top_annotation = colan,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 12),
                heatmap_legend_param = list(
                  title='z-score'
                ))


hmap
}

### Heatmap for genes in pathways

my_genes <- GO_ORA_down@result %>% 
  filter(Description=="positive regulation of leukocyte activation") %>% 
  dplyr::select(geneID) %>% 
  pull()

my_genes <- strsplit(my_genes,split='/', fixed=TRUE)
my_genes<- my_genes[[1]]

heatmap_data <- raw_fpkm %>% 
  filter(gene_name %in% my_genes) %>% 
  dplyr::select(c(1:28)) 

heatmap_data <- column_to_rownames(heatmap_data, var = 'gene_name')

heatmap_data <- heatmap_data %>% 
  dplyr::select(WT_M3, WT_M4,WT_M5,WT_M6,WT_M7,WT_M8,WT_M9, 
                WT_D7_1,WT_D7_2,WT_D7_3,WT_D7_4,WT_D7_5,WT_D7_6,WT_D7_7,
                db_M1,db_M2,db_M3,db_M4,db_M5,
                db_D7_1,db_D7_2,db_D7_3,db_D7_4,db_D7_5,db_D7_6,db_D7_7)

meta_heat <- data.frame(
  row.names = colnames(heatmap_data), 
  groups = rep(c('WT Mock', 'WT Day 7','db/db Mock', 'db/db Day 7' ), c(7,7,5,7)),
  status = factor(rep(c('WT', 'db/db'), c(14,12)), levels= c('WT', 'db/db'))
)


heat <- t(scale(t(heatmap_data)))

myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

col <- list(status = c('WT' = "pink", 'db/db' = '#abe3f3'),
            groups = c('WT Mock' = '#3D3D3D','WT Day 7' = '#e01327', 'db/db Mock'= '#8913ff', 'db/db Day 7' = '#15a354'))


colan <- HeatmapAnnotation(
  df= meta_heat,
  which = 'col',
  col = col,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm')
)

hmap <- Heatmap(heat,
                show_row_dend = TRUE,
                cluster_row_slices = TRUE,
                col = colorRamp2(myBreaks, myCol),
                top_annotation = colan,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                column_split = meta_heat$status,
                row_names_gp = gpar(fontsize = 6),
                heatmap_legend_param = list(
                  title='z-score'
                ))


hmap
