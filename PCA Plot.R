library(readxl)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(plotly)

raw_fpkm <- read_xlsx("gene_fpkm.xls.xlsx")

fpkm <- raw_fpkm %>% 
  select(c(2:27))

filtered_fpkm <- fpkm %>% 
  filter(rowSums(.) >= 10)

meta <- rep(c('WT mock', 'dbdb mock', 'WT day 7', 'dbdb day 7' ), c(7,5,7,7))

meta_data <- as.data.frame(meta)
meta_data <- meta_data %>% 
  mutate(meta = factor(
    meta, levels=c('WT mock', 'dbdb mock', 'WT day 7', 'dbdb day 7' )
  ))

rownames(meta_data) <- colnames(filtered_fpkm)

filtered_fpkm <- round(filtered_fpkm)


dds <- DESeqDataSetFromMatrix(
  countData = filtered_fpkm,
  colData = meta_data,
  design = ~1
)



dds_norm <- vst(dds)

plotPCA(
  dds_norm,
  intgroup= 'meta',
  pcsToUse = 1:2
) + scale_color_manual(values = c('blue', 'red', 'green', 'orange')) + theme_bw()
# Extract PCA coordinates
pca_data <- plotPCA(dds_norm, intgroup = 'meta', returnData = TRUE)
library(ggforce)
?plotPCA
# Plot PCA with confidence ellipses
ggplot(pca_data, aes(x = PC1, y = PC2, color = meta)) +
  geom_point(size= 3, alpha=0.75) +
  geom_mark_ellipse(aes(color = meta))+
  #stat_ellipse(aes(group = meta), type = "t", level= 0.99) +
  scale_color_manual(values = c('#3D3D3D', '#8913ff', '#e01327', '#15a354')) +
  theme_bw() +
  ylim(-13,18) +
  xlim(-18,18) +
  labs(x = paste("PC1: 39%"),
       y = paste("PC2: 21%")) +
  theme(text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5), 
        legend.position="bottom", 
        legend.title = element_blank())

### calculating the variance for each component
dds <- rlog( dds, blind = T )
rv <- rowVars(assay(dds_norm))
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(dds_norm)[select,]))

PC1=pca$x[,1]
PC2=pca$x[,2]
PC3=pca$x[,3]

percentVar <- pca$sdev^2/sum(pca$sdev^2)

components <- data.frame(meta,
                         PC1, 
                         PC2, 
                         PC3)
colnames(components) <- c('meta', 'PC1: 39%', 'PC2: 21%', 'PC3: 8%')

##### 3D PCA
fig <- plot_ly(components, x = ~`PC1: 39%`, y = ~`PC2: 21%`, z = ~`PC3: 8%`, 
               color = ~components$meta, 
               colors = c('#15a354', '#8913ff', '#e01327', '#3D3D3D'),
               marker = list(size = 6))


fig <- fig %>%
  layout(
    scene = list(bgcolor = "white")
  )


library(rgl)
install.packages("Morpho")
library(Morpho)

ellipse <- ellipse3d(cov(components[1:7,2:4]),centre=apply(components[1:7,2:4], 2, mean)) #WT mock
ellipse1<- ellipse3d(cov(components[8:12,2:4]),centre=apply(components[8:12,2:4], 2, mean)) #dbdb mock 
ellipse2 <- ellipse3d(cov(components[13:19,2:4]),centre=apply(components[13:19,2:4], 2, mean)) #WT day 7
ellipse3 <- ellipse3d(cov(components[20:26,2:4]),centre=apply(components[20:26,2:4], 2, mean)) #dbdb day 7

ellipsoid <- Morpho::quad2trimesh(ellipse)
ellipsoid1 <- Morpho::quad2trimesh(ellipse1)
ellipsoid2 <- Morpho::quad2trimesh(ellipse2)
ellipsoid3 <- Morpho::quad2trimesh(ellipse3)

# Define the color palette
color_palette <- colorRamp(c('#15a354', '#8913ff', '#e01327', '#3D3D3D'))

# Match each ellipse to the appropriate color in the palette
fig <- plot_ly(colors = color_palette) %>% 
  # Add scatter3d trace
  add_trace(data = components, 
            x = ~`PC1: 39%`, y = ~`PC2: 21%`, z = ~`PC3: 8%`,
            type = "scatter3d", mode = 'markers', 
            marker = list(size = 3),
            color = ~components$meta) %>% 
  # Add mesh3d traces for ellipses, matching scatter colors
  add_trace(x = ellipse$vb[1,], y = ellipse$vb[2,], z = ellipse$vb[3,], 
            type = 'mesh3d', alphahull = 0, opacity = 0.2,
            facecolor = rep("#3D3D3D", ncol(ellipsoid$it))) %>%  
  add_trace(x = ellipse1$vb[1,], y = ellipse1$vb[2,], z = ellipse1$vb[3,], 
            type = 'mesh3d', alphahull = 0, opacity = 0.2,
            facecolor = rep("#8913ff", ncol(ellipsoid$it))
            ) %>%  
  add_trace(x = ellipse2$vb[1,], y = ellipse2$vb[2,], z = ellipse2$vb[3,], 
            type = 'mesh3d', alphahull = 0, opacity = 0.2,
            facecolor = rep("#e01327", ncol(ellipsoid$it))
            ) %>% 
  add_trace(x = ellipse3$vb[1,], y = ellipse3$vb[2,], z = ellipse3$vb[3,], 
            type = 'mesh3d', alphahull = 0, opacity = 0.2,
            facecolor = rep("#15a354", ncol(ellipsoid$it))
            ) %>%  
  # Set background color
  layout(
    scene = list(bgcolor = "white"))

fig

