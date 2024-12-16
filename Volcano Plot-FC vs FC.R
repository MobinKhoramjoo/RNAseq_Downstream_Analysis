{
library(readxl)
library(dplyr)
library(ggplot2)
library(grid)
library(ggrepel)
}

db7vsdbmock_deg <- read_excel("dbdb_Day7vsdbdb_Mock_deg.xls.xlsx")
db7vsWT7_deg <- read_excel("dbdb_Day7vsWT_Day7_deg.xls.xlsx")
dbmockvsWTmock_deg <- read_excel("dbdb_MockvsWT_Mock_deg.xls.xlsx")
WT7vsWTmock_deg <- read_excel("WT_Day7vsWT_Mock_deg.xls.xlsx")

############# dbdb day 7 vs dbdb mock
{
  db7_mock <- db7vsdbmock_deg %>% 
    select(gene_id, log2FoldChange, pvalue, padj, gene_name)
  
  db7_mock$pvalue <- replace(db7_mock$pvalue, is.na(db7_mock$pvalue), 1)
  
  cat_db7_mock <- c() #create a vector for categories
  for (i in 1:length(row.names(db7_mock))) { #assign corresponding names to vector
    if(db7_mock$log2FoldChange[i] > 0 & db7_mock$pvalue[i] < 0.05){
      cat_db7_mock[i] <- "Upregulated"
    }
    else if(db7_mock$log2FoldChange[i] < 0 & db7_mock$pvalue[i] < 0.05){
      cat_db7_mock[i] <- "Downregulated"
    }
    else{
      cat_db7_mock[i] <- "Non-significant"
    }
  }
  
  db7_mock <- db7_mock %>% # assigning the vector to a column 
    mutate(category = cat_db7_mock)
  
  showing_label_db7_mock <- c() #create a vector for the molecules that should be shown 
  for (i in 1:length(row.names(db7_mock))) { # assigning corresponding lables
    if(db7_mock$category[i] == "Upregulated" | db7_mock$category[i] == "Downregulated"){
      showing_label_db7_mock[i] <- db7_mock$gene_name[i]
    }
    else {
      showing_label_db7_mock[i] <- NA
    }
  }
  
  db7_mock <- db7_mock %>% #assigning the created vector to a column
    mutate(label = showing_label_db7_mock)
  
  n_up <- db7_mock %>% 
    filter(category == 'Upregulated') %>% 
    count() %>% 
    pull()
  
  n_down <- db7_mock %>% 
    filter(category == 'Downregulated') %>% 
    count() %>% 
    pull()
  
  ggplot(data = db7_mock, 
         aes(x = log2FoldChange, 
             y = -log10(pvalue), 
             colour=category)) +
    geom_point(alpha=0.4, size=2) +
    geom_text_repel(aes(label= label),size = 3, force = 1.5, max.overlaps = 15) +
    scale_color_manual(values=c("blue", "grey","red"),
                       labels = c(paste("Downregulated (",n_down,")",sep = ""), "Non-significant", paste("Upregulated (",n_up,")",sep = "")))+
    geom_vline(xintercept=c(-1,1) ,lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
    labs(x=expression(log[2]~(FC)),
         y=expression(-log[10]~(adj-p)),
         title= expression(italic(dbdb)~ "day 7 vs" ~italic(dbdb)~"mock"),
         caption = expression("adj-p < 0.05    " ~ log[2]~"(FC) > 1")) +
    theme_bw()+
    theme(text = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5), 
          legend.position="bottom", 
          legend.title = element_blank())
  
}

############ dbdb day 7 vs WT day 7
{
  db7_WT7 <- db7vsWT7_deg %>% 
    select(gene_id, log2FoldChange, pvalue, padj)
  
  db7_WT7$pvalue[is.na(db7_WT7$pvalue)] <- 1
  
  cat_db7_WT7 <- c() #create a vector for categories
  for (i in 1:length(row.names(db7_WT7))) { #assign corresponding names to vector
    if(db7_WT7$log2FoldChange[i] > 0 & db7_WT7$pvalue[i] < 0.05){
      cat_db7_WT7[i] <- "Upregulated"
    }
    else if(db7_WT7$log2FoldChange[i] < 0 & db7_WT7$pvalue[i] < 0.05){
      cat_db7_WT7[i] <- "Downregulated"
    }
    else{
      cat_db7_WT7[i] <- "Non-significant"
    }
  }
  
  db7_WT7 <- db7_WT7 %>% # assigning the vector to a column 
    mutate(category = cat_db7_WT7)
  
  db7_WT7 %>% 
    filter(category == 'Upregulated') %>% 
    count()
  
  ggplot(data = db7_WT7, 
         aes(x = log2FoldChange, 
             y = -log10(pvalue), 
             colour=category)) +
    geom_point(alpha=0.4, size=3) +
    scale_color_manual(values=c("blue", "grey","red"),
                       labels = c("Downregulated (3667)", "Non-significant", "Upregulated (3918)"))+
    geom_vline(xintercept=0 ,lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)",
         title="dbdb day 7 vs WT day 7")  +
    theme_bw()+
    theme(text = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5), 
          legend.position="bottom", 
          legend.title = element_blank())
  
}

########### dbdb mock vs WT mock 
{
  dbmock_WTmock <- dbmockvsWTmock_deg %>% 
    select(gene_id, log2FoldChange, pvalue, padj)
  
  dbmock_WTmock$pvalue[is.na(dbmock_WTmock$pvalue)] <- 1
  
  cat_dbmock_WTmock <- c() #create a vector for categories
  for (i in 1:length(row.names(dbmock_WTmock))) { #assign corresponding names to vector
    if(dbmock_WTmock$log2FoldChange[i] > 0 & dbmock_WTmock$pvalue[i] < 0.05){
      cat_dbmock_WTmock[i] <- "Upregulated"
    }
    else if(dbmock_WTmock$log2FoldChange[i] < 0 & dbmock_WTmock$pvalue[i] < 0.05){
      cat_dbmock_WTmock[i] <- "Downregulated"
    }
    else{
      cat_dbmock_WTmock[i] <- "Non-significant"
    }
  }
  
  dbmock_WTmock <- dbmock_WTmock %>% # assigning the vector to a column 
    mutate(category = cat_dbmock_WTmock)
  
  dbmock_WTmock %>% 
    filter(category == 'Downregulated') %>% 
    count()
  
  ggplot(data = dbmock_WTmock, 
         aes(x = log2FoldChange, 
             y = -log10(pvalue), 
             colour=category)) +
    geom_point(alpha=0.4, size=3) +
    scale_color_manual(values=c("blue", "grey","red"),
                       labels = c("Downregulated (3634)", "Non-significant", "Upregulated (1771)"))+
    geom_vline(xintercept=0 ,lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)",
         title="dbdb mock vs WT mock")  +
    theme_bw()+
    theme(text = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5), 
          legend.position="bottom", 
          legend.title = element_blank())
  
}

########## WT day 7 vs WT mock
{
  WT7_WTmock <- WT7vsWTmock_deg %>% 
    select(gene_id, log2FoldChange, pvalue, padj, gene_name)
  
  WT7_WTmock$pvalue <- replace(WT7_WTmock$pvalue, is.na(WT7_WTmock$pvalue), 1)
  
  cat_WT7_WTmock <- c() #create a vector for categories
  for (i in 1:length(row.names(WT7_WTmock))) { #assign corresponding names to vector
    if(WT7_WTmock$log2FoldChange[i] > 0 & WT7_WTmock$pvalue[i] < 0.05){
      cat_WT7_WTmock[i] <- "Upregulated"
    }
    else if(WT7_WTmock$log2FoldChange[i] < 0 & WT7_WTmock$pvalue[i] < 0.05){
      cat_WT7_WTmock[i] <- "Downregulated"
    }
    else{
      cat_WT7_WTmock[i] <- "Non-significant"
    }
  }
  
  WT7_WTmock <- WT7_WTmock %>% # assigning the vector to a column 
    mutate(category = cat_WT7_WTmock)
  
  WT7_WTmock %>% 
    filter(category == 'Upregulated') %>% 
    count()
  
  ggplot(data = WT7_WTmock, 
         aes(x = log2FoldChange, 
             y = -log10(padj), 
             colour=category)) +
    geom_point(alpha=0.4, size=3) +
    scale_color_manual(values=c("blue", "grey","red"),
                       labels = c("Downregulated (1160)", "Non-significant", "Upregulated (1211)"))+
    geom_vline(xintercept=0 ,lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",
         y="-log10 (p-value)",
         title="WT day 7 vs WT mock")  +
    theme_bw()+
    xlim(-7.5, 7.5) +
    theme(text = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5), 
          legend.position="bottom", 
          legend.title = element_blank())
  
}

##volcano plot function:
my_volcano_plot <- function(dataframe_to_use, plot_title, repel_force = 1.5, repel_max_overlaps = 15) {
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  
  New_data <- dataframe_to_use %>% 
    select(gene_id, log2FoldChange, pvalue, padj, gene_name)
  
  New_data$pvalue <- replace(New_data$pvalue, is.na(New_data$pvalue), 1)
  
  cat <- c() # create a vector for categories
  for (i in 1:length(row.names(New_data))) { # assign corresponding names to vector
    if (New_data$log2FoldChange[i] > 0 & New_data$pvalue[i] < 0.05) {
      cat[i] <- "Upregulated"
    } else if (New_data$log2FoldChange[i] < 0 & New_data$pvalue[i] < 0.05) {
      cat[i] <- "Downregulated"
    } else {
      cat[i] <- "Non-significant"
    }
  }
  
  New_data <- New_data %>% # assigning the vector to a column 
    mutate(category = cat)
  
  showing_label <- c() # create a vector for the molecules that should be shown 
  for (i in 1:length(row.names(New_data))) { # assigning corresponding labels
    if (New_data$category[i] == "Upregulated" | New_data$category[i] == "Downregulated") {
      showing_label[i] <- New_data$gene_name[i]
    } else {
      showing_label[i] <- NA
    }
  }
  
  New_data <- New_data %>% # assigning the created vector to a column
    mutate(label = showing_label)
  
  n_up <- New_data %>% 
    filter(category == 'Upregulated') %>% 
    count() %>% 
    pull()
  
  n_down <- New_data %>% 
    filter(category == 'Downregulated') %>% 
    count() %>% 
    pull()
  
  ggplot(data = New_data, 
         aes(x = log2FoldChange, 
             y = -log10(pvalue), 
             colour = category)) +
    geom_point(alpha = 0.4, size = 2) +
    #geom_text_repel(aes(label = label), size = 3, force = repel_force, max.overlaps = repel_max_overlaps) +
    scale_color_manual(values = c("blue", "grey", "red"),
                       labels = c(paste("Downregulated (", n_down, ")", sep = ""), "Non-significant", paste("Upregulated (", n_up, ")", sep = ""))) +
    geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
    labs(x = expression(log[2]~(FC)),
         y = expression(-log[10]~(p-value)),
         title = plot_title,
         caption = expression("p-value < 0.05    " ~ log[2]~"(FC) > 0")) +
    theme_bw() +
    theme(text = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5), 
          legend.position = "bottom", 
          legend.title = element_blank())
}

# plots:
my_volcano_plot(dataframe_to_use = db7vsdbmock_deg, 
                plot_title = expression(italic(db/db)~ " Day 7 vs " ~ italic(db/db)~ " Mock"), 
                repel_force = 0.2, 
                repel_max_overlaps = 10)

my_volcano_plot(dataframe_to_use = db7vsWT7_deg, 
                plot_title = expression(italic(db/db)~ "Day 7 vs WT Day 7"), 
                repel_force = 0.2, 
                repel_max_overlaps = 10)

my_volcano_plot(dataframe_to_use = dbmockvsWTmock_deg, 
                plot_title = expression(italic(db/db)~ "Mock vs WT Mock"), 
                repel_force = 0.2, 
                repel_max_overlaps = 10)

my_volcano_plot(dataframe_to_use = WT7vsWTmock_deg, 
                plot_title = expression("WT Day 7 vs WT Mock"), 
                repel_force = 0.2, 
                repel_max_overlaps = 10)
### a library for volcano plot
install.packages('EnhancedVolcano')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
?EnhancedVolcano
EnhancedVolcano(WT7vsWTmock_deg, x="log2FoldChange", y="padj", lab= WT7vsWTmock_deg$gene_name, 
                pCutoff = 0.05, FCcutoff = 1)


###FC vs FC:
##########WT day 7 - WT mock vs dbdb day 7 vs dbdb mock
### Make a dataframe with shared genes in WT7 vs WTmock AND db7 vs db mock
shared_genes <-intersect(WT7vsWTmock_deg$gene_id,db7vsdbmock_deg$gene_id)

shared_WT7_mock <- WT7_WTmock %>% 
  filter(gene_id %in% shared_genes) %>%
  arrange(gene_id)

colnames(shared_WT7_mock) <- c('gene_id (Wt)', 'log2foldchange(WT)', 'pvalue(WT)', 'padj(WT)', 'gene_name(WT)', 'category(WT)')

shared_db7_mock <- db7_mock %>% 
  filter(gene_id %in% shared_genes) %>% 
  arrange(gene_id)

shared_db7_mock[,7] <- NULL

colnames(shared_db7_mock) <- c('gene_id (db)', 'log2foldchange(db)', 'pvalue(db)', 'padj(db)', 'gene_name(db)', 'category(db)')

FC_WT_db <- cbind(shared_WT7_mock, shared_db7_mock)

cat <- c()
for(i in 1:27336){
  if(FC_WT_db$`category(WT)`[i] != "Non-significant" & FC_WT_db$`category(db)`[i] != "Non-significant"){
    cat[i] <- "Sig in both (3665)" #Significantly changed in both
  } 
  else if(FC_WT_db$`category(WT)`[i] == "Non-significant" & FC_WT_db$`category(db)`[i] != "Non-significant"){
    cat[i] <- "Sig only in db/db Day 7 (2886)" #Significant only in dbdb
  } 
  else if(FC_WT_db$`category(WT)`[i] != "Non-significant" & FC_WT_db$`category(db)`[i] == "Non-significant"){
    cat[i] <- "Sig only in WT Day 7 (4501)" #Significant only in WT
  } 
  else if(FC_WT_db$`category(WT)`[i] == "Non-significant" & FC_WT_db$`category(db)`[i] == "Non-significant"){
    cat[i] <- "non-Sig in both (16284)" #Non-significant 
  }
}

FC_WT_db <- FC_WT_db %>% 
  mutate(Category= cat)

FC_WT_db$Category <- factor(FC_WT_db$Category, 
                            levels = c("Sig in both (3665)", "Sig only in db/db Day 7 (2886)", "Sig only in WT Day 7 (4501)", "non-Sig in both (16284)"))

FC_WT_db %>% 
  filter(Category == "Sig only in db/db Day 7 (2886)") %>% 
  count()

library(RColorBrewer)
display.brewer.pal(n = 8, name = 'Set1')
brewer.pal(n = 8, name = "Set1")

FC_WT_db %>% 
  filter(Category == 'Sig in both (3665)' &
           `log2foldchange(WT)` < 0 &
           `log2foldchange(db)` > 0) %>% 
  count()

ggplot(FC_WT_db, 
       aes(x= FC_WT_db$`log2foldchange(WT)`,
           y= FC_WT_db$`log2foldchange(db)`,
           color= FC_WT_db$Category)) +
  scale_color_manual(values = c("#4DAF4A", "#984EA3", "#377EB8", 'NA')) +
  geom_point(aes(alpha= FC_WT_db$Category), size=3) + 
  scale_alpha_manual(values = c(0.85, 0.15, 0.15, 0))+
  geom_vline(xintercept= 0,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept= 0,lty=4,col="black",lwd=0.8) +
  labs(x=expression(~log[2]~"FC( WT Day 7 vs WT Mock)"),
       y=expression(~log[2]~"FC("~italic(db/db)~" Day 7"~ " vs " ~italic(db/db)~ " Mock)"),
       title="") +
  theme_bw()+
  xlim(-4.5,4.5)+ 
  ylim(-6,9)+
  theme(text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14), 
        legend.position="bottom", 
        legend.title = element_blank()) +
  geom_text(font = "bold",size= 6, color= '#4DAF4A', x=-3.75, y=7.5, label="151") +
  geom_text(font = "bold",size= 6, color= '#4DAF4A', x=3.75, y=7.5, label="1579") +
  geom_text(font = "bold",size= 6, color= '#4DAF4A', x=-3.75, y=-5, label="1853") +
  geom_text(font = "bold",size= 6, color= '#4DAF4A', x=3.75, y=-5, label="82")
