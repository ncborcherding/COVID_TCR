library(biomaRt)

add.sex.chromosome.counts <- function(seuratObj) {
  
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Retrieve the list of genes on the Y chromosome
  y_chromosome_genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                              filters = 'chromosome_name',
                              values = 'Y',
                              mart = ensembl)
  y_chromosome_genes <- unique(y_chromosome_genes$external_gene_name)
  y_chromosome_genes <- y_chromosome_genes[y_chromosome_genes != ""]
  
  x_chromosome_genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                              filters = 'chromosome_name',
                              values = 'X',
                              mart = ensembl)
  x_chromosome_genes <- unique(x_chromosome_genes$external_gene_name)
  x_chromosome_genes <- x_chromosome_genes[x_chromosome_genes != ""]
  
  count.matrix <- seuratObj@assays$RNA@counts
  x.pos <- which(rownames(count.matrix) %in% x_chromosome_genes)
  y.pos <- which(rownames(count.matrix) %in% y_chromosome_genes)
  
  
  seuratObj$nCount_ChrY_RNA <- colSums2(seuratObj@assays$RNA@counts[y.pos,])
  seuratObj$nFeature_ChrY_RNA <- colSums(seuratObj@assays$RNA@counts[y.pos,] != 0)
  
  seuratObj$nCount_ChrX_RNA <- colSums2(seuratObj@assays$RNA@counts[x.pos,])
  seuratObj$nFeature_ChrX_RNA <- colSums(seuratObj@assays$RNA@counts[x.pos,] != 0)
  
  return(seuratObj)
}


library(patchwork)

ChrX.features <- FeaturePlot(Figure7, "nFeature_ChrX_RNA") + 
                    scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) + 
                    theme(plot.title = element_blank())
ggsave("~/Documents/ChrX.features.png", height =4, width = 4.5, dpi = 600)
        
ChrY.features <- FeaturePlot(Figure7, "nFeature_ChrY_RNA") + 
                    scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) + 
                    theme(plot.title = element_blank())
ggsave("~/Documents/Chr4.features.png", height =4, width = 4.5, dpi = 600)

Figure7$relative_nCount_ChrX <- Figure7$nCount_ChrX_RNA/Figure7$nCount_RNA*100
Figure7$relative_nCount_ChrY <- Figure7$nCount_ChrY_RNA/Figure7$nCount_RNA*100

ChrY.percent <- FeaturePlot(Figure7, "relative_nCount_ChrY") + 
                    scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) + 
                    theme(plot.title = element_blank())
ggsave("~/Documents/ChrY.percent.png", height =4, width = 4.5, dpi = 600)

ChrX.percent <- FeaturePlot(Figure7, "relative_nCount_ChrX") + 
                    scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu"))) + 
                    theme(plot.title = element_blank())
ggsave("~/Documents/ChrX.percent.png", height =4, width = 4.5, dpi = 600)

meta.data <- Figure7@meta.data[,c("seurat_clusters", "donor", "nFeature_ChrX_RNA", "nFeature_ChrY_RNA", "relative_nCount_ChrY", "relative_nCount_ChrX")]
meta.data[,3:6] <- sapply(meta.data[,3:6], scale)
meta.data <- meta.data[order(meta.data$seurat_clusters),]

my_sample_col <- data.frame(cluster = meta.data$seurat_clusters,
                            donor = meta.data$donor)
row.names(my_sample_col) <- rownames(meta.data)

mycolors <- colorRampPalette(brewer.pal(11, "Paired"))(9)
names(mycolors) <- 0:8

donor.pal <- viridis::viridis_pal(option = "H")(10)
names(donor.pal) <- unique(sort(meta.data$donor))

my_colour = list(cluster = mycolors,
                 donor = donor.pal)

pheatmap::pheatmap(t(meta.data[,3:6]), 
                   show_colnames = FALSE, 
                   annotation_col = my_sample_col, 
                   annotation_colors = my_colour, 
                   cluster_rows = FALSE)

cor(Figure7@meta.data$Frequency, Figure7$nFeature_ChrY_RNA)

