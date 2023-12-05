library(biomaRt)
seuratObj <- scRepertoire::scRep_example
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