# Pipeline para análisis de expresión diferencial con DESeq2

### Cargar las librerías necesarias
library("pheatmap")
library("DESeq2")
library("ggplot2")
library("corrplot")
library("RColorBrewer")
library("gplots")
library("vsn")
library("genefilter")
library("EnhancedVolcano")

###  ---- 1. Preparación de los datos ---- #
- `metadata`: tabla con nombres de muestras y sus condiciones (p. ej., control o tratamiento).
- `counts`: matriz de conteos crudos con genes como filas y muestras como columnas.

1) Importar datos
```R
metadata <- read.csv("metadata.csv", row.names = 1)
counts <- read.csv("counts.csv", row.names = 1)
```
2) Verificar integridad de los datos
```R
head(metadata)
head(counts)
stopifnot(all(colnames(counts) == metadata$SampleID))
```

3) Crear objeto DESeqDataSet
```R
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~Condition)
```

### ---- 2. Filtrado y normalización ---- #

1) Filtrar genes con pocos conteos
```R
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]
```

 2)Normalización
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
write.table(normalized_counts, file = "normalized_counts.txt", sep = "\t", quote = FALSE, col.names = NA)

### ---- 3. Análisis de expresión diferencial ---- #
dds <- DESeq(dds)
res <- results(dds)
write.table(res, file = "differential_expression_results.txt", sep = "\t", quote = FALSE, col.names = NA)

### Filtrar resultados por p-valor ajustado
res.05 <- results(dds, alpha = 0.05)
summary(res.05)

### ---- 4. Visualización ---- #
### Gráfico de dispersión
plotDispEsts(dds, main = "Dispersion plot")

### MA-plot
plotMA(res, main = "MA Plot")
plotMA(res.05, main = "MA Plot (p-adj < 0.05)")

### Shrinkage de Log2 Fold Changes
resLFC <- lfcShrink(dds, coef = "Condition_Treatment_vs_Control", type = "apeglm")
plotMA(resLFC, ylim = c(-2, 2), main = "MA Plot with LFC Shrinkage")

### PCA
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = c("Condition"))

### Heatmap de distancias entre muestras
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(metadata$Condition, metadata$SampleID, sep = "-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, col = colors)

### Heatmap de genes
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:50]
pheatmap(assay(vsd)[select, ], cluster_rows = TRUE, scale = "row", show_rownames = FALSE)

### Volcano plot
EnhancedVolcano(resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC),
                xlim = c(-5, 5),
                ylim = c(0, -log10(1e-6)),
                pCutoff = 0.05,
                FCcutoff = 2)

### ---- 5. Guardar resultados y gráficos ---- #
pdf("heatmap_genes.pdf")
pheatmap(assay(vsd)[select, ], cluster_rows = TRUE, scale = "row", show_rownames = FALSE)
dev.off()

pdf("pca_plot.pdf")
plotPCA(vsd, intgroup = c("Condition"))
dev.off()
