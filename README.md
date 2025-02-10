# Análisis de Expresión Diferencial en Leucemia Mieloide Aguda (LAML)

Este repositorio contiene el análisis de expresión diferencial utilizando datos de pacientes con **Leucemia Mieloide Aguda (LAML)** del proyecto TCGA. Se emplean herramientas de **bioinformática** para la selección, limpieza y visualización de datos de expresión génica.

## Contenido

- [Introducción](#introducción)
- [Requisitos](#requisitos)
- [Configuración del Entorno](#configuración-del-entorno)
- [Selección del Set de Datos](#selección-del-set-de-datos)
- [Formateo de la Información del Proyecto](#formateo-de-la-información-del-proyecto)
- [Análisis de Expresión Diferencial](#análisis-de-expresión-diferencial)
- [Visualización de Resultados](#visualización-de-resultados)
  - [Gráfico Tipo Volcán](#gráfico-tipo-volcán)
  - [MA Plot de las Muestras](#ma-plot-de-las-muestras)
  - [Heatmap de Genes de Interés](#heatmap-de-genes-de-interés)
  - [Boxplot de Expresión de un Gen](#boxplot-de-expresión-de-un-gen)
  - [PCA de las Muestras](#pca-de-las-muestras)
- [Conclusión](#conclusión)
- [Instalación](#instalación)
- [Uso](#uso)
- [Agradecimientos](#agradecimientos)

## Introducción

Este documento presenta un análisis de expresión diferencial utilizando datos de pacientes con Leucemia Mieloide Aguda (LAML) del proyecto TCGA. Se emplean herramientas de bioinformática para la selección, limpieza y visualización de datos de expresión génica.

## Requisitos

Las siguientes librerías de R son necesarias para ejecutar el análisis:

- `recount3`
- `ggplot2`
- `DESeq2`
- `pheatmap`
- `EnhancedVolcano`

## Configuración del Entorno

Primero, se deben cargar las librerías necesarias y verificar que se hayan cargado correctamente.

```r
# Cargar librerías necesarias
library(recount3)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)

# Verificar que los paquetes se hayan cargado correctamente
if (all(c("recount3", "ggplot2", "DESeq2", "pheatmap", "EnhancedVolcano") %in% loadedNamespaces())) {
  print("Se cargaron de manera exitosa los paquetes")
} else {
  print("No se pudieron cargar todos los paquetes")
}
```

## Selección del Set de Datos
Se selecciona el conjunto de datos correspondiente a LAML del repositorio recount3, específicamente del consorcio TCGA.

``` r
# Obtener la lista de proyectos disponibles de humanos en recount3
human_projects <- available_projects(organism = "human")

# Filtrar solo los proyectos pertenecientes a TCGA
tcga_projects <- human_projects[human_projects$project_home == "data_sources/tcga", ]

# Seleccionar el proyecto LAML
project_info <- subset(tcga_projects, project == "LAML")

# Crear un objeto RangedSummarizedExperiment con los datos de LAML
rse_LAML <- create_rse(project_info)

# Calcular los conteos de lectura y almacenarlos en el assay "counts"
assay(rse_LAML, "counts") <- compute_read_counts(rse_LAML)

# Mostrar el objeto con información de los datos
print(rse_LAML)
```

## Formateo de la Información del Proyecto
Para un análisis más preciso, se filtran genes de interés y se seleccionan las variables clínicas relevantes del dataset.

``` r
# Definir genes de interés para el análisis
genes_interes <- c("SPN", "RUNX1", "CEBPA", "GATA2", "SPI1", "MYB", "FLI1", "ERG", "MECOM", "TAL1", "LMO2", "LDB1", "CBFB", "GFI1", "HOXA9", "MEIS1", "KMT2A", "WT1", "EZH2", "DNMT3A", "TET2", "ASXL1", "IDH1", "IDH2", "NPM1", "FLT3", "KIT", "CSF3R", "MPL", "JAK2", "STAT5A", "STAT3", "ETV6", "NRAS", "KRAS", "PTPN11", "NF1", "CBL", "GATA1", "GATA3", "ZBTB16", "EVI5", "FOXO3", "BCL2", "BCL6", "BAX", "MCL1", "CDKN1A", "CDKN2A", "TP53", "RB1", "MDM2", "IKZF1", "DNTT", "RAG1", "RAG2", "E2A", "HHEX", "ZNF521", "PRDM16", "ARID5B", "KLF4", "KLF5", "MAFB", "IRF8", "IRF4", "NFE2", "NFE2L2", "BACH1", "BACH2", "EP300", "CREBBP", "CBX5", "SUZ12", "SMARCA4", "SMARCB1", "CTCF", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FOXP1", "FOXP3", "NOTCH1", "NOTCH2", "DLL1", "JAG1", "HES1", "HEY1", "SOCS1", "SOCS3", "PPARG", "NCOR1", "NCOR2", "RXRA", "VDR", "MBD2", "TGFBR1", "SMAD3", "SMAD4", "CD3D", "MYC", "RELA", "NFKB1", "NFKB2", "BCOR", "CD3E", "CD3G", "CEBPB", "CEBPD", "CEBPG", "STAT2", "STAT4", "STAT6", "SOCS2", "SOCS4", "SOCS5", "SOCS6", "SOCS7", "SMAD1", "SMAD2", "SMAD5", "SMAD6", "SMAD7", "TGFB1", "TGFB2", "TGFB3", "TNF", "TNFRSF1A", "TNFRSF1B", "TNFAIP3", "TNIP1", "BIRC2", "BIRC3", "BIRC5", "XIAP", "FAS", "FASLG", "TRAF1", "TRAF2", "TRAF3", "TRAF6", "NLRP3", "NLRP1", "CASP1", "CASP3", "CASP7", "CASP8", "CASP9", "CASP10", "BAK1", "BID", "BAD", "BBC3", "MALT1", "CARD11", "CARD9", "NOD1", "NOD2", "MYD88", "TICAM1", "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "TLR10", "DOK1", "DOK2", "DOK3", "DOK4", "DOK5", "DOK6", "SH2B1", "SH2B2", "SH2B3", "CBLB", "CBL2", "UBASH3A", "UBASH3B", "LCP2", "LAT", "FYB", "GRAP", "GRB2", "GAB1", "GAB2", "GAB3", "SHC1", "SHC2", "SHC3", "SHC4", "CRKL", "CRK", "NCK1", "NCK2", "VAV1", "VAV2", "VAV3", "DOCK2", "DOCK8", "ITK", "BTK", "TXK", "TEC", "LCK", "FYN", "HCK", "LYN", "BLK", "YES1", "SYK", "ZAP70", "CSK", "PTK2", "PTK2B", "FER", "FES", "FGR", "EPHA1", "EPHA2", "EPHA3", "EPHA4", "EPHA5", "EPHA6", "EPHA7", "EPHA8", "EPHB1", "EPHB2", "EPHB3", "EPHB4", "EPHB6", "KITLG", "FLT1", "FLT4", "KDR", "PDGFRA", "PDGFRB", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "EGFR", "ERBB2", "ERBB3", "ERBB4", "INSR", "IGF1R", "IGF2R", "MET", "RON", "AXL", "MERTK", "TYRO3", "TEK", "TIE1", "ROR1", "ROR2", "ALK", "ROS1", "NTRK1", "NTRK2", "NTRK3", "DDR1", "DDR2", "EPHA10", "EPHB10","STK11", "MTOR", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "AKT1", "AKT2", "AKT3", "PTEN", "PDPK1", "PDPK2", "RAC1", "RAC2", "RAC3", "RHOA", "RHOB", "RHOC", "CDC42", "ARHGEF1", "ARHGEF2", "ARHGEF3", "ARHGEF4", "ARHGEF5", "ARHGEF6", "ARHGEF7", "ARHGEF8", "ARHGEF9", "ARHGEF10", "ARHGEF11", "DOCK1", "DOCK3", "DOCK4", "DOCK5", "DOCK6", "DOCK7", "DOCK9", "DOCK10", "DOCK11", "DOCK12", "DOCK13", "DOCK14", "DOCK15", "DOCK16", "DOCK17", "DOCK18", "DOCK19", "DOCK20", "DOCK21", "DOCK22", "DOCK23", "DOCK24", "DOCK25", "DOCK26", "DOCK27")

# Filtrar solo los genes de interés en el conjunto de datos
rse_LAML2 <- rse_LAML[which(rowData(rse_LAML)$gene_name %in% genes_interes), ]

# Filtrar columnas relevantes en la metadata del proyecto
columnas_interes <- c("tcga.gdc_cases.diagnoses.age_at_diagnosis",
                      "tcga.gdc_cases.diagnoses.vital_status",
                      "tcga.gdc_cases.samples.sample_type")
colData(rse_LAML2) <- colData(rse_LAML)[, columnas_interes]
```

## Análisis de Expresión Diferencial
Se lleva a cabo un análisis de expresión diferencial utilizando DESeq2, que permite identificar genes diferencialmente expresados en función del estado vital de los pacientes.

``` r
# Asegurar que 'counts' sea la primera matriz en assays
assays(rse_LAML2) <- assays(rse_LAML2)[c("counts", "raw_counts")]

# Filtrar muestras sin NA en vital_status
rse_LAML2 <- rse_LAML2[, !is.na(colData(rse_LAML2)$tcga.gdc_cases.diagnoses.vital_status)]

# Convertir datos a objeto DESeq2
dds <- DESeqDataSet(rse_LAML2, design = ~ tcga.gdc_cases.diagnoses.vital_status)
dds <- DESeq(dds)  # Ajuste del modelo

# Obtener resultados del análisis diferencial
res <- results(dds)
res$gene <- rowData(rse_LAML2)$gene_name[match(rownames(res), rownames(rse_LAML2))]
res <- res[order(res$padj), ]
```

## Visualización de Resultados
### Gráfico Tipo Volcán
Un gráfico tipo volcán permite visualizar los genes más significativamente diferencialmente expresados en función del log2FoldChange y el valor de p-value.

``` r
EnhancedVolcano(res,
                lab = res$gene,
                x = "log2FoldChange",
                y = "pvalue",
                title = "Expresión diferencial en LAML",
                pCutoff = 0.05)
```

### MA Plot de las Muestras
Un MA plot permite visualizar la expresión diferencial de los genes al comparar la media de expresión normalizada (Eje X) con el log2FoldChange (Eje Y).

``` r
pltMA(res, ylim = c(-5, 5), cex = 1)
```

### Heatmap de Genes de Interés
Se genera un heatmap para visualizar la expresión de los genes de interés en las distintas muestras.

``` r
# Transformación de varianza estabilizada
vsd <- varianceStabilizingTransformation(dds)

# Asegurar que los rownames de vsd sean los nombres de los genes
rownames(vsd) <- rowData(vsd)$gene_name

# Crear la matriz de expresión
vsd_matrix <- assay(vsd)

genes_interes_heatmap <- c("SPN", "RUNX1", "CEBPA", "GATA2", "SPI1", "MYB", "FLI1", "ERG")

# Filtrar solo los genes de interés
genes_presentes <- genes_interes_heatmap[genes_interes_heatmap %in% rownames(vsd_matrix)]
vsd_matrix <- vsd_matrix[genes_presentes, , drop = FALSE]

# Verificar que la matriz tenga suficientes genes antes de hacer el heatmap
if (length(genes_presentes) == 0) {
    stop("Ninguno de los genes de interés está presente en la matriz de expresión.")
} else if (length(genes_presentes) == 1) {
    pheatmap::pheatmap(vsd_matrix, show_colnames = FALSE, cluster_rows = FALSE)
} else {
    pheatmap::pheatmap(vsd_matrix, show_colnames = FALSE)
}
```

### Boxplot de Expresión de un Gen
Se visualiza la expresión del gen SPN por estado vital en un boxplot.

``` r
# Extraer expresión de un gen específico
SPN_expr <- as.vector(assay(vsd)["SPN", ])

df_boxplot <- data.frame(
  vital_status = colData(vsd)$tcga.gdc_cases.diagnoses.vital_status,
  expression = SPN_expr
)

# Generar boxplot
ggplot(df_boxplot, aes(x = vital_status, y = expression)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Expresión de SPN por estado vital")
```

### PCA de las Muestras
Se realiza un Análisis de Componentes Principales (PCA) para explorar la variabilidad entre muestras.

``` r
vsd2 <- varianceStabilizingTransformation(dds, blind = TRUE)
pcaData <- plotPCA(vsd2, intgroup = "tcga.gdc_cases.diagnoses.vital_status", returnData = TRUE)

ggplot(pcaData, aes(PC1, PC2, color = tcga.gdc_cases.diagnoses.vital_status)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA de las muestras")
```

## Conclusión
Este análisis proporciona una visión detallada de la expresión diferencial en Leucemia Mieloide Aguda, ayudando a identificar genes clave que podrían ser relevantes en la progresión de la enfermedad y posibles objetivos terapéuticos.

## Instalación
Para instalar las librerías necesarias, puede utilizar el siguiente comando en R:

``` r
install.packages(c("recount3", "ggplot2", "DESeq2", "pheatmap", "EnhancedVolcano"))
```

## Uso
Para ejecutar el análisis, simplemente ejecute el script analisis_expresion.R en su entorno de R. Asegúrese de que todas las librerías necesarias estén instaladas y cargadas correctamente.

## Agradecimientos
Agradecemos al proyecto TCGA por proporcionar los datos utilizados en este análisis y a los desarrolladores de las librerías de R utilizadas en este estudio.
