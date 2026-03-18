library(Seurat)
library(dplyr)
library(tidyverse)
library(future)
library(ggplot2)
library(arrow)
library(sf)
library(writexl)
library(grid)
library(patchwork)

packageVersion("Seurat")

xenium.obj <- LoadXenium(
  "INSERT FOLDER PATH",
  assay = "Xenium"
)

xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

# QC plots
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

# Filtering thresholds
xenium.obj <- subset(
  xenium.obj, 
  subset = nFeature_Xenium > 10 & nFeature_Xenium < 480 & 
    nCount_Xenium > 10 & nCount_Xenium < 3000
)

VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

# ---------------------------
# SCTransform on merged object for clustering (as in Methods)
# ---------------------------
num_pcs <- 15
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", verbose = FALSE)
num_features <- nrow(xenium.obj)
dims_to_use <- 1:num_pcs

xenium.obj <- RunPCA(xenium.obj, npcs = num_pcs, features = rownames(xenium.obj), verbose = FALSE)
xenium.obj <- RunUMAP(xenium.obj, dims = dims_to_use)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = dims_to_use)
xenium.obj <- FindClusters(xenium.obj, resolution = 4.3)

DimPlot(xenium.obj, group.by = "seurat_clusters", label = TRUE)

#----------------------------
# Cluster Annotation
#----------------------------

custom_gene_sets <- list(
 
  Basal_epidermal = c("POSTN", "COL17A1", "DST"),
  Non_epidermal_keratinocytes = c("FABP5", "SOSTDC1", "TNC", "PTN", "GJB2"),
  Suprabasal_keratinocytes = c("DSG1", "DSC1", "KRT2", "GATA3"),

  Melanocytes = c("MLANA", "TYR", "PMEL", "TRPM1"),

  Sebocytes = c("FADS2", "DHCR7", "INSIG1"),
  Eccrine_secretory_duct = c("DCD", "KRT8", "KRT19"),
  Endothelial_cells = c("AQP1", "VWF", "SPARCL1", "CD93"),
  Mural = c("ACTA2", "MYL9", "MYH11"),
  
  Fibroblasts = c("COL1A1", "FBN1", "FN1", "COL6A1"),

  Immune0 = c("LYZ", "TYROBP", "ITGAX", "AIF1"),
  Immune1 = c("LYZ", "CD68", "CD163", "TYROBP"),
  Immune2 = c("ITGAX", "CPVL", "CLEC9A", "IRF8", "WDFY4"),
  Immune3 = c("ITGAX", "CPVL", "CLEC10A", "IRF4", "CD1C"),
  Immune4 = c("CD207", "CD1A", "FCER1A"),
  Immune5 = c("C1QA", "C1QB", "CD14"),
  Immune6 = c("TPSAB1", "CPA3", "GATA2", "IL1R1", "TPSB2"),
  Immune7 = c("CD3D", "CD3E", "CD2", "CD52"),
  Immune8 = c("NKG7", "GNLY", "NCR1", "NCAM1"),
  Immune9 = c("CD79A", "MZB1", "TNFRSF17", "CD19"),
  Immune10 = c("LTB", "LCK", "PTPRCAP"),
  Immune11 = c("TNFRSF1B", "TGFB1", "COTL1", "CTLA4")

)

filtered_sets <- lapply(custom_gene_sets, function(g) g[g %in% rownames(xenium.obj)])
filtered_sets <- filtered_sets[sapply(filtered_sets, length) > 0]

# -----------------------------
# Add module scores
# -----------------------------
for (ct in names(filtered_sets)) {
  xenium.obj <- AddModuleScore(xenium.obj, features = list(filtered_sets[[ct]]), name = ct, ctrl = 5)
}

score_names <- paste0(names(filtered_sets), "1")

# -----------------------------
# Cluster annotation by max module score
# -----------------------------
cluster_scores <- lapply(score_names, function(s) {
  xenium.obj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarize(mean_score = mean(.data[[s]], na.rm = TRUE)) %>%
    mutate(score_name = s)
}) %>% bind_rows()

auto_ann <- cluster_scores %>%
  group_by(seurat_clusters) %>%
  slice_max(mean_score, n = 1) %>%
  ungroup() %>%
  mutate(score_name = gsub("1$", "", score_name)) %>%
  select(seurat_clusters, score_name) %>%
  deframe()

xenium.obj <- RenameIdents(xenium.obj, auto_ann)
xenium.obj$CellType <- Idents(xenium.obj)

xenium.obj$CellType <- as.character(xenium.obj$CellType)

#Combine Immune clusters
immune_clusters <- grep("^Immune", xenium.obj$CellType, value = TRUE)
xenium.obj$CellType[xenium.obj$CellType %in% immune_clusters] <- "Immune"
xenium.obj$CellType <- factor(xenium.obj$CellType)
Idents(xenium.obj) <- (xenium.obj$CellType)

# -----------------------------
# Add custom cell cluster colours
# -----------------------------
custom_colors <- c(
  "Basal_epidermal" = "#6495ED",
  "Suprabasal_keratinocytes" = "#00ff7f",
  "Non_epidermal_keratinocytes" = "#CC79A7",

  "Melanocytes" = "#FF1493", 
  
  "Fibroblasts" = "#228B22",
    
  "Sebocytes" = "#BA55D3",	
  "Eccrine_secretory_duct" = "#cd5c5c",
  "Endothelial_cells" = "#9400D3",
  "Mural" = "#87CEEB",

  "Immune" = "#FFD700"
)

DimPlot(
  xenium.obj,
  reduction = "umap",
  label = FALSE,
  label.size = 4,
  cols = custom_colors
) + theme(legend.position = "right")


# -----------------------------
# Hightlight cell clusters spatially
# -----------------------------

seg_file <- "INSERT cell_boundaries.parquet PATH"

seg_data <- read_parquet(seg_file)

seg_sf <- seg_data %>%
  group_by(cell_id, label_id) %>%
  summarise(
    geometry = st_sfc(
      st_polygon(list(
        rbind(
          cbind(vertex_x, vertex_y),
          cbind(vertex_x[1], vertex_y[1])
        )
      ))
    ),
    .groups = "drop"
  ) %>%
  st_as_sf()

# Add cell_id column based on rownames
cell_annotations <- xenium.obj@meta.data %>%
  mutate(cell_id = rownames(.)) %>%
  select(cell_id, CellType)

# Now join with the polygons
seg_annot <- seg_sf %>%
  left_join(cell_annotations, by = "cell_id")

# Create robust fill_group for Immune highlighting
seg_annot <- seg_annot %>%
  mutate(
    fill_group = ifelse(grepl("Immune", CellType, ignore.case = TRUE),
                        "Immune", "Other")
  )


# Remove cells with NA CellType
seg_annot_filtered <- seg_annot %>%
  filter(!is.na(CellType))

# Make CellType a factor with the same order as custom_colors
seg_annot_filtered$CellType <- factor(seg_annot_filtered$CellType, levels = names(custom_colors))

# Plot all cell types with custom colors
p_all <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = CellType), color = NA) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  labs(title = "All Cell Types") +
  guides(fill = guide_legend(ncol = 1))

# Display plot
#p_all


# -----------------------------
# Basal Keratinocytes
# -----------------------------
seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Basal_epidermal",
                        "Basal_epidermal", "Other")
  )

Basal <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Basal_epidermal" = "#6495ED", "Other" = "grey85")) +
  theme_void() +
  labs(title = "Basal epidermal (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Basal

# -----------------------------
# Suprabasal Keratinocytes
# -----------------------------
seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Suprabasal_keratinocytes",
                        "Suprabasal_keratinocytes", "Other")
  )

Suprabasal_keratinocytes <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Suprabasal_keratinocytes" = "#00ff7f", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Suprabasal_keratinocytes (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Suprabasal_keratinocytes

# -----------------------------
# Non_epidermal Keratinocytes
# -----------------------------

seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Non_epidermal_keratinocytes",
                        "Non_epidermal_keratinocytes", "Other")
  )

Non_epidermal_keratinocytes <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Non_epidermal_keratinocytes" = "#CC79A7", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Non_epidermal_keratinocytes (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Non_epidermal_keratinocytes

# -----------------------------
# Melanocytes
# -----------------------------

seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Melanocytes",
                        "Melanocytes", "Other")
  )

Melanocytes <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Melanocytes" = "#FF1493", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Melanocytes (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Melanocytes

# -----------------------------
# Sebocytes 
# -----------------------------

seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Sebocytes",
                        "Sebocytes", "Other")
  )

Sebocytes <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Sebocytes" = "#BA55D3", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Sebocytes (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Sebocytes

# -----------------------------
# Eccrine Secretory Duct Cells
# -----------------------------

seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Eccrine_secretory_duct",
                        "Eccrine_secretory_duct", "Other")
  )

Eccrine_secretory_duct <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Eccrine_secretory_duct" = "#cd5c5c", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Eccrine_secretory_duct (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Eccrine_secretory_duct


# -----------------------------
# Endothelial
# -----------------------------

seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Endothelial_cells",
                        "Endothelial_cells", "Other")
  )

Endothelial_cells <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Endothelial_cells" = "#9400D3", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Endothelial_cells (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Endothelial_cells


# -----------------------------
# Mural
# -----------------------------

seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Mural",
                        "Mural", "Other")
  )

Mural <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Mural" = "#87CEEB", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Mural (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Mural

# -----------------------------
# Fibroblasts
# -----------------------------

seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Fibroblasts",
                        "Fibroblasts", "Other")
  )

Fibroblasts <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Fibroblasts" = "#228B22", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Fibroblasts (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Fibroblasts

# -----------------------------
# Subsetting of Fibroblasts
# -----------------------------

fibro <- subset(xenium.obj, idents = "Fibroblasts")
DefaultAssay(fibro) <- "Xenium"

fibro <- SCTransform(fibro, assay = "Xenium", verbose = FALSE)
fibro <- RunPCA(fibro, npcs = num_pcs, verbose = FALSE)
fibro <- FindNeighbors(fibro, dims = 1:num_pcs)
fibro <- FindClusters(fibro, resolution = 0.8)
fibro <- RunUMAP(fibro, dims = 1:num_pcs)

DimPlot(fibro, label = TRUE) + ggtitle("Fibroblast Subclusters")

fibro_gene_sets <- list(
  Fb1       = c("COL18A1", "COL23A1", "APCDD1", "PTGDS"),
  Fb2       = c("MMP2", "COL12A1", "THBS2", "MFAP5"))

# Keep only genes present in the dataset
fibro_gene_sets_present <- lapply(fibro_gene_sets, function(x) x[x %in% rownames(fibro)])

expr_mat <- GetAssayData(fibro, layer = "data")

for (gene_set_name in names(fibro_gene_sets_present)) {
  genes <- fibro_gene_sets_present[[gene_set_name]]
  score_vector <- colMeans(expr_mat[genes, , drop = FALSE])
  fibro@meta.data[[paste0(gene_set_name, "_score")]] <- score_vector[colnames(fibro)]
}

score_cols <- paste0(names(fibro_gene_sets_present), "_score")
fibro$subtype_annotation <- apply(fibro@meta.data[, score_cols], 1, function(x) names(fibro_gene_sets_present)[which.max(x)])

cluster_subtypes <- fibro@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(dominant_subtype = names(sort(table(subtype_annotation), decreasing = TRUE))[1])

# Map cluster-level dominant subtype back to cells
fibro$cluster_subtype <- cluster_subtypes$dominant_subtype[match(fibro$seurat_clusters, cluster_subtypes$seurat_clusters)]

#Plotting fibroblast subcluster
DimPlot(fibro, group.by = "cluster_subtype", label = TRUE) +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank()
  )

xenium.obj$Fibroblast_Subtype_Cluster <- NA
xenium.obj$Fibroblast_Subtype_Cluster[Cells(fibro)] <- fibro$cluster_subtype

#---------------------------------------
# Spatial locaiton of fibroblast subsets
#---------------------------------------

cell_annotations <- xenium.obj@meta.data %>%
  mutate(cell_id = rownames(.)) %>%
  select(cell_id, Fibroblast_Subtype_Cluster)

# Join polygons + subtype data
seg_annot <- seg_sf %>%
  left_join(cell_annotations, by = "cell_id")

# Create final plotting group
seg_annot <- seg_annot %>%
  mutate(
    fibro_subtype = case_when(
      Fibroblast_Subtype_Cluster == "Fb1"       ~ "Fb1",
      Fibroblast_Subtype_Cluster == "Fb2"       ~ "Fb2",
      TRUE ~ "Non-fibroblast"
    )
  )

# Spatial plot
p_fibro_subtypes <- ggplot(seg_annot) +
  geom_sf(aes(fill = fibro_subtype), color = NA) +
  scale_y_reverse() +
  scale_fill_manual(values = c(
    "Fb1"       = "green",
    "Fb2"       = "red",
    "Non-fibroblast"  = "grey80"
  )) +
  theme_void() +
  labs(title = "Fibroblast Subtypes in Spatial Context") +
  guides(fill = guide_legend(ncol = 1))

p_fibro_subtypes

# -----------------------------
# Immune
# -----------------------------

seg_annot_filtered <- seg_annot_filtered %>%
  mutate(
    fill_basal = ifelse(CellType == "Immune",
                        "Immune", "Other")
  )

Immune <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = fill_basal), color = NA) +
  scale_fill_manual(values = c("Immune" = "#FFD700", "Other" = "grey85")) +
  theme_void() +
  labs(title = " Immune (highlighted)") +
  guides(fill = guide_legend(ncol = 1))

Immune
