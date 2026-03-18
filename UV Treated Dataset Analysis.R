#----------------------------------------
# Load Packages
#----------------------------------------
library(Seurat)
library(dplyr)
library(tidyverse)
library(future)
library(ggplot2)
library(arrow)   
library(sf)      
library(writexl)
library(grid)   

#----------------------------------------
# Load Xenium dataset
#----------------------------------------
xenium.obj <- LoadXenium(
  "C:/Users/mikek/OneDrive - The University of Manchester/Xenium/Custom pan V2 OCT 2025/20251003__115535__BPT_QNU66D_03102025/output-XETG00416__0077735__SITE_C__20251003__115655",
  fov = "fov", assay = "Xenium"
)

# Filter empty cells
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

# QC plots
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

# Filtering thresholds
xenium.obj <- subset(
  xenium.obj, 
  subset = nFeature_Xenium > 5 & nFeature_Xenium < 480 & 
    nCount_Xenium > 10 & nCount_Xenium < 3000
)

VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

#----------------------------------------
# Normalization, PCA, UMAP, clustering
#----------------------------------------
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", verbose = FALSE)

num_features <- nrow(xenium.obj)
num_pcs <- min(15, num_features)
dims_to_use <- 1:num_pcs

xenium.obj <- RunPCA(xenium.obj, npcs = num_pcs, features = rownames(xenium.obj), verbose = FALSE)
xenium.obj <- RunUMAP(xenium.obj, dims = dims_to_use)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = dims_to_use)
xenium.obj <- FindClusters(xenium.obj, resolution = 1.0)
DimPlot(xenium.obj, group.by = "seurat_clusters", label = TRUE)

ElbowPlot(xenium.obj, ndims = 50)

#----------------------------------------
# Cell annotation list
#----------------------------------------

custom_gene_sets <- list(
  Melanocytes = c("MLANA", "TYR", "PMEL", "TRPM1"),
  Basal_epidermal = c("POSTN", "COL17A1", "DST"),
  Non_basal_epidermal = c("DSG1", "DSC1", "KRT2", "GATA3"),
  Non_epidermal_keratinocytes = c("FABP5", "SOSTDC1", "TNC", "PTN", "GJB2"),
  Sebaceuos_gland = c("FADS2", "DHCR7", "INSIG1"),
  Eccrine_gland = c("DCD", "KRT8", "KRT19"),
  Endothelial_cells = c("AQP1", "VWF", "SPARCL1", "CD93"),
  Mural = c("ACTA2", "MYL9", "MYH11"),
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
  Immune11 = c("TNFRSF1B", "TGFB1", "COTL1", "CTLA4"),
  Immune12 = c("MPO", "ELANE"),
  Fibroblasts = c("COL1A1", "FBN1", "FN1", "COL6A1")
)

#----------------------------------------
# Cluster annotation
#----------------------------------------
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
# 3) Cluster annotation by max module score
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
immune_clusters <- grep("^Immune", xenium.obj$CellType, value = TRUE)
xenium.obj$CellType[xenium.obj$CellType %in% immune_clusters] <- "Immune"
xenium.obj$CellType <- factor(xenium.obj$CellType)
Idents(xenium.obj) <- (xenium.obj$CellType)

DimPlot(xenium.obj, reduction = "umap", label = TRUE, label.size = 4)

#----------------------------------------
# Colour palette change
#----------------------------------------
levels(Idents(xenium.obj))

custom_colors <- c(
  
  "Basal_epidermal" = "#6495ED",
  "Non_basal_epidermal" = "#00ff7f",
  "Non_epidermal_keratinocytes" = "#CC79A7",
  "Melanocytes" = "#FF1493", 
  "Fibroblasts" = "#228B22",
  "Sebaceuos_gland" = "#BA55D3",	
  "Eccrine_gland" = "#cd5c5c",
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

#----------------------------------------
# Cell proportions
#----------------------------------------

cell_type_counts <- xenium.obj@meta.data %>%
  group_by(CellType) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

ggplot(cell_type_counts, aes(x = "All Cells", y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

write_xlsx(cell_type_counts, ".xlsx") # change filename here


#----------------------------------------
# colour change
#----------------------------------------
custom_colors <- c(
  "Basal_epidermal" = "#6495ED",
  "Non_basal_epidermal" = "#00ff7f",
  "Non_epidermal_keratinocytes" = "#CC79A7",
  "Melanocytes" = "#FF1493", 
  "Fibroblasts" = "#228B22",
  "Sebaceuos_gland" = "#BA55D3",	
  "Eccrine_gland" = "#cd5c5c",
  "Endothelial_cells" = "#9400D3",
  "Mural" = "#87CEEB",
  "Immune" = "#FFD700"
)

# Set factor levels to match custom_colors order
cell_type_counts$CellType <- factor(cell_type_counts$CellType, levels = names(custom_colors))

ggplot(cell_type_counts, aes(x = "All Cells", y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = custom_colors) +   # custom colors
  scale_y_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),                         
    legend.position = "right",                            
    legend.justification = "center",                      
    legend.margin = margin(l = -15, unit = "pt")          
  ) +
  guides(fill = guide_legend(keyheight = unit(1.0, "cm")))


#----------------------------------------------------
# Load boundaires parquet file for spatial assessment
#----------------------------------------------------

seg_file <- "C:/Users/mikek/OneDrive - The University of Manchester/Xenium/Custom pan V2 OCT 2025/20251003__115535__BPT_QNU66D_03102025/output-XETG00416__0077735__SITE_C__20251003__115655/cell_boundaries.parquet"

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

cell_annotations <- xenium.obj@meta.data %>%
  mutate(cell_id = rownames(.)) %>%
  select(cell_id, CellType)

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

# -------------------------------------------
# Plotting spatial image with annotated cells
# -------------------------------------------
p_all <- ggplot(seg_annot_filtered) +
  geom_sf(aes(fill = CellType), color = NA) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  labs(title = "All Cell Types") +
  guides(fill = guide_legend(ncol = 1))

p_all

#Saving plot
setwd("") # change working directory here

ggsave(".png", p_all, width = 8, height = 6, dpi = 600) # change title here

#--------------------------------
# Spatial of immune cells
#--------------------------------

seg_annot <- seg_annot %>%
  mutate(
    fill_group = ifelse(grepl("Immune", CellType, ignore.case = TRUE),
                        "Immune", "Other")
  )

# (B) Highlight Immune in yellow, others in grey
p_Immune <- ggplot(seg_annot) +
  geom_sf(aes(fill = fill_group), color = NA) +
  scale_fill_manual(values = c("Other" = "grey80", "Immune" = "#FFD700")) +
  theme_void() +
  labs(title = "Immune Highlighted")

# Display Immune plot
p_Immune

#save
ggsave(".png", p_Immune, width = 8, height = 6, dpi = 600)

#--------------------------------
#Immune subset
#--------------------------------

Immune_subset <- subset(xenium.obj, idents = "Immune")
DefaultAssay(Immune_subset) <- "Xenium"
Immune_subset <- SCTransform(Immune_subset, assay = "Xenium", verbose = FALSE)
num_pcs <- min(15, nrow(Immune_subset))
Immune_subset <- RunPCA(Immune_subset, npcs = num_pcs, verbose = FALSE)
Immune_subset <- FindNeighbors(Immune_subset, dims = 1:num_pcs)
Immune_subset <- FindClusters(Immune_subset, resolution = 1.0)
Immune_subset <- RunUMAP(Immune_subset, dims = 1:num_pcs)

DimPlot(Immune_subset, label = TRUE) + ggtitle("Immune Subclusters")

Immune_gene_sets <- list(
  Macrophage = c("CD163", "C1QA", "C1QB", "CD14", "C5AR1"),
  Mast_cell = c("TPSAB1", "CPA3", "IL1RL1", "TPSB2"),
  DC1 = c("ITGAX", "CPVL", "CLEC9A", "WDFY4", "XCR1"),
  DC2 = c("ITGAX", "CPVL", "CLEC10A", "CD1C"),
  Langerhans = c("CD207", "CD1A"),
  T_cell = c("CD2", "CD3E", "CD3D", "CD52"),
  NK_cell = c("GNLY", "NKG7", "CTSW", "NCAM1"),
  B_cell = c("CD79A", "MZB1", "TNFRSF17", "CD19"),
  Neutrophil = c("MPO", "ELANE"))

Immune_gene_sets_present <- lapply(Immune_gene_sets, function(x) x[x %in% rownames(Immune_subset)])

expr_mat <- GetAssayData(Immune_subset, slot = "data")

for (gene_set_name in names(Immune_gene_sets_present)) {
  genes <- Immune_gene_sets_present[[gene_set_name]]
  score_vector <- colMeans(expr_mat[genes, , drop = FALSE])
  Immune_subset@meta.data[[paste0(gene_set_name, "_score")]] <- score_vector[colnames(Immune_subset)]
}

score_cols <- paste0(names(Immune_gene_sets_present), "_score")
Immune_subset$subtype_annotation <- apply(Immune_subset@meta.data[, score_cols], 1, function(x) names(Immune_gene_sets_present)[which.max(x)])

cluster_subtypes <- Immune_subset@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(dominant_subtype = names(sort(table(subtype_annotation), decreasing = TRUE))[1])

# Map cluster-level dominant subtype back to cells
Immune_subset$cluster_subtype <- cluster_subtypes$dominant_subtype[match(Immune_subset$seurat_clusters, cluster_subtypes$seurat_clusters)]

DimPlot(Immune_subset, group.by = "subtype_annotation", label = TRUE, repel = TRUE, label.size = 6) +
  ggtitle("Immune Subtype Annotation (Per Cell)") +
  theme(legend.position = "bottom") +
  xlab(NULL) + ylab(NULL)

xenium.obj$Immune_subset <- NA
xenium.obj$Immune_subset[Cells(Immune_subset)] <- Immune$subtype_annotation

#--------------------------------
# Immune subset spatial
#--------------------------------

# Extract immune subtype annotation per cell
cell_annotations <- Immune_subset@meta.data %>%
  mutate(cell_id = rownames(.)) %>%
  select(cell_id, subtype_annotation)

seg_annot <- seg_sf %>%
  left_join(cell_annotations, by = "cell_id")

seg_annot <- seg_annot %>%
  mutate(
    immune_subtype = ifelse(is.na(subtype_annotation),
                            "Non-immune",
                            subtype_annotation)
  )

p_immune_subtypes <- ggplot(seg_annot) +
  geom_sf(aes(fill = immune_subtype), color = NA) +
  scale_y_reverse() +
  scale_fill_manual(values = c(
    "Non-immune"  = "grey80",
    "B_cell"      = "red",
    "DC2"         = "gold",
    "DC1"         = "lightblue",
    "Macrophage"  = "darkblue",
    "T_cell"      = "pink",
    "Neutrophil"  = "yellow",
    "Langerhans"  = "orange",
    "Mast_cell"   = "darkgreen",
    "NK_cell"     = "purple"
  )) +
  theme_void() +
  labs(title = "Immune Cell Subsets in Spatial Context") +
  guides(fill = guide_legend(ncol = 1))

p_immune_subtypes

ggsave(".png",
       p_immune_subtypes,
       width = 8, height = 6, dpi = 600) # change title here

#-------------------------------------------------------
# Single cell type spatial analysis - Langerhans example
#-------------------------------------------------------

Langerhans <- seg_annot %>%
  filter(immune_subtype == "Langerhans")

ggplot(seg_annot) +
  geom_sf(fill = "grey90", color = NA) +  # background
  geom_sf(data = seg_annot %>% filter(immune_subtype == "Langerhans"),
          fill = "red", color = NA) +
  scale_y_reverse() +
  theme_void() +
  labs(title = "Langerhans Highlighted")

ggsave("Langerhans.png",
       width = 8, height = 6, dpi = 600)