# -----------------------------
# Load Packages
# -----------------------------
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
library(extrafont)

# -----------------------------
# Load Xenium files
# -----------------------------
paths <- c(
"INSERT CELL FEATURE MATRIX PATHS"
)

# Read Xenium files
process_xenium_sample <- function(path, sample_name) {
  x <- Read10X(data.dir = path)
  if (is.list(x)) x <- x[[1]]  # Xenium sometimes returns a list
  
  obj <- CreateSeuratObject(
    counts = x,
    assay = "Xenium",
    project = sample_name
  )
  
  # QC: keep cells with ≥10 detected genes
  obj <- subset(obj, subset = nFeature_Xenium >= 10)
  
  obj$sample <- sample_name
  return(obj)
}

sample_names <- c("D1921", "D1923", "D1924", "D1925", "D1922")

# Read + merge all samples in one go
objs <- lapply(seq_along(paths), function(i) {
  s_name <- sample_names[i]
   obj <- process_xenium_sample(paths[i], s_name)
  obj <- RenameCells(obj, add.cell.id = s_name)
  return (obj)
})

all_combined_data <- Reduce(merge, objs)

# ---------------------------
# Preprocessing: SCTransform
# ---------------------------

num_pcs <- 15

all_combined_data <- subset(
  all_combined_data,
  subset =
    nFeature_Xenium > 10 & nFeature_Xenium < 480 &
    nCount_Xenium  > 10 & nCount_Xenium  < 3000
)

all_combined_data <- SCTransform(
  all_combined_data,
  assay = "Xenium",
  layer = "counts",
  verbose = FALSE
)


# ---------------------------
# Preprocessing: PCA
# ---------------------------

all_combined_data <- RunPCA(
  all_combined_data,
  assay = "SCT",
  npcs = num_pcs,
  verbose = FALSE
)

# ---------------------------
# Clustering: UMAP
# ---------------------------

all_combined_data <- RunUMAP(
  all_combined_data,
  dims = 1:num_pcs
)

all_combined_data <- FindNeighbors(
  all_combined_data,
  dims = 1:num_pcs
)

all_combined_data <- FindClusters(
  all_combined_data,
  resolution = 1.0
)

DimPlot(
  all_combined_data,
  reduction = "umap",
  label = TRUE
) + ggtitle("Clusters (SCT-based)")


# ---------------------------
# Elbow Plot
# ---------------------------

ElbowPlot(all_combined_data, ndims = 50)


# ---------------------------
# Cell Annotation List
# ---------------------------

cell_annotation_gene_sets <- list(

  Melanocytes = c("MLANA", "TYR", "PMEL", "TRPM1"),
  Basal_epidermal = c("POSTN", "COL17A1", "DST"),
  Non_epidermal_keratinocytes = c("FABP5", "SOSTDC1", "TNC", "PTN", "GJB2"),
  Suprabasal_epidermal = c("DSG1", "DSC1", "KRT2", "GATA3"),

  Sebocytes = c("FADS2", "DHCR7", "INSIG1"),
  Eccrine_duct_cells = c("DCD", "KRT8", "KRT19"),
  Endothelial_cells = c("AQP1", "VWF", "SPARCL1", "CD93"),
  Mural = c("ACTA2", "MYL9", "MYH11"),

  Fibroblasts = c("COL1A1", "FBN1", "FN1", "COL6A1"),

  Immune0 = c("LYZ", "TYROBP", "ITGAX", "AIF1"),
  Immune1 = c("LYZ", "CD68", "CD163", "TYROBP"),
  Immune2 = c("ITGAX", "CPVL", "CLEC9A", "IRF8", "WDFY4"),
  Immune3 = c("ITGAX", "CPVL", "CLEC10A", "IRF4", "CD1C"),
  Immune4 = c("CD207", "CD1A", "FCER1A"),
  Immune5 = c("C1QA", "C1QB", "CD14", "C5AR1"),
  Immune6 = c("TPSAB1", "CPA3", "GATA2", "IL1R1", "TPSB2"),
  Immune7 = c("CD3D", "CD3E", "CD2", "CD52"),
  Immune8 = c("NKG7", "GNLY", "NCR1", "NCAM1"),
  Immune9 = c("CD79A", "MZB1", "TNFRSF17", "CD19"),
  Immune10 = c("LTB", "LCK", "PTPRCAP"),
  Immune11 = c("TNFRSF1B", "TGFB1", "COTL1", "CTLA4"),
  Immune12 = c("MPO", "ELANE")


)

# -----------------------------
# Filter genes 
# -----------------------------
filtered_sets <- lapply(cell_annotation_gene_sets, function(g) g[g %in% rownames(all_combined_data)])
filtered_sets <- filtered_sets[sapply(filtered_sets, length) > 0]

# -----------------------------
# Add module scores
# -----------------------------
for (ct in names(filtered_sets)) {
  all_combined_data <- AddModuleScore(all_combined_data, features = list(filtered_sets[[ct]]), name = ct, ctrl = 5)
}

score_names <- paste0(names(filtered_sets), "1")

# -----------------------------
# Cluster annotation by max module score
# -----------------------------
cluster_scores <- lapply(score_names, function(s) {
  all_combined_data@meta.data %>%
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


all_combined_data <- RenameIdents(all_combined_data, auto_ann)
all_combined_data$CellType <- Idents(all_combined_data)

all_combined_data$CellType <- as.character(all_combined_data$CellType)


# -----------------------------
# Combine the Immune cells
# -----------------------------

immune_clusters <- grep("^Immune", all_combined_data$CellType, value = TRUE)
all_combined_data$CellType[all_combined_data$CellType %in% immune_clusters] <- "Immune"
all_combined_data$CellType <- factor(all_combined_data$CellType)
Idents(all_combined_data) <- all_combined_data$CellType


# -----------------------------
# Plot Annotated clusters
# -----------------------------
DimPlot(all_combined_data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# ---------------------------
# Assign custom cell cluster colours
# ---------------------------
custom_colors <- c(

  #Keratinocytes
  "Basal_epidermal" = "#6495ED",
  "Suprabasal_epidermal" = "#00ff7f",
  "Non_epidermal_keratinocytes" = "#CC79A7",
  
  #Melanocytes
  "Melanocytes" = "#FF1493", 
  
  #Fibroblasts
  "Fibroblasts" = "#228B22",
  
  #Appendage and vasculature
  "Sebocytes" = "#BA55D3",	
  "Eccrine_duct_cells" = "#cd5c5c",
  "Endothelial_cells" = "#9400D3",
  "Mural" = "#87CEEB",
  
  #Immune cells
  "Immune" = "#FFD700"
)

# ---------------------------
# Plot new UMAP with custom cell cluster colours
# ---------------------------

DimPlot(
  all_combined_data,
  reduction = "umap",
  label = FALSE,
  label.size = 4,
  cols = custom_colors
) + theme(legend.position = "right")

# ---------------------------
# Cell Proportions
# ---------------------------
cell_type_counts <- all_combined_data@meta.data %>%
  group_by(CellType) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

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


#write_xlsx(cell_type_counts, "*INSERT PATH*")

# ---------------------------
# Dotplots
# ---------------------------

all_combined_data <- RenameCells(all_combined_data, add.cell.id = "sample")
Idents(all_combined_data) <- all_combined_data$CellType


# Define marker genes
marker_genes <- c( 
                  # Keratinocytes
                  "COL17A1", "DST", "POSTN", 
                  "DSG1", "DSC1", "KRT2", "GATA3",
                  "FABP5", "SOSTDC1", "TNC", "PTN", "GJB2",
                  
                  # Melanocytes
                  "MLANA", "PMEL", "TYR", "TRPM1",
                  
                  # Appendages & Vasculature
                  "FADS2", "DHCR7", "INSIG1",
                  "KRT18", "KRT8", "DCD",               
                  "ACTA2", "MYL9", "MYH11",
                  "AQP1", "VWF", "SPARCL1","CD93",

                  # Fibroblasts
                  "COL1A1", "FBN1", "FN1", "COL6A1", "LUM",
                  
                  # Immune
                  "LYZ", "ITGAX", "AIF1", "TYROBP", "CD68", "CD14",
                  "CD207", "TPSAB1", 
                  "CD2")


celltype_order <- rev(c(
                        "Basal_epidermal",
                        "Suprabasal_epidermal",
                        "Non_epidermal_keratinocytes",
                      
                        "Melanocytes",
                       
                        "Sebocytes",
                        "Eccrine_duct_cells",
                        "Mural",
                        "Endothelial_cells",

                        "Fibroblasts",
                      
                        "Immune"
                      ))

all_combined_data$CellType <- factor(all_combined_data$CellType, levels = celltype_order)
Idents(all_combined_data) <- all_combined_data$CellType


# All genes Dot plot
DotPlot(all_combined_data, features = marker_genes) +
  RotatedAxis() +
  scale_color_gradientn(colors = c("white", "lightblue", "darkblue")) +
  scale_size(range = c(0,5)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "right"
  )


# ---------------------------
# DotPlot for Epidermal Layers
# ---------------------------

epidermal_celltype_order <- c(
  "Basal_epidermal",
  "Suprabasal_epidermal",
  "Non_epidermal_keratinocytes",
  "Melanocytes"
)

# 2. Subset the object to only these types
plot_subset <- subset(all_combined_data, idents = epidermal_celltype_order)

# 3. Ensure the factor levels match the desired order (reversed for the Y-axis)
plot_subset$CellType <- factor(plot_subset$CellType, levels = rev(epidermal_celltype_order))
Idents(plot_subset) <- plot_subset$CellType

# 4. Define ONLY the genes relevant to these four groups
relevant_genes <- c(
  "COL17A1", "DST", "POSTN",         # Basal
  "DSG1", "DSC1", "KRT2", "GATA3",   # Suprabasal
  "FABP5", "SOSTDC1", "TNC", "PTN",  # Non-epidermal Kerat.
  "MLANA", "PMEL", "TYR", "TRPM1"    # Melanocytes
)

# 5. Generate the Plot
epi <- DotPlot(
  plot_subset, 
  features = relevant_genes, 
  assay = "SCT", 
 # cols = c("white", "lightblue", "darkblue"),
  dot.min = 0.05,
) + scale_color_gradientn(
    colors = c("white", "lightblue", "darkblue"),
    name = "Average expression"
  ) +  RotatedAxis() +
  scale_size(range = c(0, 5), name = "Percentage expression") +
  theme(
    # Base text settings
    text = element_text(family = "Arial", size = 18),
    
    # Axis settings
    axis.text.x = element_text(family = "Arial", size = 18, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    axis.text.y = element_text(family = "Arial", size = 18),
    axis.title = element_text(family = "Arial", size = 18),
    
    # Legend settings
    legend.text = element_text(family = "Arial", size = 10),
    legend.title = element_text(family = "Arial", size = 10),
    legend.key.size = unit(0.3, "cm"), # Smaller legend icons to match 5pt text
    
    # Border
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(color = "Average expression", x = "Markers", y = "Cell Type")

# 6. Save as TIFF,JPG or PNG 
ggsave(
  filename = "Epidermal_Markers_Dot_plot.tiff",
  plot = epi,
  device = "tiff",
  dpi = 600,
  compression = "lzw",
  width = 8.5,
  height = 6,
  units = "cm"
)

# ---------------------------
# DotPlot for Appendages
# ---------------------------

# 1. Define the specific cell types and their order
appendagetype_order <- c(
                        "Sebocytes",
                        "Eccrine_duct_cells",
                        "Mural",
                        "Endothelial_cells"
                        
)

# 2. Subset the object to only these types
plot_subset <- subset(all_combined_data, idents = appendagetype_order)

# 3. Ensure the factor levels match the desired order
plot_subset$CellType <- factor(plot_subset$CellType, levels = rev(appendagetype_order))
Idents(plot_subset) <- plot_subset$CellType

# 4. Define ONLY the genes relevant to these four groups
appendage_relevant_genes <- c(
                 "FADS2", "DHCR7", "INSIG1",
                  "KRT18", "KRT8", "DCD",               
                  "ACTA2", "MYL9", "MYH11",
                  "AQP1", "VWF", "SPARCL1","CD93"
)

# 5. Generate the Plot
app <- DotPlot(
  plot_subset, 
  features = appendage_relevant_genes, 
  assay = "SCT", 
 # cols = c("white", "lightblue", "darkblue"),
  dot.min = 0.05,
) + scale_color_gradientn(
    colors = c("white", "lightblue", "darkblue"),
    name = "Average expression"
  ) +  RotatedAxis() +
  scale_size(range = c(0, 5), name = "Percentage expression") +
  theme(
    # Base text settings
    text = element_text(family = "Arial", size = 18),
    
    # Axis settings
    axis.text.x = element_text(family = "Arial", size = 18, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    axis.text.y = element_text(family = "Arial", size = 18),
    axis.title = element_text(family = "Arial", size = 18),
    
    # Legend settings
    legend.text = element_text(family = "Arial", size = 10),
    legend.title = element_text(family = "Arial", size = 10),
    legend.key.size = unit(0.3, "cm"), # Smaller legend icons to match 5pt text
    
    # Border
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(color = "Average expression", x = "Markers", y = "Cell Type")

# 6. Save as TIFF, JPG or PNG
ggsave(
  filename = "Appendage_Markers_Dot_plot.tiff",
  plot = app,
  device = "tiff",
  dpi = 600,
  compression = "lzw",
  width = 8.5,
  height = 6,
  units = "cm"
)


# ---------------------------
# Keratinocyte Only Cluster Highlighting
# ---------------------------

# Create a plotting label that keeps 4 groups, and collapses everything else to "Other"
all_combined_data$highlight_group <- ifelse(
  all_combined_data$CellType %in% c(
    "Basal_epidermal",
    "Non_epidermal_keratinocytes",
    "Suprabasal_epidermal"
   # "Melanocytes"
  ),
  as.character(all_combined_data$CellType),
  "Other"
)

# Set order so legend is nice and consistent
all_combined_data$highlight_group <- factor(
  all_combined_data$highlight_group,
  levels = c(
    "Basal_epidermal",
    "Non_epidermal_keratinocytes",
    "Suprabasal_epidermal",
    "Other"
  )
)

DimPlot(
  all_combined_data,
  reduction = "umap",
  group.by = "highlight_group"
) +
  scale_color_manual(values = c(
    "Basal_epidermal"             = "#6495ED",
    "Non_epidermal_keratinocytes" = "#CC79A7",
    "Suprabasal_epidermal"        = "#00ff7f",
    "Other"                       = "grey95"
  )) +
  ggtitle("Main Epidermal Cells") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# ---------------------------
# Feature and Dimplots per cell type 
# ---------------------------

# Basal Keratinocytes
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Basal_epidermal", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#6495ED", "no" = "grey95")) +
  ggtitle("Basal_epidermal") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )

Basal_Keratinocytes <-c("POSTN", "COL17A1", "DST")
# Generate FeaturePlots for each gene
bkp <- FeaturePlot(
  all_combined_data,
  features = Basal_Keratinocytes,
  reduction = "umap",
  cols = c("grey95", "#6495ED")  # grey (low) → green (high)
) &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
bkp[[1]] + bkp[[2]] + bkp[[3]]  + patchwork::plot_layout(ncol = 2)


# Suprabasal Keratinocytes
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Suprabasal_epidermal", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#00ff7f", "no" = "grey95")) +
  ggtitle("Suprabasal_epidermal") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )


Suprabasal_Keratinocytes <-c("DSC1", "DSG1", "KRT2", "GATA3")
# Generate FeaturePlots for each gene
dkp <- FeaturePlot(
  all_combined_data,
  features = Suprabasal_Keratinocytes,
  reduction = "umap",
  cols = c("grey95", "#00ff7f")  # grey (low) → green (high)
) 
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
dkp[[1]] + dkp[[2]] + dkp[[3]] + dkp[[4]] + patchwork::plot_layout(ncol = 2)

# Non-epidermal Keratinocytes
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Non_epidermal_keratinocytes", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#CC79A7", "no" = "grey95")) +
  ggtitle("Non_epidermal_keratinocytes") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )

Non_Epi_Keratinocytes <-c("FABP5", "SOSTDC1", "TNC", "PTN", "GJB2")
# Generate FeaturePlots for each gene
nekp <- FeaturePlot(
  all_combined_data,
  features = Non_Epi_Keratinocytes,
  reduction = "umap",
  cols = c("grey95", "#CC79A7")  # grey (low) → green (high)
) &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
nekp[[1]] + nekp[[2]] + nekp[[3]] + nekp[[4]] + nekp[[5]] + patchwork::plot_layout(ncol = 3)


# Melanocytes
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Melanocytes", "yes", "no")


DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#FF1493", "no" = "grey95")) +
  ggtitle("Melanocytes") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )


melanocytes <- c("MLANA", "PMEL", "TYR","TRPM1")
# Generate FeaturePlots for each gene
mp <- FeaturePlot(
  all_combined_data,
  features = melanocytes,
  reduction = "umap",
  cols = c("grey95", "#FF1493")  # grey (low) → green (high)
) &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
mp[[1]] + mp[[2]] + mp[[3]] + mp[[4]] + patchwork::plot_layout(ncol = 2)


# Sebocytes
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Sebocytes", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#BA55D3", "no" = "grey95")) +
  ggtitle("Sebocytes") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )

Sebocytes <- c("FADS2", "DHCR7", "INSIG1")
# Generate FeaturePlots for each gene
sp <- FeaturePlot(
  all_combined_data,
  features = Sebocytes,
  reduction = "umap",
  cols = c("grey95", "#BA55D3")  # grey (low) → green (high)
) &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
sp[[1]] + sp[[2]] + sp[[3]]  + patchwork::plot_layout(ncol = 2)

# Eccrine Secretory duct cells
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Eccrine_duct_cells", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#cd5c5c", "no" = "grey95")) +
  ggtitle("Eccrine_duct_cells") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )

Eccrine_duct_cells <- c("DCD", "KRT8", "KRT19")
# Generate FeaturePlots for each gene
egp <- FeaturePlot(
  all_combined_data,
  features = Eccrine_duct_cells,
  reduction = "umap",
  cols = c("grey95", "#cd5c5c")  # grey (low) → green (high)
) &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
egp[[1]] + egp[[2]] + egp[[3]] + patchwork::plot_layout(ncol = 2)


# Mural cells
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Mural", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#87CEEB", "no" = "grey95")) +
  ggtitle("Mural") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )


Mural <- c("ACTA2", "MYL9", "MYH11")
# Generate FeaturePlots for each gene
mup <- FeaturePlot(
  all_combined_data,
  features = Mural,
  reduction = "umap",
  cols = c("grey95", "#87CEEB")  # grey (low) → green (high)
) &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
mup[[1]] + mup[[2]] + mup[[3]]  + patchwork::plot_layout(ncol = 2)

# Endothelial_cells
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Endothelial_cells", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#9400D3", "no" = "grey95")) +
  ggtitle("Endothelial_cells") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )

Endothelial <- c("AQP1", "VWF", "SPARCL1", "CD93")
# Generate FeaturePlots for each gene
endop <- FeaturePlot(
  all_combined_data,
  features = Endothelial,
  reduction = "umap",
  cols = c("grey95", "#9400D3")  # grey (low) → green (high)
) &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
endop[[1]] + endop[[2]] + endop[[3]] + endop[[4]] + patchwork::plot_layout(ncol = 2)


# Fibroblasts 
all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Fibroblasts", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#228B22", "no" = "grey95")) +
  ggtitle("Fibroblasts") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )

Fibroblasts <- c("COL1A1", "FBN1", "FN1", "COL6A1")

# Generate FeaturePlots for each gene
fib <- FeaturePlot(
  all_combined_data,
  features = Fibroblasts,
  reduction = "umap",
  cols = c("grey95", "#228B22")  # grey (low) → green (high)
) &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.3),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 11)
  )
# Arrange the four plots in a 2×2 grid
fib[[1]] + fib[[2]] + fib[[3]] + fib[[4]] + patchwork::plot_layout(ncol = 2)

# ---------------------------
# Fibroblast cluster subset
# ---------------------------

fibro <- subset(all_combined_data, idents = "Fibroblasts")
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

# ---------------------------
# Clustering of fibroblast subset
# ---------------------------

DimPlot(
  fibro,
  group.by = "cluster_subtype",
  label = TRUE,
  cols = c(
    "Fb1" = "green",
    "Fb2" = "red"
  )
) +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank()
  )


# ---------------------------
# Immune Cluster Highlight
# ---------------------------

all_combined_data$highlight <- ifelse(all_combined_data$CellType == "Immune", "yes", "no")

DimPlot(
  all_combined_data, 
  reduction = "umap", 
  group.by = "highlight"
) +
  scale_color_manual(values = c("yes" = "#FFD700", "no" = "grey95")) +
  ggtitle("Immune") +
  theme_void() +                # remove axes, ticks, grids, and background
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)  # center title (optional)
  )


# ---------------------------
# Immune Subtyping
# ---------------------------
Immune_subset <- subset(all_combined_data, idents = "Immune")
DefaultAssay(Immune_subset) <- "Xenium"

# Preprocessing and UMAP
Immune_subset <- SCTransform(Immune_subset, assay = "Xenium", verbose = FALSE)

Immune_subset <- RunPCA(Immune_subset, npcs = num_pcs, verbose = FALSE)

Immune_subset <- FindNeighbors(Immune_subset, dims = 1:num_pcs)

Immune_subset <- FindClusters(Immune_subset, resolution = 2.0)

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

# Keep only genes present in the dataset
Immune_gene_sets_present <- lapply(Immune_gene_sets, function(x) x[x %in% rownames(Immune_subset)])

expr_mat <- GetAssayData(Immune_subset, layer = "data")

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

immune_colors <- c(
  Macrophage   = "#56B4E9",
  Mast_cell    = "#F0E442",
  DC1          = "#1B9E77",
  DC2          = "#984EA3",
  Langerhans   = "#E69F00",
  T_cell       = "#66A61E",
  NK_cell      = "#A65628",
  B_cell       = "#F781BF"
)

#----------------------------------
# Clustering of immune subset
#----------------------------------
DimPlot(Immune_subset, group.by = "subtype_annotation", label = TRUE, label.size = 5) +
  ggtitle("Immune Subtype Annotation (Per Cell)") +
  theme(legend.position = "bottom")

all_combined_data$Immune_Subtype_Cell <- NA
all_combined_data$Immune_Subtype_Cell[Cells(Immune)] <- Immune$subtype_annotation


#----------------------------------------------
# Merged Immune DotPlot
#----------------------------------------------

immune_markers <- c(
  # Myeloid
  "CD163", "C1QA", "C1QB",  "CD14", "C5AR1",         # Macrophage
  "CLEC9A", "WDFY4", "XCR1", "CPVL", "ITGAX",        # DC1
  "CLEC10A", "CD1C", "FCER1A",                         # DC2
  "CD207", "CD1A", #Langerhans
  "TPSAB1", "CPA3", "GATA2", "IL1RL1", # mast cell
  # Lymphoid
  "CD3D", "CD3E", "CD2", "CD52", "CTLA4",                  # T cells
  "NKG7", "GNLY", "NCAM1", "CTSW", "GZMK"                   # NK cells
)

immune_markers <- immune_markers[immune_markers %in% rownames(Immune_subset)]

Immune_plot <- subset(
  Immune_subset,
  subtype_annotation %in% c(
    "Macrophage",
    "DC1",
    "DC2",
    "Langerhans",
    "Mast_cell",
    "T_cell",
    "NK_cell"
  )
)

Immune_plot$subtype_annotation <- factor(
  Immune_plot$subtype_annotation,
  levels = c(
    "Macrophage",
    "DC1",
    "DC2",
    "Langerhans",
    "Mast_cell",
    "T_cell",
    "NK_cell"
  )
)

Immune_plot$subtype_annotation <- factor(
  Immune_plot$subtype_annotation,
  levels = rev(levels(Immune_plot$subtype_annotation))
)

guides(
  color = guide_colorbar(order = 2),
  size  = guide_legend(order = 1)
)


imm <- DotPlot(
  Immune_plot, 
  features = immune_markers, 
  group.by = "subtype_annotation",
  dot.min = 0.05,
) + scale_color_gradientn(
    colors = c("white", "lightblue", "darkblue"),
    name = "Average expression"
  ) +  RotatedAxis() +
  scale_size(range = c(0, 5), name = "Percentage expression") +
  theme(
    # Base text settings
    text = element_text(family = "Arial", size = 18),
    
    # Axis settings
    axis.text.x = element_text(family = "Arial", size = 18, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
    axis.text.y = element_text(family = "Arial", size = 18),
    axis.title = element_text(family = "Arial", size = 18),
    
    # Legend settings
    legend.text = element_text(family = "Arial", size = 10),
    legend.title = element_text(family = "Arial", size = 10),
    legend.key.size = unit(0.3, "cm"), # Smaller legend icons to match 5pt text
    
    # Border
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(color = "Average expression", x = "Markers", y = "Cell Type") +ggtitle("DotPlot: Immune Subtype Marker Expression")
imm
# 6. Save as TIFF,JPG or PNG 
ggsave(
  filename = "Immune_Markers_Dot_plot.tiff",
  plot = imm,
  device = "tiff",
  dpi = 600,
  compression = "lzw",
  width = 8.5,
  height = 6,
  units = "cm"
)

#------------------------------------
# Single Individual Analysis for Representative Images
#------------------------------------

d1922 <- subset(all_combined_data, subset = sample == "D1922")

DimPlot(d1922, reduction = "umap", group.by = "CellType", label = FALSE)

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
cell_annotations <- d1922@meta.data %>%
  mutate(cell_id = sub("^sample_D1922_", "", rownames(.))) %>%
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
p_all



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

fibro <- subset(d1922, idents = "Fibroblasts")
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

d1922$Fibroblast_Subtype_Cluster <- NA
d1922$Fibroblast_Subtype_Cluster[Cells(fibro)] <- fibro$cluster_subtype

#---------------------------------------
# Spatial locaiton of fibroblast subsets
#---------------------------------------

cell_annotations <- d1922@meta.data %>%
  mutate(cell_id = sub("^sample_D1922_", "", rownames(.))) %>%

  
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
