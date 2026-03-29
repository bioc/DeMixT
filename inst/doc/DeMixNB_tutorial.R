## ----example1-1, eval=FALSE---------------------------------------------------
# # --------------------------------------
# # Step 0: Install and load packages
# # --------------------------------------
# # Install packages
# install.packages('Seurat' ,'hdf5r', 'Rfast2')
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("glmGamPoi")
# 
# # Load library
# library(Seurat)
# library(DeMixT)
# library(dplyr)
# library(ggplot2)
# 
# # --------------------------------------
# # Step 1: Load and Preprocess Data
# # --------------------------------------
# # Load the LUSC spatial transcriptomics dataset (Seurat object)
# 
# data_dir <- "~/lusc_public_data"
# 
# spatial_img <- Read10X_Image(
#   image.dir = file.path(data_dir, "spatial"),
#   image.name = "tissue_lowres_image.png",
#   filter.matrix = TRUE
# )
# lusc_seurat <- Load10X_Spatial(data.dir = file.path(data_dir,
#                                                     "filtered_feature_bc_matrix/"),
#                                filename = "filtered_feature_bc_matrix.h5",
#                                assay = "Spatial",slice = "Lung_P14",
#                                filter.matrix = TRUE,
#                                to.upper = FALSE,
#                                image = spatial_img)
# 
# # Display initial dataset dimensions
# message("Before quality control, the dataset contains ",
#         ncol(lusc_seurat), " spots and ",
#         nrow(lusc_seurat), " genes.")

## ----example1-2, eval=FALSE---------------------------------------------------
# # --------------------------------------
# # Step 2: Filter Genes
# # --------------------------------------
# # Remove unwanted genes (e.g., mitochondrial, pseudogenes, etc.)
# # Regular expressions identify genes to exclude based on their names
# unwanted_genes <- grepl("^(MT-|MIR-|LINC|A[A-Z]|B[A-Z]|
#                         C[A-Z]|Rik|Gm|Vmn|Olfr|Trav)",
#                         rownames(lusc_seurat))
# lusc_seurat <- lusc_seurat[!unwanted_genes, ]
# 
# # Filter genes with low expression (detected in fewer than 10 spots)
# count_data <- GetAssayData(lusc_seurat, assay = "Spatial",layer = "counts")
# gene_detection_counts <- rowSums(count_data > 0)
# genes_to_keep <- names(gene_detection_counts[gene_detection_counts >= 10])
# 
# # --------------------------------------
# # Step 3: Load and Process Cell Metadata
# # --------------------------------------
# # Load cell metadata from CSV file
# cell_metadata <- read.csv(file.path(data_dir,
#                                     "LUSC_spot_annotation.csv"))
# 
# # Filter valid barcodes and classify cells as tumor or non-tumor
# cell_metadata <- cell_metadata %>%
#   filter(!is.na(Barcode) & Barcode != "" & Barcode %in% colnames(count_data)) %>%
#   mutate(class = ifelse(class2 == "t", "tumor", "non_tumor"))
# 
# # Summarize spot-level cell information (e.g., coordinates, cell type proportions)
# barcode_summary <- cell_metadata %>%
#   group_by(Barcode) %>%
#   summarise(
#     x_coor = mean(x_fullres),
#     y_coor = mean(y_fullres),
#     tumor_proportion = sum(class == "tumor") / n(),
#     stromal_proportion = sum(class2 == "f") / n(),
#     immune_proportion = sum(class2 == "l") / n(),
#     non_tumor_proportion = 1 - tumor_proportion,
#     number_of_cells = n()
#   )
# 
# # --------------------------------------
# # Step 4: Quality Control (QC) on Spots
# # --------------------------------------
# # Subset Seurat object to keep spots with sufficient counts and features
# lusc_seurat <- subset(
#   x = lusc_seurat,
#   subset = (nCount_Spatial >= 600) & (nFeature_Spatial >= 300),
#   features = genes_to_keep
# )
# 
# # Update barcode summary to match filtered Seurat object
# barcode_summary <- cell_metadata %>%
#   filter(Barcode %in% colnames(lusc_seurat)) %>%
#   group_by(Barcode) %>%
#   summarise(
#     x_coor = mean(x_fullres),
#     y_coor = mean(y_fullres),
#     tumor_proportion = sum(class == "tumor") / n(),
#     stromal_proportion = sum(class2 == "f") / n(),
#     immune_proportion = sum(class2 == "l") / n(),
#     non_tumor_proportion = 1 - tumor_proportion,
#     tumor_cell = sum(class == "tumor"),
#     number_of_cells = n()
#   )
# 
# # Filter spots with at least 5 cells and tumor proportion between 0.05 and 0.95
# barcode_summary_filtered <- barcode_summary %>%
#   filter(number_of_cells >= 5 & (tumor_proportion == 0 |
#                                    (tumor_proportion >= 0.05
#                                     & tumor_proportion <= 0.95)))
# 
# # Display post-QC dataset dimensions
# count_matrix <- GetAssayData(lusc_seurat, assay = "Spatial", layer = "counts")
# message("After quality control, the dataset contains ",
#         nrow(barcode_summary_filtered), " spots and ",
#         nrow(count_matrix), " genes.")

## ----example1-3, eval=FALSE---------------------------------------------------
# # --------------------------------------
# # Step 5: Classify Spots
# # --------------------------------------
# # Identify mixed (tumor + non-tumor) and reference (non-tumor) spots
# barcode_mixed <- barcode_summary_filtered %>%
#   filter(tumor_proportion != 0) %>%
#   pull(Barcode)
# barcode_reference <- barcode_summary_filtered %>%
#   filter(tumor_proportion == 0) %>%
#   pull(Barcode)
# 
# # Create data frames for mixed and reference spots
# cell_proportions <- barcode_summary_filtered %>%
#   select(Barcode, tumor_proportion, stromal_proportion,
#          immune_proportion, non_tumor_proportion, number_of_cells)
# 
# mixed_spots <- cell_proportions %>%
#   filter(Barcode %in% barcode_mixed) %>%
#   mutate(Type = "Candidate mixed spot")
# reference_spots <- cell_proportions %>%
#   filter(Barcode %in% barcode_reference) %>%
#   mutate(Type = "Candidate reference spot")
# 
# # Combine spot data
# combined_spots <- rbind(mixed_spots, reference_spots)
# 
# # --------------------------------------
# # Step 6: Subset Seurat Object
# # --------------------------------------
# # Create Seurat objects for reference and mixed spots
# reference_seurat <- lusc_seurat[, colnames(lusc_seurat) %in% barcode_reference]
# mixed_seurat <- lusc_seurat[, colnames(lusc_seurat) %in% barcode_mixed]
# combined_seurat <- lusc_seurat[, colnames(lusc_seurat) %in% combined_spots$Barcode]
# 
# # Add sample type metadata
# combined_seurat$sample_type <- ifelse(
#   colnames(combined_seurat) %in% barcode_reference, "Normal",
#   ifelse(colnames(combined_seurat) %in% barcode_mixed, "Tumor Mixed", NA)
# )
# 
# # --------------------------------------
# # Step 7: Visualize Spot Classifications
# # --------------------------------------
# # Plot spatial distribution of normal and tumor-mixed spots
# spatial_plot <- SpatialDimPlot(
#   combined_seurat,
#   group.by = "sample_type",
#   cols = c("Normal" = "blue", "Tumor Mixed" = "red"),
#   pt.size.factor = 1.5,
#   interactive = FALSE,
#   crop = FALSE
# )
# 
# # --------------------------------------
# # Step 8: Prepare Data for DeMixNB
# # --------------------------------------
# # Extract count matrices for tumor and normal spots
# genes_shared <- intersect(rownames(reference_seurat), rownames(mixed_seurat))
# tumor_counts <- GetAssayData(mixed_seurat, layer = "counts")[genes_shared, ]
# normal_counts <- GetAssayData(reference_seurat, layer = "counts")[genes_shared, ]
# 
# # --------------------------------------
# # Step 9: SCTransform and Dimensionality Reduction
# # --------------------------------------
# # Normalize and transform data using SCTransform
# combined_seurat <- SCTransform(combined_seurat, assay = "Spatial", verbose = FALSE)
# 
# # Run PCA, clustering, and UMAP
# combined_seurat <- RunPCA(combined_seurat, assay = "SCT", verbose = FALSE)
# combined_seurat <- FindNeighbors(combined_seurat, reduction = "pca", dims = 1:50)
# combined_seurat <- FindClusters(combined_seurat, verbose = FALSE, resolution = 0.8)
# combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:50)
# 
# # --------------------------------------
# # Step 10: Identify Spatially Variable Genes
# # --------------------------------------
# # Select top 1000 spatially variable genes using Moran's I
# # This step will take a while to process
# # User can choose number of genes as needed
# combined_seurat <- FindSpatiallyVariableFeatures(
#   combined_seurat,
#   assay = "SCT",
#   features = VariableFeatures(combined_seurat)[1:1000],
#   selection.method = "moransi"
# )
# 
# # Extract Moran's I results and sort by observed Moran's I
# moransi_results <- combined_seurat@assays$SCT@meta.features %>%
#   na.exclude() %>%
#   arrange(desc(MoransI_observed))
# top_genes <- rownames(moransi_results)[1:1000]
# 
# 
# # --------------------------------------
# # Step 11: Run DeMixNB
# # --------------------------------------
# # Define gene list sizes and prepare gene lists
# # For illustration purpose, this tutorial only uses gene size of 1000.
# 
# gene_list_sizes <- c(1000)
# 
# # Initialize deconvolution results matrix
# deconvolution_results <- matrix(nrow = length(barcode_mixed), ncol = 0)
# 
# 
# for (i in 1:length(gene_list_sizes)) {
#   gl <- top_genes[1:gene_list_sizes[i]]
#   filtered_mat <- count_matrix[gl, ]  # Subset count matrix for selected genes
#   data.N1 <- as.matrix(filtered_mat[, barcode_reference])  # Normal spot counts
#   data.Y <- as.matrix(filtered_mat[, barcode_mixed])
#   matrix_data <- as.matrix(filtered_mat[, barcode_mixed])
#   data.Y <- SummarizedExperiment(assays = list(counts = matrix_data))
#   matrix_data <- as.matrix(filtered_mat[, barcode_reference])
#   data.N1 <- SummarizedExperiment(assays = list(counts = matrix_data))
#   test_simulation <- DeMixNB(data.Y = data.Y, data.N1 = data.N1, niter = 30, tol = 1e-3)
# 
#   ngene.selected <- gene_list_sizes[i]
#   name <- paste0("DemixNB_LUSC_", ngene.selected, ".RData")
#   save(test_simulation, file = paste0('~/Deconvolution_result/', name))
# 
#   # Get pi results and add gene count prefix to column names
# 
#   pi_results <- test_simulation$pi_t_summary
#   colnames(pi_results) <- paste0(ngene.selected, "_", colnames(pi_results))
#   deconvolution_results <- cbind(deconvolution_results, pi_results)
# }
# 
# 
# # --------------------------------------
# # Step 12: Visualize DeMixNB Results
# # --------------------------------------
# # Subset Seurat object for mixed spots and add tumor proportion estimates
# consistent_rows_1000 <- deconvolution_results$`1000_consistency` == "consistent"
# 
# # Select the 500 gene results specifically
# tumor_proportions <- data.frame(
#   sample.id = barcode_mixed[consistent_rows_1000],
#   PiT = deconvolution_results$`1000_pi_run1`[consistent_rows_1000]
# )
# 
# mixed_seurat <- subset(lusc_seurat, cells = tumor_proportions$sample.id)
# matched_indices <- match(Cells(mixed_seurat), tumor_proportions$sample.id)
# tumor_proportions <- tumor_proportions[matched_indices, ]
# mixed_seurat$PiT <- tumor_proportions$PiT
# 
# # Plot tumor proportion estimates spatially
# tumor_proportion_plot <- SpatialFeaturePlot(
#   mixed_seurat,
#   features = "PiT",
#   pt.size.factor = 1.5,
#   crop = FALSE
# ) + theme(legend.position = "top")
# 
# # Save or display plots
# ggsave("~/spatial_plot.png", spatial_plot)
# ggsave("~/tumor_proportion_plot.png", tumor_proportion_plot)
# 

## ----plot1, echo=FALSE, fig.align='center'------------------------------------
knitr::include_graphics(path = paste0("spatial_plot.png"))

## ----plot2, echo=FALSE, fig.align='center'------------------------------------
knitr::include_graphics(path = paste0("tumor_proportion_plot.png"))

## ----session-info-------------------------------------------------------------
sessionInfo()

