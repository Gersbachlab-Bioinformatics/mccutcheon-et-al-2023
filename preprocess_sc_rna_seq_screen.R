# Load necessary libraries
library(Seurat)
library(argparse)
library(DoubletFinder)

# Set up argument parser
parser <- ArgumentParser(description = 'Process single cell RNA-seq data with Seurat')
parser$add_argument('--input_dir', required = TRUE, help = 'Directory containing the input files')
parser$add_argument('--output_dir', required = TRUE, help = 'Directory to save the output files')
parser$add_argument('--donor', required = TRUE, help = 'Donor identifier')
parser$add_argument('--crispr_app', required = TRUE, help = 'CRISPR application type')

args <- parser$parse_args()

# Load data
RNA_guide_umis <- Read10X(args$input_dir)

# Split mRNA and guide counts into two separate Seurat Objects
RNA_umis <- RNA_guide_umis[["Gene Expression"]]
guide_umis <- RNA_guide_umis[["CRISPR Guide Capture"]]

# Create a Seurat Object with multimodal info: mRNA and gRNA
multimodal <- CreateSeuratObject(counts = RNA_umis)

# Add another independent assay ("guides") to store gRNA info
multimodal[["gRNAs"]] <- CreateAssayObject(counts = guide_umis)

# Add appropriate metadata to Seurat Object
multimodal$donor <- args$donor
multimodal$CRISPR_application <- args$crispr_app

# Define the percentage of mitochondrial and ribosomal genes for each cell
multimodal <- PercentageFeatureSet(multimodal, "^MT-", col.name = "percent_mito")
multimodal <- PercentageFeatureSet(multimodal, "^RP[SL]", col.name = "percent_ribo")

# QC filtering
multimodal <- subset(multimodal, subset = nFeature_RNA > 200 & percent_mito < 20 & percent_ribo > 5)

# Normalize mRNA counts
multimodal <- NormalizeData(multimodal)

# Run scaling, variable gene selection, PCA, and UMAP for visualization
multimodal <- FindVariableFeatures(multimodal, verbose = FALSE)
multimodal <- ScaleData(multimodal, vars.to.regress = c("nFeature_RNA", "percent_mito"), verbose = FALSE)
multimodal <- RunPCA(multimodal, verbose = FALSE, npcs = 20)
multimodal <- RunUMAP(multimodal, dims = 1:10, verbose = FALSE)

# Run DoubletFinder
nExp <- round(ncol(multimodal) * 0.08)  # expect 8% doublets
multimodal <- DoubletFinder::doubletFinder_v3(multimodal, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# Extract the correct column name for DF prediction
DF.name <- colnames(multimodal@meta.data)[grepl("DF.classification", colnames(multimodal@meta.data))]

# Remove all predicted doublets from our data
multimodal <- multimodal[, multimodal@meta.data[, DF.name] == "Singlet"]
print(dim(multimodal))

# Save Seurat object for subsequent merging with other donors
dir.create(args$output_dir, showWarnings = FALSE)
saveRDS(multimodal, file.path(args$output_dir, paste0("CRISPRi_", args$donor, "_qc.rds")))
