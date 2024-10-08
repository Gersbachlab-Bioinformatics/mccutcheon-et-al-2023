#!/usr/bin/env Rscript
options(future.globals.maxSize = 32000 * 1024^2)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))

# Set up argument parser
parser <- ArgumentParser(description = 'Process single cell RNA-seq data with Seurat')
parser$add_argument("-j", "--job", type = "integer", required = TRUE, help = "Integer used as index in the list of gRNAs")
parser$add_argument("--rds-files", nargs = '+', required = TRUE, help = "List of RDS files to be merged")
parser$add_argument("--output-basename", required = TRUE, help = "Basename used to create output file (args$output_basename + args$de_method + args$job)")
parser$add_argument("--min-umi-thres", type = "integer", default = 5, help = "Minimum UMI threshold (default: 5)")

args <- parser$parse_args()
job <- args$job
min_umi_thres <- args$min_umi_thres

# Load Seurat objects from the provided RDS files
seurat_objects <- lapply(args$rds_files, readRDS)
multimodal <- Reduce(function(x, y) merge(x, y), seurat_objects)
rm(seurat_objects)

# Set desired assay for clustering
DefaultAssay(multimodal) <- 'RNA'
multimodal <- NormalizeData(multimodal)

# Uncomment the following lines if you need to perform these steps
# multimodal <- FindVariableFeatures(multimodal, selection.method = 'vst')
# multimodal <- ScaleData(multimodal, features = all.genes)
# multimodal <- RunPCA(multimodal)
# multimodal <- FindNeighbors(multimodal, dims = 1:15)
# multimodal <- FindClusters(multimodal, resolution = 0.8)
# multimodal <- RunUMAP(multimodal, dims = 1:15)

DefaultAssay(multimodal) <- 'gRNAs'
gRNA_lib_size <- nrow(multimodal)
cell.indices <- vector("list", gRNA_lib_size)

for (i in 1:gRNA_lib_size) {
  cells <- multimodal@assays[["gRNAs"]]@counts[i, ]
  temp.df <- as.data.frame(which(cells > min_umi_thres))
  row.names(temp.df) <- NULL
  cell.indices[[i]] <- temp.df
}
names(cell.indices) <- row.names(multimodal)
cells_with_targeting_gRNA <- unlist(cell.indices[1:32])
cells_with_nt_gRNA <- unlist(cell.indices[33:40])
cells_with_only_nt <- setdiff(cells_with_nt_gRNA, cells_with_targeting_gRNA)
nt_masks <- grepl("^NT-", row.names(multimodal[['gRNAs']]))

# Assign gRNA index from command line argument
grna_ix <- as.integer(args$job)

# Perform differential analysis
tryCatch({
  Idents(object = multimodal) <- 0
  Idents(object = multimodal, cells = which(as.vector(multimodal[['gRNAs']][grna_ix, ] >= min_umi_thres))) <- 1
  Idents(object = multimodal, cells = cells_with_only_nt) <- 2
  protein.de.markers <- FindMarkers(
    multimodal, 
    ident.1 = 1, 
    ident.2 = 2, 
    test.use = 'MAST', 
    assay = 'RNA', 
    features = all.genes,
    logfc.threshold = 0,
    min.pct = 0,
    min.diff.pct = -Inf,
    min.cells.feature = 0,
    min.cells.group = 0,
    pseudocount.use = 0.001,
    verbose = FALSE
  )
  protein.de.markers['grna'] <- row.names(multimodal[['gRNAs']])[grna_ix]
  protein_de_tests <- protein.de.markers
}, error = function(e) {
  NA
})

protein_de_tests <- list(protein_de_tests)
protein_de_tests_all <- data.table::rbindlist(protein_de_tests[!is.na(protein_de_tests)])
protein_de_tests_all[, target_gene := unlist(lapply(protein_de_tests[!is.na(protein_de_tests)], rownames))]

write.table(protein_de_tests_all, 
            paste(args$output_basename, args$de_method, 
                  paste0("umi_thres_", min_umi_thres), 
                  args$job, "txt", sep = "."), 
            sep = '\t', quote = FALSE, row.names = FALSE)
