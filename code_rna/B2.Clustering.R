# ml R/4.2.0
# R

################################################
### BiNGS Bulk RNA - Pipeline - Transcriptomics
################################################

### Load Libraries ###
.libPaths(c("/hpc/packages/minerva-rocky9/rpackages/4.2.0/site-library", "/hpc/packages/minerva-rocky9/rpackages/bioconductor/3.15"))
library(DESeq2)
library(fst)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(rtracklayer)

### Detect script location ###
script_path <- dirname(normalizePath(sys.frame(1)$ofile))
base_dir <- dirname(script_path)

### Define Directories ###
raw_dir <- file.path(base_dir, "data_rna/raw")
preprocessed_dir <- file.path(base_dir, "data_rna/preprocessed")
processed_dir <- file.path(base_dir, "data_rna/processed")
compiled_data_dir <- file.path(processed_dir, "compiled_data")
figures_dir <- file.path(processed_dir, "figures")
analysis_dir <- file.path(processed_dir, "analysis")
anno_dir <- file.path(base_dir, "supporting_files/annotation")

### Create Required Directories ###
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "clustering"), recursive = TRUE, showWarnings = FALSE)

#################################################
### Detect organism automatically from folders
#################################################

detected_organisms = list.dirs(anno_dir, recursive = FALSE, full.names = FALSE)

if (length(detected_organisms) == 1) {
    organism <- detected_organisms[1]
    message(paste0("✅ Detected organism: ", organism))
} else if (length(detected_organisms) == 0) {
    stop("❌ No organism folders detected under annotation/. Please prepare annotation files first.")
} else {
    message("⚠️ Multiple organisms detected:")
    for (org in detected_organisms) message("- ", org)
    
    repeat {
        cat("\nPlease enter the organism you want to use exactly as listed above: ")
        organism <- tolower(trimws(readLines(con = stdin(), n = 1)))
        if (organism %in% detected_organisms) {
            message(paste0("✅ Using organism: ", organism))
            break
        } else {
            message("❌ Invalid input. Try again.")
        }
    }
}

#################################################
### Load annotation GTF
#################################################

if (organism == "homo_sapiens") {
    subfolder = "grch38_gencode_36"
    gtf_file <- file.path(anno_dir, organism, subfolder, "data", "gencode.v36.annotation.gtf")
    gtf_url <- "https://ulukag01.dmz.hpc.mssm.edu/supporting_files/annotation/homo_sapiens/grch38_gencode_36/data/gencode.v36.annotation.gtf"
} else if (organism == "mus_musculus") {
    subfolder = "grcm38_gencode_M25"
    gtf_file <- file.path(anno_dir, organism, subfolder, "data", "gencode.vM25.annotation.gtf")
    gtf_url <- "https://ulukag01.dmz.hpc.mssm.edu/supporting_files/annotation/mus_musculus/grcm38_gencode_M25/data/gencode.vM25.annotation.gtf"
} else {
    cat("\n⚠️ You picked a custom organism.\nPlease provide full path to a GTF file (*.annotation.gtf): ")
    gtf_file <- trimws(readLines(con = stdin(), n = 1))
    gtf_url <- NULL
}

if (!file.exists(gtf_file) && !is.null(gtf_url)) {
    message("⬇️ Downloading GTF file from: ", gtf_url)
    dir.create(dirname(gtf_file), recursive = TRUE, showWarnings = FALSE)
    download.file(url = gtf_url, destfile = gtf_file, mode = "wb")
}

gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data) %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  distinct(gene_id, .keep_all = TRUE)

#################################################
### Load compiled data
#################################################

sample_metadata <- fst::read_fst(file.path(compiled_data_dir, "dtl_sample_metadata.fst"))
dtl_salmon_gene_counts <- fst::read_fst(file.path(compiled_data_dir, "dtl_salmon_gene_counts.fst"))
dtl_salmon_gene_counts_medianratios <- fst::read_fst(file.path(compiled_data_dir, "dtl_salmon_gene_counts_medianratios.fst"))

# Prepare Counts Matrix
if ("gene_id" %in% colnames(dtl_salmon_gene_counts_medianratios)) {
  rownames(dtl_salmon_gene_counts_medianratios) <- dtl_salmon_gene_counts_medianratios$gene_id
  dtl_salmon_gene_counts_medianratios$gene_id <- NULL
}

dtl_counts <- as.matrix(round(dtl_salmon_gene_counts_medianratios))
mode(dtl_counts) <- "integer"

# Prepare Sample Metadata
coldata <- sample_metadata[, c("condition", "sequencing_batch")]
rownames(coldata) <- sample_metadata$sample_id

if (!all(rownames(coldata) == colnames(dtl_counts))) {
    stop("❌ Error: rownames(coldata) must match colnames(counts)!")
} else {
    message("🎯 Sample metadata and counts matrix match.")
}

#################################################
### Create DESeq2 Object and PCA
#################################################

dds <- DESeqDataSetFromMatrix(countData = dtl_counts, colData = coldata, design = ~ condition)
vsd <- vst(dds, blind = TRUE)

# PCA colored by condition
pcaData <- plotPCA(vsd, ntop = 1000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_condition <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = 2, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_light()

pdf(file.path(figures_dir, "clustering", "preqc_pca_colored_by_condition.pdf"), width = 8, height = 8)
print(pca_condition)
dev.off()
message("✅ PCA (condition) saved at: ", file.path(figures_dir, "clustering", "preqc_pca_colored_by_condition.pdf"))

# PCA colored by sequencing batch
if ("sequencing_batch" %in% colnames(coldata) && length(unique(coldata$sequencing_batch)) > 1) {
  pca_batch <- ggplot(pcaData, aes(x = PC1, y = PC2, color = sequencing_batch, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = 2, size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_light()
  
  pdf(file.path(figures_dir, "clustering", "preqc_pca_colored_by_batch.pdf"), width = 8, height = 8)
  print(pca_batch)
  dev.off()
  message("✅ PCA (batch) saved at: ", file.path(figures_dir, "clustering", "preqc_pca_colored_by_batch.pdf"))
} else {
  message("👉 Skipping PCA colored by batch: only one unique batch detected.")
}

# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(figures_dir, "clustering", "preqc_distance_matrix.pdf"), width = 8, height = 8)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
dev.off()
message("✅ Distance matrix heatmap saved at: ", file.path(figures_dir, "clustering", "preqc_distance_matrix.pdf"))

### Ask user if they want to remove low quality samples and redo PCA ###
cat("\nWould you like to remove low-quality samples and re-make PCA/Distance matrix plots? (yes/no): ")
remove_samples_answer <- tolower(trimws(readLines(con = stdin(), n = 1)))

if (remove_samples_answer == "yes") {
  
  cat("\nHere are the available sample IDs:\n")
  print(colnames(vsd))
  
  cat("\nPlease enter the sample IDs to remove (comma-separated like IFE_1,HF_3): ")
  low_quality_samples_input <- trimws(readLines(con = stdin(), n = 1))
  low_quality_samples <- unlist(strsplit(low_quality_samples_input, ","))
  low_quality_samples <- trimws(low_quality_samples)
  
  # Filter vsd and coldata
  vsd_filtered <- vsd[, !(colnames(vsd) %in% low_quality_samples)]
  coldata_filtered <- coldata[!(rownames(coldata) %in% low_quality_samples), ]
  
  ### PCA Plots (Post-QC) ###
  pcaData_postqc <- plotPCA(vsd_filtered, ntop = 1000, returnData = TRUE)
  percentVar_postqc <- round(100 * attr(pcaData_postqc, "percentVar"))
  
  # PCA colored by condition (Post-QC)
  pca_condition_postqc <- ggplot(pcaData_postqc, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = 2, size = 3) +
    xlab(paste0("PC1: ", percentVar_postqc[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_postqc[2], "% variance")) +
    theme_light()
  
  pdf(file.path(figures_dir, "clustering", "postqc_pca_colored_by_condition.pdf"), width = 8, height = 8)
  print(pca_condition_postqc)
  dev.off()
  message("✅ Post-QC PCA plot (condition) saved:")
  print(file.path(figures_dir, "clustering", "postqc_pca_colored_by_condition.pdf"))
  
  # PCA colored by batch (Post-QC, if applicable)
  if ("sequencing_batch" %in% colnames(coldata_filtered) && length(unique(coldata_filtered$sequencing_batch)) > 1) {
    pca_batch_postqc <- ggplot(pcaData_postqc, aes(x = PC1, y = PC2, color = sequencing_batch, label = name)) +
      geom_point(size = 3) +
      geom_text(vjust = 2, size = 3) +
      xlab(paste0("PC1: ", percentVar_postqc[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar_postqc[2], "% variance")) +
      theme_light()
    
    pdf(file.path(figures_dir, "clustering", "postqc_pca_colored_by_batch.pdf"), width = 8, height = 8)
    print(pca_batch_postqc)
    dev.off()
    message("✅ Post-QC PCA plot (batch) saved:")
    print(file.path(figures_dir, "clustering", "postqc_pca_colored_by_batch.pdf"))
  } 
  
  ### Sample Distance Matrix (Post-QC) ###
  sampleDists_postqc <- dist(t(assay(vsd_filtered)))
  sampleDistMatrix_postqc <- as.matrix(sampleDists_postqc)
  
  rownames(sampleDistMatrix_postqc) <- colnames(vsd_filtered)
  colnames(sampleDistMatrix_postqc) <- colnames(vsd_filtered)
  
  pdf(file.path(figures_dir, "clustering", "postqc_distance_matrix.pdf"), width = 8, height = 8)
  pheatmap(sampleDistMatrix_postqc,
           clustering_distance_rows = sampleDists_postqc,
           clustering_distance_cols = sampleDists_postqc,
           col = colors,
           display_numbers = FALSE,
           number_color = "black",
           fontsize_number = 6)
  dev.off()
  message("✅ Post-QC sample distance matrix saved:")
  print(file.path(figures_dir, "clustering", "postqc_distance_matrix.pdf"))
  
  # Replace for further steps
  vsd <- vsd_filtered
  coldata <- coldata_filtered
  
} 

### Ask user if they want batch correction ###
if ("sequencing_batch" %in% colnames(coldata) && length(unique(coldata$sequencing_batch)) > 1) {
  
  cat("\nWould you like to apply batch correction before PCA? (yes/no): ")
  batch_correct_answer <- tolower(trimws(readLines(con = stdin(), n = 1)))
  
  if (batch_correct_answer == "yes") {
    
    ### Perform Batch Correction using limma::removeBatchEffect ###
    library(limma)
    
    message("Performing Batch Correction using limma::removeBatchEffect")
    vsd_corrected_matrix <- removeBatchEffect(assay(vsd),
                                              batch = coldata$sequencing_batch)
    
    # Save batch corrected counts
    batch_corrected_fst <- file.path(compiled_data_dir, "dtl_salmon_gene_counts_vsd_batchcorrected.fst")
    batch_corrected_csv <- file.path(compiled_data_dir, "dtl_salmon_gene_counts_vsd_batchcorrected.csv")
    
    fst::write_fst(as.data.frame(vsd_corrected_matrix), batch_corrected_fst)
    write.csv(as.data.frame(vsd_corrected_matrix), batch_corrected_csv, row.names = TRUE)
    
    message("✅ Batch corrected counts saved as FST:")
    print(batch_corrected_fst)
    message("✅ Batch corrected counts saved as CSV:")
    print(batch_corrected_csv)
    
    ### PCA Plot after Batch Correction ###
    vsd_corrected_for_pca <- vsd
    assay(vsd_corrected_for_pca) <- vsd_corrected_matrix
    
    pcaData_batchcorrected <- plotPCA(vsd_corrected_for_pca, ntop = 1000, returnData = TRUE)
    percentVar_batchcorrected <- round(100 * attr(pcaData_batchcorrected, "percentVar"))
    
    pca_batchcorrected_plot <- ggplot(pcaData_batchcorrected, aes(x = PC1, y = PC2, color = condition, label = name)) +
      geom_point(size = 3) +
      geom_text(vjust = 2, size = 3) +
      xlab(paste0("PC1: ", percentVar_batchcorrected[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar_batchcorrected[2], "% variance")) +
      theme_light()
    
    pdf(file.path(figures_dir, "clustering", "postqc_pca_batchcorrected_colored_by_condition.pdf"), width = 8, height = 8)
    print(pca_batchcorrected_plot)
    dev.off()
    
    message("✅ Batch-corrected PCA plot saved:")
    print(file.path(figures_dir, "clustering", "postqc_pca_batchcorrected_colored_by_condition.pdf"))
    
    pca_batchcorrected_plot_batch <- ggplot(pcaData_batchcorrected, aes(x = PC1, y = PC2, color = sequencing_batch, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = 2, size = 3) +
    xlab(paste0("PC1: ", percentVar_batchcorrected[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_batchcorrected[2], "% variance")) +
    theme_light()
  
    pdf(file.path(figures_dir, "clustering", "postqc_pca_batchcorrected_colored_by_batch.pdf"), width = 8, height = 8)
    print(pca_batchcorrected_plot_batch)
    dev.off()
    message("✅ Batch-corrected PCA plot (batch) saved:")
    print(file.path(figures_dir, "clustering", "postqc_pca_batchcorrected_colored_by_batch.pdf"))
    cat("\n⚠️ This plot is just to see if you have batch effect that needs to be corrected. The C and D scripts are not capable of running batch corrected analysis.")

  } else {
    message("Skipping batch corrected PCA for now as per user request.")
  }
} 