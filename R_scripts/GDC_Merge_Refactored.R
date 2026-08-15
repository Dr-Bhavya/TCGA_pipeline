library(GDCRNATools)

# ============================================================================
# FUNCTION: gdcRNAMerge
# Description: Aggregates individual sample quantification files from GDC into 
#              a single expression matrix based on the specified data type.
# Parameters:
#   - metadata: Dataframe containing GDC file manifests and sample mappings
#   - path: Directory path pointing to the raw GDC data downloads
#   - data.type: Target sequencing profile ('RNAseq', 'pre-miRNAs', 'miRNAs')
#   - organized: Boolean flag specifying if directories are structured or flat
# ============================================================================
gdcRNAMerge <- function(metadata, path, data.type, organized = FALSE) {
    if (organized == TRUE) {
        filenames <- file.path(path, metadata$file_name, fsep = .Platform$file.sep)
    } else {
        filenames <- file.path(path, metadata$file_id, metadata$file_name, fsep = .Platform$file.sep)
    }
    
    if (data.type == "RNAseq") {
        message("############### Merging RNAseq data ################\n### This step may take a few minutes ###\n")
        rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) read.table(fl, sep = "\t")$V4))
        rownames(rnaMatrix) <- read.table(filenames[1], sep = "\t")$V1
        colnames(rnaMatrix) <- rownames(metadata)
        rnaMatrix <- rnaMatrix[-c(1:5), ]
        message(paste("Number of samples: ", ncol(rnaMatrix), "\n", sep = ""), paste("Number of genes: ", nrow(rnaMatrix), "\n", sep = ""))
        return(rnaMatrix)
    } else if (data.type == "pre-miRNAs") {
        message("############### Merging pre-miRNAs data ################\n### This step may take a few minutes ###\n")
        rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) read.delim(fl)$read_count))
        rownames(rnaMatrix) <- read.delim(filenames[1])$miRNA_ID
        colnames(rnaMatrix) <- metadata$sample
        message(paste("Number of samples: ", ncol(rnaMatrix), "\n", sep = ""), paste("Number of genes: ", nrow(rnaMatrix), "\n", sep = ""))
        return(rnaMatrix)
    } else if (data.type == "miRNAs") {
        message("############### Merging miRNAs data ###############\n")
        mirMatrix <- lapply(filenames, function(fl) cleanMirFun(fl))
        mirs <- rownames(mirbase)
        mirMatrix <- do.call("cbind", lapply(mirMatrix, function(expr) expr[mirs]))
        rownames(mirMatrix) <- mirbase$v21[match(mirs, rownames(mirbase))]
        colnames(mirMatrix) <- metadata$sample
        mirMatrix[is.na(mirMatrix)] <- 0
        message(paste("Number of samples: ", ncol(mirMatrix), "\n", sep = ""), paste("Number of miRNAs: ", nrow(mirMatrix), "\n", sep = ""))
        return(mirMatrix)
    } else {
        return("error !!!")
    }
}

# ============================================================================
# FUNCTION: cleanMirFun
# Description: Sub-routine to parse microRNA expression profiles, filtering 
#              out untranslated precursors to isolate mature sequences.
# Parameters:
#   - fl: Target file path to parse
# ============================================================================
cleanMirFun <- function(fl) {
    expr <- read.table(fl, header = TRUE, stringsAsFactors = FALSE)
    expr <- expr[startsWith(expr$miRNA_region, "mature"), ]
    expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
    mirs <- unlist(lapply(strsplit(expr$Group.1, ",", fixed = TRUE), function(mir) mir[2]))
    expr <- expr[, -1]
    names(expr) <- mirs
    return(expr)
}