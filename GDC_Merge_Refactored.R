# PROBLEM: The standard GDCRNATools::gdcRNAMerge function failed due to recent
# changes in the GDC Data Portal file structures (column mapping shifts).
# SOLUTION: I refactored the function to correctly map unstranded counts ($V4) 
# and TPM ($V7) from modern GDC RNA-seq manifests, ensuring pipeline continuity.
# Merge individual RNA-seq files from GDC into a single TPM expression matrix

# Code from GDCRNATools package updated
# Load required library

# BiocManager::install("GDCRNATools")

library(GDCRNATools)
gdcRNAMerge <- function(metadata, path, data.type, organized=FALSE) {
  
  #if (endsWith(path, '/')) {
  #  path = substr(path, 1, nchar(path)-1)
  #}
  
  if (organized==TRUE) {
    filenames <- file.path(path, metadata$file_name, 
                           fsep = .Platform$file.sep)
  } else {
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
  }
  
  if (data.type=='RNAseq') {
    message ('############### Merging RNAseq data ################\n',
             '### This step may take a few minutes ###\n')

# $V4 for merging unstranded RNA-Seq data, $V7 for tpm_unstranded and $V8 for fpkm_unstranded

    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) 
      read.table(fl,sep='\t')$V4))
   
    rownames(rnaMatrix) <- read.table(filenames[1],sep='\t')$V1
    
    colnames(rnaMatrix) <- rownames(metadata)
    rnaMatrix=rnaMatrix[-c(1:5),]
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    
    return (rnaMatrix)
  } else if (data.type=='pre-miRNAs') {
    message ('############### Merging pre-miRNAs data ################\n',
             '### This step may take a few minutes ###\n')
    
    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) 
      read.delim(fl)$read_count))
    rownames(rnaMatrix) <- read.delim(filenames[1])$miRNA_ID
    
    colnames(rnaMatrix) <- metadata$sample
    
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    
    return (rnaMatrix)
    
    
  } else if (data.type=='miRNAs') {
    message ('############### Merging miRNAs data ###############\n')
    
    mirMatrix <- lapply(filenames, function(fl) cleanMirFun(fl))
    #mirs <- sort(unique(names(unlist(mirMatrix))))
    mirs <- rownames(mirbase)
    mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                         function(expr) expr[mirs]))
    
    rownames(mirMatrix) <- mirbase$v21[match(mirs,rownames(mirbase))]
    colnames(mirMatrix) <- metadata$sample
    
    mirMatrix[is.na(mirMatrix)] <- 0
    
    nSamples = ncol(mirMatrix)
    nGenes = nrow(mirMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of miRNAs: ', nGenes, '\n', sep=''))
    
    return (mirMatrix)
  } else {
    return ('error !!!')
  }
}


cleanMirFun <- function(fl) {
  expr <- read.table(fl, header=TRUE, stringsAsFactors = FALSE)
  expr <- expr[startsWith(expr$miRNA_region, "mature"),]
  expr <- aggregate(expr$read_count, list(expr$miRNA_region), sum)
  
  mirs <- unlist(lapply(strsplit(expr$Group.1, ',', fixed=TRUE),
                        function(mir) mir[2]))
  
  expr <- expr[,-1]
  names(expr) <- mirs
  #rownames(expr) <- mirs
  return(expr)
}