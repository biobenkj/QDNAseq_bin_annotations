##generate hg38 bin annotations for QDNAseq
library(QDNAseq)
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)

#set bin sizes
binSizes <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000)

#generate bin annotations
binAnnots.hg38 <- mclapply(binSizes, function(x) {
  bins <- QDNAseq::createBins(BSgenome.Hsapiens.UCSC.hg38,
                              binSize = x)
  return(bins)
}, mc.cores = ifelse(length(binSizes) <= 16, length(binSizes), 16))

names(binAnnots.hg38) <- paste0("hg38_", binSizes, "kbp")

saveRDS(binAnnots.hg38,
        file = "QDNAseq_hg38_all_bin_annots_precomputed.rds")

##generate mm9 bin annotations for QDNAseq
library(QDNAseq)
library(Biobase)
library(BSgenome.Mmusculus.UCSC.mm9)
library(parallel)

#set bin sizes
binSizes <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000)

#generate bin annotations
binAnnots.mm9 <- mclapply(binSizes, function(x) {
  bins <- QDNAseq::createBins(BSgenome.Mmusculus.UCSC.mm9,
                              binSize = x)
  return(bins)
}, mc.cores = ifelse(length(binSizes) <= 16, length(binSizes), 16))

names(binAnnots.mm9) <- paste0("mm9_", binSizes, "kbp")

saveRDS(binAnnots.mm9,
        file = "QDNAseq_mm9_all_bin_annots_precomputed.rds")

##generate mm10 bin annotations for QDNAseq
library(QDNAseq)
library(Biobase)
library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)

#set bin sizes
binSizes <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000)

#generate bin annotations
binAnnots.mm10 <- mclapply(binSizes, function(x) {
  bins <- QDNAseq::createBins(BSgenome.Mmusculus.UCSC.mm10,
                              binSize = x)
  return(bins)
}, mc.cores = ifelse(length(binSizes) <= 16, length(binSizes), 16))

names(binAnnots.mm10) <- paste0("mm10_", binSizes, "kbp")

saveRDS(binAnnots.mm10,
        file = "QDNAseq_mm10_all_bin_annots_precomputed.rds")

