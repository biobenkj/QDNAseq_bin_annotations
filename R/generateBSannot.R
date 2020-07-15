##generate BS bins in QDNAseq annotation(ish) form
library(GenomicRanges)
library(Biostrings)
library(BiocParallel)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)

#set bin sizes
binSizes <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000)

#generate bin annotations for hg38
binAnnots.hg38 <- mclapply(binSizes, function(x) {
  bins <- createBSBins(BSgenome.Hsapiens.UCSC.hg38, bin.size = x, blacklist = "ENCFF356LFX.bed")
  return(bins)
}, mc.cores = ifelse(length(binSizes) <= 16, length(binSizes), 16))

names(binAnnots.hg38) <- paste0("hg38_", binSizes, "kbp")

saveRDS(binAnnots.hg38, file = "QDNAseq_hg38_all_BS_bin_annots_precomputed.rds")

#generate bin annotations for hg19
library(BSgenome.Hsapiens.UCSC.hg19)
binAnnots.hg19 <- mclapply(binSizes, function(x) {
  bins <- createBSBins(BSgenome.Hsapiens.UCSC.hg19, bin.size = x)
  return(bins)
}, mc.cores = ifelse(length(binSizes) <= 16, length(binSizes), 16))

names(binAnnots.hg19) <- paste0("hg19_", binSizes, "kbp")

saveRDS(binAnnots.hg19,
        file = "QDNAseq_hg19_all_BS_bin_annots_precomputed.rds")

