#' Compute effective GC content for bisulfite sequencing experiments
#' 
#' This is mainly used as a way to compute moderated GC content for
#' CNV sketching with QDNAseq from bisulfite sequencing data
#'
#' @param obj Input CpG BSseq object
#' @param bins A GRanges object of bin annotations generated from createBSBins
#' @param bin.size The size of the bin annotations (DEFAULT: 30kb)
#' @param nome Whether we are running in NOMe-seq mode
#' @param nome.obj The GpC methylation BSseq object
#'
#' @return A matrix of moderated effective GC-content for each sample
#' @export
#'
#' @examples
#' 

getEffectiveGC <- function(obj, bins, bin.size = 3e4,
                           nome = FALSE, nome.obj = NULL) {
  ## Simple wrapper function to compute the _effective_ gc content
  ## from bisulfite sequencing experiments
  
  # This expects a bsseq object as input
  if (!is(obj, "BSseq")) {
    stop("Input for this function must be a BSseq object.")
  }
  
  # Check the input if we are in NOMe-seq mode
  if (nome) {
    message("NOMe-seq mode.")
    if (is.null(nome.obj)) stop("Need to specify the GpC BSseq object as nome.obj")
    if (!is.null(nome.obj) & !is(nome.obj, "BSseq")) {
      stop("The input nome.obj needs to be a BSseq object of GpC methylation.")
    }
  }
  
  # Pull in the bin annotations
  if (is.null(unique(genome(obj)))) {
    message("No genome found. You need to assign a genome to your bsseq object.")
    stop("As an example: genome(bsseq_obj) <- 'hg19'")
  }
  my.genome <- unique(genome(obj))
  supported.genomes <- c("hg19", "hg38",
                         "mm9", "mm10")
  if (!my.genome %in% supported.genomes) {
    message(my.genome, " is currently not supported.")
    message("Currently supported genomes are: ", supported.genomes)
    stop("Please file an issue here: https://github.com/trichelab/biscuiteer/issues")
  }
  
  ## Temporarily pull bins from user supplied argument
  ## Eventually want to use getBSBinAnnotations(my_genome, my_bin_size)
  if (!is(bins, "GRanges")) stop("The bins need to be a GRanges object.")
  bins <- sort(bins)
  
  ## Summarize methylation level
  ms <- .summarizeMethylation(obj = obj,
                              bins = bins,
                              genome = my.genome)
  
  if (nome) {
    nome.ms <- .summarizeMethylation(obj = nome.obj,
                                     bins = bins,
                                     genome = my.genome)
    ## Compute effective GC
    eff_gc <- .computeEffectiveGC(ms, bins = bins,
                                  nome = TRUE,
                                  meth_summarized.nome = nome.ms)
    return(eff_gc)
  }
  
  ## Compute effective GC
  eff_gc <- .computeEffectiveGC(ms, bins = bins, nome = FALSE)
  return(eff_gc)
}

## Stub for now until pre-computed annotations are uploaded
getBSBinAnnotations <- function() {
  ## This function is a wrapper for biscuiteerDataGet()
  ## TODO: upload the annotations to AnnotationHub()
  return(NULL)
}

## Helper function to summarize methylation levels within a bin
.summarizeMethylation <- function(obj, bins, genome) {
  if (!is(obj, "BSseq")) stop("Input needs to be a BSseq object.")
  if (!is(bins, "GRanges") & is(bins, "AnnotatedDataFrame")) {
    ## This is making the assumption that the input is from QDNAseq bins
    bins <- makeGRangesFromDataFrame(obj, keep.extra.columns = TRUE,
                                     seqnames.field = "chromosome",
                                     start.field = "start",
                                     end.field = "end",
                                     ignore.strand = TRUE)
    seqlevels(bins) <- as.character(unique(seqnames(bins)))
    genome(bins) <- genome
  }
  ## Make sure they are sorted
  bins <- sort(bins)
  
  ## Filter regions that don't have at least 1 read in them
  ## ahead of time
  ## NOTE: this *shouldn't* have to happen if imported with readBiscuit
  ## NOTE: may need to remove this check later
  cov.filt <- getCoverage(obj, type = "Cov", what = "perBase")
  cov.filt <- cov.filt > 0
  obj <- obj[cov.filt,]
  
  ## Subset the bins to regions/chromosomes covered
  bins <- subsetByOverlaps(bins, obj)
  rnames <- as.character(granges(bins))
  
  ## Summarize
  M <- getCoverage(obj, regions = bins,
                   type = "M", what = "perRegionTotal")
  U <- getCoverage(obj, regions = bins, type = "Cov",
                   what = "perRegionTotal") - M
  meth_summarized <- M / (M + U)
  rownames(meth_summarized) <- rnames
  colnames(meth_summarized) <- colnames(obj)
  return(meth_summarized)
}

## Helper function to compute the effective GC content from either
## CpG or GpC binned methylation for CNV sketching
.computeEffectiveGC <- function(meth_summarized,
                                bins, nome = FALSE,
                                meth_summarized.nome = NULL) {
  ## compute the effective GC content for CNV sketching off of CpG coverage
  ## bsgc[i] + (%mCpG[i, j] * cggc[i]) + (%mGpC[i, j] * gcgc[i])
  if (!all(c("bsgc", "cpg_gc") %in% names(mcols(bins)))) {
    stop("Cannot calculate effective GC without bsgc and cpg_gc.")
  }
  
  ## intersect with summarized bins
  meth_summarized.gr <- as(rownames(meth_summarized), "GRanges")
  mcols(meth_summarized.gr) <- meth_summarized
  meth_summarized.gr <- sort(meth_summarized.gr)
  
  bins.sub <- subsetByOverlaps(bins, meth_summarized.gr)
  bins.sub <- sort(bins.sub)
  
  if (nome) {
    meth_summarized.nome.gr <- as(rownames(meth_summarized.nome),
                                  "GRanges")
    mcols(meth_summarized.nome.gr) <- meth_summarized.nome
    meth_summarized.nome.gr <- sort(meth_summarized.nome.gr)
    ## sync up the bins
    meth_summarized.nome.gr <- subsetByOverlaps(meth_summarized.nome.gr,
                                                bins.sub)
    ## sanity check
    stopifnot(all(as.character(granges(meth_summarized.gr)) == as.character(granges(meth_summarized.nome.gr))))
    eff_gc <- mcols(bins.sub)$bsgc +
      (mcols(meth_summarized.gr)$X * mcols(bins.sub)$cpg_gc) +
      (mcols(meth_summarized.nome.gr)$X * mcols(bins.sub)$gpc_gc)
    return(eff_gc)
  }
  eff_gc <- mcols(bins.sub)$bsgc +
    (as.matrix(mcols(meth_summarized.gr)) * mcols(bins.sub)$cpg_gc)
  return(eff_gc)
}

#### Testing #####
bin.annots <- readRDS("~/git_repos/QDNAseq_bin_annotations/data/bs/QDNAseq_hg38_all_BS_bin_annots_precomputed.rds")
## Fix the rownames
# bin.annots <- lapply(bin.annots, function(x) {
#   rownames(x) <- paste0(x@data$chromosome, ":", x@data$start, "-", x@data$end)
#   return(x)
# })
# saveRDS(bin.annots, file = "~/git_repos/QDNAseq_bin_annotations/data/bs/QDNAseq_hg38_all_BS_bin_annots_precomputed.rds")
bin.30kb.hg38 <- bin.annots[[paste0("hg38_", "30", "kbp")]]

library(biscuiteer)
orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
                        package="biscuiteer")
orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
                        package="biscuiteer")
bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
                    merged = FALSE)
## swap the genome
genome(bisc) <- "hg38"

## summarize
bin.30kb.hg38.gr <- as(rownames(bin.30kb.hg38), "GRanges")
mcols(bin.30kb.hg38.gr) <- bin.30kb.hg38@data
bisc.example.effgc <- getEffectiveGC(bisc, bins = bin.30kb.hg38.gr)
