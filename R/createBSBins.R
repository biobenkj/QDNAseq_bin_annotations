#' This function creates QDNAseq-like annotations for use with copy-number sketching
#'
#' @param bsgenome Input BSgenome object
#' @param bin.size The bin-size in kilobase-pairs
#' @param chrs The chromosomes to operate on (default is all standard chromosomes)
#'
#' @return A data frame with columns chromosome, start, end, non-N base frequency within the bin, CpG GC content, and GpC GC content.
#' 
#' @import BSgenome
#' @import GenomeInfoDb
#' @import Biostrings
#' @import GenomicRanges
#' @import BiocParallel
#' 
#' @export
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' hg38_bins<- createBSBins(BSgenome.Hsapiens.UCSC.hg38, bin.size = 30, chrs = "chr22")

createBSBins <- function(bsgenome, bin.size, chrs = NULL,
                         nome = FALSE, parallel = FALSE, cores = 2) {
  #make sure the input is what we expect
  if (!is(bsgenome, "BSgenome")) stop("Need a BSgenome to create the bins.")
  if (!is(bin.size, "numeric") & bin.size <= 0) stop("Bin size needs to be an integer value greater than 0.")
  
  #now let's reinvent the wheel and round off a few edges for BS-seq
  #this is adapted from QDNAseq::createBins
  
  #chrs to work on
  #just do standard chromosomes for now
  if (is.null(chrs)) chrs <- standardChromosomes(bsgenome)
  
  #make sure they are all there...
  if (!all((chrs %in% seqnames(bsgenome)))) stop("Supplied chromosomes not found in the BSgenome object.")
  
  #now what are we working on?
  genome_build <- unique(genome(bsgenome))
  if (!(genome_build %in% c("hg19", "hg38", "mm9", "mm10"))) {
    stop("Only supports human and mouse genomes for now...")
  }
  
  #bin the genome
  message("Binning the genome at ", bin.size, " kb")
  #convert to kbp similar to QDNAseq::createBins
  bin.size <- bin.size * 1e3
  genome_bins <- tileGenome(seqlengths = seqlengths(bsgenome)[chrs],
                            tilewidth = bin.size,
                            cut.last.tile.in.chrom = TRUE)
  
  #are we running in parallel?
  if (parallel) {
    bpparam <- MulticoreParam(workers = cores)
  } else {
    bpparam <- SerialParam()
  }
  
  #overlap with the bsgenome object
  #this will return a stringset object that will correspond to the bases within each bin
  genome_bins_with_sequence <- BSgenome::getSeq(bsgenome, genome_bins)
  #add names; this will be helpful for matching things up later
  names(genome_bins_with_sequence) <- as.character(granges(genome_bins))
  
  #compute the number of *non-N* bases within each bin
  message("Computing bin-level non-N base content.")
  base_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
    base_freq <- alphabetFrequency(bp)
    #get the summed frequency of A,C,G,T
    base_sum <- sum(base_freq[1:4])
    base_bin_freq <- (base_sum / nchar(bp))
    return(base_bin_freq)
  }, BPPARAM = bpparam))
  
  #bisulfite convert each bin
  message("Computing bin-level bisulfite-converted GC content.")
  gc_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
    bs_convert <- chartr("C", "T", bp)
    gc_freq <- (sum(alphabetFrequency(bp)[2:3]) / nchar(bp))
    return(gc_freq)
  }, BPPARAM = bpparam))
  
  
  #compute the frequency/bin (e.g. GC content) for CpG methylation
  message("Computing bin-level CpG GC content.")
  if (nome) {
    #need HCG context and not just CpG
    cpg_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
      acg_freq <- trinucleotideFrequency(bp)["ACG"]
      tcg_freq <- trinucleotideFrequency(bp)["TCG"]
      cpg_freq <-  (acg_freq + tcg_freq) / nchar(bp)
      return(cpg_freq)
    }, BPPARAM = bpparam))
  } else {
    cpg_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
      cpg_freq <- (dinucleotideFrequency(bp)["CG"] / nchar(bp))
      return(cpg_freq)
    }, BPPARAM = bpparam))
  }
  
  #also compute the GpC context (ignoring GpCpG since it's ambiguous)
  #this is for NOMe-seq
  if (nome) {
    message("Computing bin-level GpC GC content.")
    gpc_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
      gpc_freq <- sum(trinucleotideFrequency(bp)["GCA"], trinucleotideFrequency(bp)["GCT"])
      gcg_freq <- trinucleotideFrequency(bp)["GCG"]
      delta <- gpc_freq - gcg_freq
      delta_freq <- (delta / nchar(bp))
      return(delta_freq)
    }, BPPARAM = bpparam))
  }
  
  #compute mappability from bismap: https://bismap.hoffmanlab.org/
  #use k = 100 of multi-read mappability estimates
  #only support human and mouse for now...
  message("Computing MLE mappability estimates for ", genome_build, " at ", bin.size / 1e3, " kb")
  bs_mappability <- .computeMappability(genome_build, genome_bins)
  
  #reconstruct the data frame for output
  if (nome) {
    bin.annots <- data.frame(chromosome = as.character(seqnames(genome_bins)),
                             start = start(genome_bins),
                             end = end(genome_bins),
                             bases = as.numeric(base_freq_bin),
                             bsgc = as.numeric(gc_freq_bin),
                             cpg_gc = as.numeric(cpg_freq_bin),
                             gpc_gc = as.numeric(gpc_freq_bin),
                             mappability = as.numeric(bs_mappability),
                             stringsAsFactors = FALSE)
  } else {
    bin.annots <- data.frame(chromosome = as.character(seqnames(genome_bins)),
                             start = start(genome_bins),
                             end = end(genome_bins),
                             bases = as.numeric(base_freq_bin),
                             bsgc = as.numeric(gc_freq_bin),
                             cpg_gc = as.numeric(cpg_freq_bin),
                             mappability = as.numeric(bs_mappability),
                             stringsAsFactors = FALSE)
  }
  
  message("Done.")
  return(bin.annots)
}

.computeMappability <- function(gb, bins) {
  #read in the appropriate genome file and compute the mappability
  #at the desired resolution - bin.size
  #these are massive...
  gb.k100 <- readRDS(paste0("~/secondary_projects/triche/ben_projects/resources/umap_and_bismap/bismap/", gb,
                            "/", gb, "_k100_bp_res_multi_read_granges.rds"))
  #compute the MLE of the bisulfite mappability
  #find overlaps
  ids <- findOverlaps(gb.k100, bins, select = "first")
  #create a vector of zeros
  zvec <- rep(0, length(bins))
  mle_map <- tapply(as.matrix(mcols(gb.k100)$score),
                    INDEX = ids,
                    FUN = mean)
  zvec[as.numeric(names(mle_map))] <- mle_map
  return(zvec)
}
