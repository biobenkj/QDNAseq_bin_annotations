#' This function creates QDNAseq-like annotations for use with copy-number sketching
#'
#' @param bsgenome Input BSgenome object
#' @param bin.size The bin-size in kilobase-pairs
#' @param chrs The chromosomes to operate on (default is all standard chromosomes)
#' @param nome Whether the data are from NOMe-seq to compute GpC GC content
#' @param parallel Whether to run in parallel
#' @param cores How many cores to use for running in parallel
#'
#' @return A data frame with columns chromosome, start, end, non-N base frequency within the bin, CpG GC content, and GpC GC content.
#' 
#' @import BSgenome
#' @import GenomeInfoDb
#' @import Biostrings
#' @import GenomicRanges
#' @import BiocParallel
#' @importFrom QDNAseq vmesg
#' 
#' @export
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' hg38_bins<- createBSBins(BSgenome.Hsapiens.UCSC.hg38, bin.size = 30, chrs = "chr22")

createBSBins <- function(bsgenome, bin.size, chrs = NULL, blacklist = NULL,
                         nome = FALSE, parallel = FALSE, cores = 2) {
  #make sure the input is what we expect
  if (!is(bsgenome, "BSgenome")) stop("Need a BSgenome to create the bins.")
  if (!is(bin.size, "numeric") & bin.size <= 0) stop("Bin size needs to be an integer value greater than 0.")
  #make sure the blacklist file exists
  if (!is.null(blacklist)) stopifnot(file.exists(blacklist))
  #stop if blacklist null
  if (is.null(blacklist)) stop("Must include a BED file to compute blacklist overlap of bins.")
  
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
  
  #GC content, unconverted by bisulfite
  message("Computing bin-level GC content.")
  wgsgc_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
    gc_freq <- (sum(alphabetFrequency(bp)[2:3]) / sum(alphabetFrequency(bp)[1:4])) #index 2 and 3 are C and G, respectively
    return(gc_freq)
  }, BPPARAM = bpparam))
  
  #bisulfite convert each bin
  message("Computing bin-level bisulfite-converted GC content.")
  gc_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
    bs_convert <- chartr("C", "T", bp)
    gc_freq <- (sum(alphabetFrequency(bs_convert)[2:3]) / sum(alphabetFrequency(bs_convert)[1:4]))
    return(gc_freq)
  }, BPPARAM = bpparam))
  
  
  #compute the frequency/bin (e.g. GC content) for CpG methylation
  #to do this, we need to mask off di- or tri-nucleotides of interest to
  #a non-IUPAC letter, bisulfite convert, turn masked nucleotides back
  message("Computing bin-level bisulfite-converted CpG GC content.")
  if (nome) {
    #need HCG context and not just CpG
    message("NOMe-seq mode - computing bin-level HCG GC content.")
    cpg_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
      #short circuit for fully masked regions
      if (alphabetFrequency(bp)["N"] == nchar(bp)) {
        cpg_freq <- 0
        return(cpg_freq)
      }
      acg_freq <- trinucleotideFrequency(bp)["ACG"]
      tcg_freq <- trinucleotideFrequency(bp)["TCG"]
      nome_cpg_freq <-  acg_freq + tcg_freq
      #bisulfite convert
      bp_convert <- chartr("C", "T", bp)
      #compute GC content
      cpg_freq <- ((alphabetFrequency(bp_convert)[3] + nome_cpg_freq) / sum(alphabetFrequency(bp_convert)[1:4]))
      #check if we divided by 0
      #can happen in fully masked regions
      if (is.nan(cpg_freq)) cpg_freq <- 0
      return(cpg_freq)
    }, BPPARAM = bpparam))
  } else {
    cpg_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
      #short circuit for fully masked regions
      if (alphabetFrequency(bp)["N"] == nchar(bp)) {
        cpg_freq <- 0
        return(cpg_freq)
      }
      cg_freq <- dinucleotideFrequency(bp)["CG"]
      #bisulfite convert
      bp_convert <- chartr("C", "T", bp)
      #compute GC content
      cpg_freq <- ((alphabetFrequency(bp_convert)[3] + cg_freq) / sum(alphabetFrequency(bp_convert)[1:4]))
      #check if we divided by 0
      if (is.nan(cpg_freq)) cpg_freq <- 0
      return(cpg_freq)
    }, BPPARAM = bpparam))
  }
  
  #also compute the GpC context (ignoring GpCpG since it's ambiguous)
  #this is for NOMe-seq
  if (nome) {
    message("NOMe-seq mode - computing bin-level bisulfite-converted GCH GC content.")
    gpc_freq_bin <- unlist(bplapply(genome_bins_with_sequence, function(bp) {
      #effect of methylating all GpCpH dinuleotides
      gpc_freq <- sum(trinucleotideFrequency(bp)["GCA"], trinucleotideFrequency(bp)["GCT"])
      gcg_freq <- trinucleotideFrequency(bp)["GCG"] #ambiguous GpCpG context
      delta <- gpc_freq - gcg_freq
      bp_convert <- chartr("C", "T", bp)
      #compute GC content
      delta_freq <- ((alphabetFrequency(bp_convert)[3] + delta) / sum(alphabetFrequency(bp_convert)[1:4]))
      return(delta_freq)
    }, BPPARAM = bpparam))
  }
  
  #compute the percent overlap with ENCODE blacklist regions
  #use Anshul's new hg38 blacklist regions - https://www.encodeproject.org/files/ENCFF356LFX/
  #mm10 and other mouse builds coming later...
  #make use of QDNAseq's built-in functions
  #well we modified it a bit - could be refactored later if needed
  bins.in <- data.frame(chromosome = as.character(seqnames(genome_bins)),
                        start = start(genome_bins),
                        end = end(genome_bins),
                        bases = as.numeric(base_freq_bin),
                        gc = as.numeric(wgsgc_freq_bin))
  rownames(bins.in) <- as.character(granges(genome_bins))
  blklst.ovlp <- .calculateBlacklist(bins = bins.in, blacklist)

  #compute mappability from bismap: https://bismap.hoffmanlab.org/
  #use k = 100 of multi-read mappability estimates
  #only support human and mouse for now...
  message("Computing MLE mappability estimates for ", genome_build, " at ", bin.size / 1e3, " kb")
  bs_mappability <- .computeMappability(genome_build, genome_bins)
  
  #TODO gather and compute the residuals filter similar to QDNAseq but using the blood WGBS samples from blueprint
  
  if (nome) {
    bin.annots <- data.frame(chromosome = as.character(seqnames(genome_bins)),
                             start = start(genome_bins),
                             end = end(genome_bins),
                             bases = as.numeric(base_freq_bin)*100,
                             gc = as.numeric(wgsgc_freq_bin)*100,
                             bsgc = as.numeric(gc_freq_bin)*100,
                             cpg_gc = as.numeric(cpg_freq_bin)*100, #methylated HCG GC content
                             gpc_gc = as.numeric(gpc_freq_bin)*100, #methylated GCH GC content
                             mappability = as.numeric(bs_mappability)*100,
                             blacklist = as.numeric(blklst.ovlp),
                             stringsAsFactors = FALSE)
  } else {
    bin.annots <- data.frame(chromosome = as.character(seqnames(genome_bins)),
                             start = start(genome_bins),
                             end = end(genome_bins),
                             bases = as.numeric(base_freq_bin)*100, #percentage of non-N bases
                             gc = as.numeric(wgsgc_freq_bin)*100, #GC content, no conversion by bisulfite
                             bsgc = as.numeric(gc_freq_bin)*100, #all C converted by bisulfite
                             cpg_gc = as.numeric(cpg_freq_bin)*100, #methylated CpGs GC content
                             mappability = as.numeric(bs_mappability)*100, #bismap mappability
                             blacklist = as.numeric(blklst.ovlp), #percent bases overlap in blacklist region
                             stringsAsFactors = FALSE)
  }
  
  message("Done.")
  return(bin.annots)
}

.computeMappability <- function(gb, bins) {
  #read in the appropriate genome file and compute the mappability
  #at the desired resolution - bin.size
  #these are massive...
  gb.k100 <- readRDS(paste0("/secondary/projects/triche/ben_projects/resources/umap_and_bismap/bismap/", gb,
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

#This is a modified function from QDNAseq
.calculateBlacklist <- function(bins, bedFiles, ...,
                               verbose=getOption("QDNAseq::verbose", TRUE)) {
  
  ## Detect defunct usage of argument 'ncpus'
  if ("ncpus" %in% names(list(...))) {
    .Defunct(msg="Argument 'ncpus' of calculateBlacklist() is defunct. Use future::plan() instead.")
  }
  
  oopts <- options("QDNAseq::verbose"=verbose)
  on.exit(options(oopts))
  
  QDNAseq:::vmsg("Calculating overlaps per bin with BED files \n    ", paste(bedFiles,
                                                                   collapse="\n    "), "\n    ...", appendLF=FALSE)
  
  beds <- list()
  for (bed in bedFiles)
    beds[[bed]] <- read.table(bed, sep="\t", as.is=TRUE)
  combined <- beds[[1L]]
  if (length(beds) >= 2L)
    for (i in 2:length(beds))
      combined <- rbind(combined, beds[[i]])
  combined <- combined[, 1:3]
  colnames(combined) <- c("chromosome", "start", "end")
  #combined$chromosome <- sub("^chr", "", combined$chromosome)
  combined <- combined[combined$chromosome %in% unique(bins$chromosome), ]
  combined <- combined[!is.na(combined$chromosome), ]
  combined$start <- combined$start + 1
  ## define correct sorting order of chromosomes as the order in which they
  ## are in the bins
  chromosomes <- unique(bins$chromosome)
  chromosomeOrder <- factor(combined$chromosome, levels=chromosomes,
                            ordered=TRUE)
  combined <- combined[order(chromosomeOrder, combined$start), ]
  joined <- data.frame()
  prev <- combined[1L,]
  # Sanity check
  stopifnot(nrow(combined) >= 2L);
  for (i in 2:nrow(combined)) {
    if (combined[i, "chromosome"] != prev$chromosome ||
        combined[i, "start"] > (prev$end + 1)) {
      joined <- rbind(joined, prev)
      prev <- combined[i,]
    } else {
      prev$end <- max(prev$end, combined[i, "end"])
    }
  }
  joined <- rbind(joined, prev)
  overlap.counter <- function(x, joined) {
    chr <- x["chromosome"]
    start <- as.integer(x["start"])
    end <- as.integer(x["end"])
    overlaps <- joined[joined$chromosome == chr &
                         joined$start      <= end &
                         joined$end        >= start, ]
    bases <- 0
    for (i in rownames(overlaps))
      bases <- bases + min(end, overlaps[i, "end"]) -
      max(start, overlaps[i, "start"]) + 1
    bases / (end - start + 1) * 100
  }
  blacklist <- apply(bins, MARGIN=1L, FUN=overlap.counter, joined=joined)
  QDNAseq:::vmsg()
  blacklist
}
