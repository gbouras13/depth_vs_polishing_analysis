#load library
# BiocManager::install("Gviz")
# BiocManager::install("GenomicRanges")
# BiocManager::install("GenomicAlignments")
library(GenomicAlignments)
library(GenomicRanges)
library(Gviz) 

fasta_file <- "ATCC_25922_plasmid_4.fasta"

bam <- "high_conf_alignment_sorted.bam"

# minimap no secondary alignments 
# https://github.com/lh3/minimap2/issues/98
# https://github.com/ivanek/Gviz/issues/21

.import.bam.alignments.ignore.secondary <- function(file, selection) {
  indNames <- c(sub("\\.bam$", ".bai", file), paste(file, "bai", sep = "."))
  index <- NULL
  for (i in indNames) {
    if (file.exists(i)) {
      index <- i
      break
    }
  }
  if (is.null(index)) {
    stop(
      "Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t",
      "library(Rsamtools)\n\tindexBam(\"", file, "\")"
    )
  }
  pairedEnd <- parent.env(environment())[["._isPaired"]]
  if (is.null(pairedEnd)) {
    pairedEnd <- TRUE
  }
  flag <- parent.env(environment())[["._flag"]]
  if (is.null(flag)) {
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
  }
  flag <- bamFlagAND(flag, scanBamFlag(isSecondaryAlignment = FALSE))
  bf <- BamFile(file, index = index, asMates = pairedEnd)
  param <- ScanBamParam(which = selection, what = scanBamWhat(), tag = "MD", flag = flag)
  reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets)) scanBam(bf, param = param)[[1]] else list()
  md <- if (is.null(reads$tag$MD)) rep(as.character(NA), length(reads$pos)) else reads$tag$MD
  if (length(reads$pos)) {
    layed_seq <- sequenceLayer(reads$seq, reads$cigar)
    region <- unlist(bamWhich(param), use.names = FALSE)
    ans <- stackStrings(layed_seq, start(region), end(region), shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
    names(ans) <- seq_along(reads$qname)
  } else {
    ans <- DNAStringSet()
  }
  return(GRanges(
    seqnames = if (is.null(reads$rname)) character() else reads$rname,
    strand = if (is.null(reads$strand)) character() else reads$strand,
    ranges = IRanges(start = reads$pos, width = reads$qwidth),
    id = if (is.null(reads$qname)) character() else reads$qname,
    cigar = if (is.null(reads$cigar)) character() else reads$cigar,
    mapq = if (is.null(reads$mapq)) integer() else reads$mapq,
    flag = if (is.null(reads$flag)) integer() else reads$flag,
    md = md, seq = ans,
    isize = if (is.null(reads$isize)) integer() else reads$isize,
    groupid = if (pairedEnd) if (is.null(reads$groupid)) integer() else reads$groupid else seq_along(reads$pos),
    status = if (pairedEnd) {
      if (is.null(reads$mate_status)) factor(levels = c("mated", "ambiguous", "unmated")) else reads$mate_status
    } else {
      rep(
        factor("unmated", levels = c("mated", "ambiguous", "unmated")),
        length(reads$pos)
      )
    }
  ))
}



options(ucscChromosomeNames=FALSE)
seqset <- readDNAStringSet(fasta_file)

sTrack <- SequenceTrack(seqset,
                        genome = "plasmid_4",
                        chromosome = "plasmid_4")

# https://github.com/ivanek/Gviz/issues/16
alTrack <- AlignmentsTrack(
  bam, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "plasmid_4",
  chromosome = "plasmid_4",
  isPaired = TRUE,
  fill.coverage= "white",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion="#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= TRUE,
  name="True Positive", 
  type=c("coverage", "pileup"), stacking = "squish")




# grDevices::dev.size("in")
plot_dir <- "Plots"
dir.create(plot_dir )

high_conf <- plotTracks(  c(sTrack,alTrack),
             panel.only = TRUE,
             from = 45, 
             to = 55, 
             chromosome = "plasmid_4", 
             genome = "plasmid_4", 
             # transcriptAnnotation = "transcript",
             background.title="black",
             col.axis="white",
             col.title= "white")

jpeg(file = paste(plot_dir,"high_conf.jpeg" , sep = "/" ),units = "in", 
     width = 6, 
     height = 5, bg = "transparent", res = 600)

plotTracks(  c(sTrack,alTrack),
             panel.only = TRUE,
             from = 45, 
             to = 55, 
             chromosome = "plasmid_4", 
             genome = "plasmid_4", 
             # transcriptAnnotation = "transcript",
             background.title="black",
             col.axis="white",
             col.title= "white")

dev.off()




###########################
# diff
####################33

bam <- "diff_alignment_sorted.bam"


# https://github.com/ivanek/Gviz/issues/16
alTrack <- AlignmentsTrack(
  bam, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "plasmid_4",
  chromosome = "plasmid_4",
  isPaired = TRUE,
  fill.coverage= "white",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion="#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= TRUE,
  name="False Positive", 
  type=c("coverage", "pileup"), stacking = "squish")




jpeg(file = paste(plot_dir,"diff.jpeg" , sep = "/" ),units = "in", 
     width = 6, 
     height = 5, bg = "transparent", res = 600)



plotTracks(  c(sTrack,alTrack),
             panel.only = TRUE,
             from = 45, 
             to = 55, 
             chromosome = "plasmid_4", 
             genome = "plasmid_4", 
             # transcriptAnnotation = "transcript",
             background.title="black",
             col.axis="white",
             col.title= "white")

dev.off()




###########################
# ratio
####################33

bam <- "ratio_alignment_sorted.bam"


# https://github.com/ivanek/Gviz/issues/16
alTrack <- AlignmentsTrack(
  bam, importFunction = .import.bam.alignments.ignore.secondary,
  genome = "plasmid_4",
  chromosome = "plasmid_4",
  isPaired = TRUE,
  fill.coverage= "white",
  col.coverage= "#091833",
  lwd.coverage=2,
  alpha.reads= 1,
  alpha.mismatch=1,
  col.reads= "#c0c2ce",
  fill.reads= "#c0c2ce",
  showIndels = TRUE,
  col.deletion="#f20505",
  col.gap="white",
  col.insertion= "#c0c2ce",
  lwd.deletions= 20,
  lwd.gap= 0,
  showMismatches= TRUE,
  name="Ratio", 
  type=c("coverage", "pileup"), stacking = "squish")




jpeg(file = paste(plot_dir,"ratio" , sep = "/" ),units = "in", 
     width = 6, 
     height = 5, bg = "transparent", res = 600)



plotTracks(  c(sTrack,alTrack),
             panel.only = TRUE,
             from = 45, 
             to = 55, 
             chromosome = "plasmid_4", 
             genome = "plasmid_4", 
             # transcriptAnnotation = "transcript",
             background.title="black",
             col.axis="white",
             col.title= "white")

dev.off()

