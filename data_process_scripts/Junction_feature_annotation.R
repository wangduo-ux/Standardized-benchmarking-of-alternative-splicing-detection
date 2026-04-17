library(GenomicRanges)
library(dplyr)
library(rtracklayer)

## ------------------------------------------------------------
## Annotate isoform length and exon number for each junction
## ------------------------------------------------------------

## Construct GRanges object for junctions (if needed elsewhere)

junction_gr <- GRanges(
  seqnames = junction$seqnames,          # Chromosome name
  ranges   = IRanges(
    start = junction$start,
    end   = junction$end
  ),
  strand   = junction$strand
)

## Import transcript annotation (GTF format)

gtf_file <- "~/Homo_sapiens.GRCh38.109.gtf"  # input a gene annotation file
transcripts <- import(gtf_file, format = "gtf")

## Extract exon features only
exons <- transcripts[transcripts$type == "exon"]

## Function: derive splice junctions from transcript exons

get_transcript_junctions <- function(exons) {
  
  ## Lists for storing junction GRanges and metadata
  junction_list  <- list()
  metadata_list  <- list()
  
  ## Split exons by transcript
  exons_by_tx <- split(exons, exons$transcript_id)
  
  processed_count <- 0
  tx_num <- 0
  
  for (tx in names(exons_by_tx)) {
    
    tx_exons <- exons_by_tx[[tx]]
    processed_count <- processed_count + 1
    
    ## Report progress every 1000 transcripts
    if (processed_count %% 1000 == 0) {
      message("Processing ", processed_count, " transcripts...")
    }
    
    ## Reverse exon order for negative-strand transcripts
    if (as.character(strand(tx_exons[1])) == "-") {
      tx_exons <- rev(tx_exons)
    }
    
    exon_widths <- width(tx_exons)
    total_length <- sum(exon_widths)
    exon_num <- length(tx_exons)
    
    ## Only transcripts with more than one exon can form junctions
    if (length(tx_exons) > 1) {
      
      start_positions <- start(tx_exons)
      end_positions   <- end(tx_exons)
      
      ## Junction genomic coordinates
      start_junctions <- end_positions[-length(end_positions)] + 1
      end_junctions   <- start_positions[-1] - 1
      
      ## Keep valid junctions only
      valid_junctions <- start_junctions < end_junctions
      
      junction_coords <- data.frame(
        start_junction = start_junctions[valid_junctions],
        end_junction   = end_junctions[valid_junctions]
      )
      
      if (nrow(junction_coords) > 0) {
        
        tx_num <- tx_num + 1
        
        ## Relative junction position within transcript
        junction_pos <- cumsum(exon_widths)[seq_along(tx_exons) - 1]
        junction_pos <- junction_pos[valid_junctions]
        pos <- junction_pos / total_length
        
        ## Adjust position for negative-strand transcripts
        if (as.character(strand(tx_exons[1])) == "-") {
          pos <- 1 - pos
        }
        
        ## Construct GRanges for junctions
        junction <- GRanges(
          seqnames = rep(seqnames(tx_exons[1]), nrow(junction_coords)),
          ranges   = IRanges(
            start = junction_coords$start_junction,
            end   = junction_coords$end_junction
          ),
          strand   = rep(strand(tx_exons[1]), nrow(junction_coords))
        )
        
        ## Store junctions and corresponding metadata
        junction_list[[tx_num]] <- junction
        metadata_list[[tx_num]] <- data.frame(
          transcript_id = rep(tx, nrow(junction_coords)),
          pos           = pos,
          length        = rep(total_length, nrow(junction_coords)),
          exon_num      = rep(exon_num, nrow(junction_coords))
        )
      }
    }
  }
  
  ## Combine all junctions into a single GRanges object
  junctions <- do.call(c, junction_list)
  metadata  <- do.call(rbind, metadata_list)
  
  mcols(junctions) <- metadata
  
  return(junctions)
}

## Generate transcript-level junction annotations

transcript_junctions <- get_transcript_junctions(exons)


## ------------------------------------------------------------
## Annotate GC content for each junction
## ------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)

## Function: calculate GC content around splice junctions

gc_content <- function(junction) {
  
  genome <- BSgenome.Hsapiens.UCSC.hg38
  
  ## Remove non-canonical chromosomes (e.g., contigs with '.')
  junction <- junction[!grepl("\\.", junction$seqnames), ]
  
  ## Construct GRanges object for junctions
  junction_gr <- GRanges(
    seqnames = paste0("chr", junction$seqnames),
    ranges   = IRanges(
      start = junction$start,
      end   = junction$end
    ),
    strand   = junction$strand
  )
  

  ## Internal function: compute GC content for a single junction

  calc_gc_content <- function(junction, genome,
                              upstream = 25,
                              downstream = 25) {
    
    ## Define upstream and downstream regions
    start_upstream   <- start(junction) - upstream
    end_upstream     <- start(junction) - 1
    start_downstream <- end(junction) + 1
    end_downstream   <- end(junction) + downstream
    
    pos_upstream <- GRanges(
      seqnames = seqnames(junction),
      ranges   = IRanges(start = start_upstream, end = end_upstream)
    )
    
    pos_downstream <- GRanges(
      seqnames = seqnames(junction),
      ranges   = IRanges(start = start_downstream, end = end_downstream)
    )
    
    ## Extract genomic sequences
    seq_upstream   <- getSeq(genome, pos_upstream)
    seq_downstream <- getSeq(genome, pos_downstream)
    
    ## Calculate GC content
    gc_upstream   <- sum(letterFrequency(seq_upstream, letters = c("G", "C"))) /
      width(seq_upstream)
    gc_downstream <- sum(letterFrequency(seq_downstream, letters = c("G", "C"))) /
      width(seq_downstream)
    
    gc <- (gc_upstream + gc_downstream) / 2
    
    return(gc)
  }
  
 
  ## Calculate GC content for all junctions

  gc_values <- sapply(seq_along(junction_gr), function(i) {
    calc_gc_content(junction_gr[i], genome)
  })
  
  ## Append GC content to junction metadata
  junction$gc <- gc_values
  
  return(junction)
}


## ------------------------------------------------------------
## Annotate mappability for each junction
## ------------------------------------------------------------

library(GenomicRanges)
library(dplyr)
library(IRanges)

## The mappability of each genomic position is calculated using GenMap (v1.3.0),for example:
## 1       9998    9999    0.00107991
## 1       9999    10000   0.000876424
## 1       10000   10001   0.000714796
## 1       10001   10002   0.000727273
## 1       10002   10003   0.000735294
## 1       10003   10004   0.00071582
mappability_file <- "~/output of GenMap"

## Function: calculate junction flanking mappability

calc_junction_mappability <- function(junction_df,
                                      mappability_file,
                                      flank = 50) {
  
  ## Load mappability data
  map_df <- read.table(
    mappability_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(map_df) <- c("chr", "start", "end", "mappability")
  
  ## Ensure chromosome naming consistency
  map_df$chr <- ifelse(
    grepl("^chr", map_df$chr),
    map_df$chr,
    paste0("chr", map_df$chr)
  )
  
  ## Convert to GRanges
  map_gr <- GRanges(
    seqnames = map_df$chr,
    ranges   = IRanges(start = map_df$start, end = map_df$end),
    score    = map_df$mappability
  )
  

  ## Prepare junction GRanges
  junction_df$chr <- ifelse(
    grepl("^chr", junction_df$chr),
    junction_df$chr,
    paste0("chr", junction_df$chr)
  )
  
  junction_gr <- GRanges(
    seqnames = junction_df$chr,
    ranges   = IRanges(start = junction_df$start, end = junction_df$end),
    strand   = junction_df$strand
  )
  

  ## Define upstream & downstream
  upstream_gr <- GRanges(
    seqnames = seqnames(junction_gr),
    ranges   = IRanges(
      start = start(junction_gr) - flank,
      end   = start(junction_gr) - 1
    ),
    strand = strand(junction_gr)
  )
  
  downstream_gr <- GRanges(
    seqnames = seqnames(junction_gr),
    ranges   = IRanges(
      start = end(junction_gr) + 1,
      end   = end(junction_gr) + flank
    ),
    strand = strand(junction_gr)
  )
  

  ## Overlap with mappability

  calc_mean_map <- function(query_gr, map_gr) {
    hits <- findOverlaps(query_gr, map_gr)
    tapply(
      map_gr$score[subjectHits(hits)],
      queryHits(hits),
      mean,
      na.rm = TRUE
    )
  }
  
  upstream_map   <- calc_mean_map(upstream_gr, map_gr)
  downstream_map <- calc_mean_map(downstream_gr, map_gr)
  

  ## Assemble result
  result <- junction_df
  
  result$map_upstream   <- upstream_map[as.character(seq_len(nrow(result)))]
  result$map_downstream <- downstream_map[as.character(seq_len(nrow(result)))]
  
  ## Average of upstream & downstream
  result$map_mean <- rowMeans(
    cbind(result$map_upstream, result$map_downstream),
    na.rm = TRUE
  )
  
  return(result)
}

## ------------------------------------------------------------
## Example usage
## ------------------------------------------------------------
# junction <- read.table("junction.txt", header = TRUE)
# mappability <- "genome.mappability.bed"
# out <- calc_junction_mappability(junction, mappability, flank = 50)


