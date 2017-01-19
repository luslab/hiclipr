# =============================================================================
# Get transcriptome sequence from data.table
# =============================================================================

#' Get sequence
#'
#' Get sequence from annotated hybrid data.table
#'
#' @param hybrid.dt hybrid data.table
#'
#' @return hybrid data.table with sequences
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @export

GetSequence <- function(hybrid.dt) {
  
  # First set keys
  hybrid.dt[, ensembl_gene_id := L_seqnames]
  setkey(longest.rna.dt, ensembl_gene_id)
  setkey(hybrid.dt, ensembl_gene_id)
  
  # Then merge  (L == R as intragenic) and get sequences. All on +ve strand as intragenic and because of alignment settings
  seq.dt <- longest.rna.dt[, .(ensembl_gene_id, sequence)][hybrid.dt]
  seq.dt[, L_sequence := substr(sequence, start = L_start, stop = L_end)]
  seq.dt[, R_sequence := substr(sequence, start = R_start, stop = R_end)]
  
  # Remove extra temporary columns
  seq.dt[, c("ensembl_gene_id", "sequence") := NULL]
  
  return(seq.dt)
  
}

# =============================================================================
# Run RNA Hybrid using info from hybrid data table
# =============================================================================

#' Run RNA hybrid
#'
#' Run RNA Hybrid using info from hybrid data table
#'
#' @param hybrid.dt hybrid data.table
#' @param file RNAhybrid output filename
#'
#' @return annotated hybrid data.table with sequences
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @export

RunRNAHybrid <- function(hybrid.dt, file) {
  
  results.dir <- "results/rnahybrid"
  # need to rm file first, otherwise appended!
  if(file.exists(file.path(results.dir, file))) {
    cat("The file already exists. Aborting.\n")
  } else {
    
    # First create FASTA for each hybrid read and pass to RNAhybrid in turn
    # (Cannot pass in bulk as otherwise compares every combination!)
    hybrid.dt[, .(id, L_sequence, R_sequence)]
    hybrid.dt[, fastaL := paste0(">", id, "_L\n", L_sequence)]
    hybrid.dt[, fastaR := paste0(">", id, "_R\n", R_sequence)]
    
    for (n in 1:nrow(hybrid.dt)) {
      
      writeLines(hybrid.dt$fastaL[n], file.path(results.dir, "RNAhybridL.fa"))
      writeLines(hybrid.dt$fastaR[n], file.path(results.dir, "RNAhybridR.fa"))
      
      cmd <- paste0("RNAhybrid -c -s 3utr_human -m 1000 -n 1000 -q ", results.dir, "/RNAhybridL.fa -t ", results.dir, "/RNAhybridR.fa >> ", results.dir, "/", file)
      system(cmd)
      
    }
  }
}

# =============================================================================
# Get start positions for each arm and width for longest duplex
# =============================================================================

#' Find Duplex starts
#'
#' Get start positions for each arm and width for longest duplex. Called by FindDuplex
#'
#' @param RNAhybrid RNA Hybrid row
#' @return list of duplex starts and length
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @export

FindDuplexPositions <- function(RNAhybrid) {
  
  # If RunRNAHybrid used then:
  # t == target == R
  # q == query == L
  
  position.index <- as.integer(RNAhybrid[7]) # get position index
  
  output.dt <- data.table(t.u = strsplit(RNAhybrid[8], "")[[1]],
                          t.p = strsplit(RNAhybrid[9], "")[[1]],
                          q.p = strsplit(RNAhybrid[10], "")[[1]],
                          q.u = strsplit(RNAhybrid[11], "")[[1]]) # reproduces the RNAhybrid output as a dataframe with one nucleotide per cell
  output.dt[t.p != " ", duplex := 1] # identify duplex as "1"
  output.dt[is.na(duplex), duplex := 0][]
  output.dt[t.u != " " | t.p != " ", target.index := position.index:(position.index + .N - 1)] # target index (allowing for kinks in query); -1 as 1-based coordinate system
  output.dt[q.u != " " | q.p != " ", query.index := .N:1][] # query index (allowing for kinks in target)
  
  duplex <- paste(output.dt[, duplex], collapse = "") # collapses the duplex column into a string
  duplex.l <- gregexpr("1*", duplex) # gets runs of 1
  longest.duplex <- max(attr(duplex.l[[1]], "match.length")) # gets longest run of 1
  duplex.index <- min(duplex.l[[1]][attr(duplex.l[[1]], "match.length") == longest.duplex]) # identifies which position has longest run of 1, min means selects first...
  
  target.start <- output.dt[duplex.index, target.index]
  query.start <- output.dt[(duplex.index + longest.duplex - 1), query.index] # Need to get from other end as 3' to 5'; -1 as 1-based
  
  results <- list(target.start = target.start,
                  query.start = query.start,
                  length = longest.duplex)
  
  return(results)
  
}

# =============================================================================
# Find duplices by using duplex position to update datatable
# =============================================================================

#' Find Duplex starts
#'
#' Get start positions for each arm and width for longest duplex.
#'
#' @param hybrid.dt hybrid data.table
#' @param RNAhybrid.out path to RNA Hybrid output
#' @return data.table of duplexes
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @export

FindDuplex <- function(hybrid.dt, RNAhybrid.out) {
  
  results.dir <- "results/rnahybrid"
  
  # Load RNAhybrid.out
  RNAhybrid <- read.table(file.path(results.dir, RNAhybrid.out), sep = ":", header = FALSE, stringsAsFactors = FALSE)
  RNAhybrid <- data.table(RNAhybrid)
    if(nrow(RNAhybrid[V5 == 0]) > 0) {
    print(paste("RNAhybrid did not find a duplex for hybrid", gsub("_R", "", RNAhybrid[V5 == 0, V1])))
  }
  RNAhybrid <- RNAhybrid[V5 != 0]

  # Get duplex positions
  duplex.positions <- apply(RNAhybrid, 1, FindDuplexPositions)
  duplex.positions <- rbindlist(duplex.positions) # convert from a list of lists (with target.start, query.start and length) to a datatable
  duplex.positions[, id := gsub("_R", "", RNAhybrid$V1)] # add id from RNAhybrid.out
  setkey(duplex.positions, id)
  
  # Merge with non-redundant hybrid datatable
  setkey(hybrid.dt, id)
  results.dt <- merge(hybrid.dt, duplex.positions)
  
  # Adjust positions of duplex with RNA hybrid positions
  results.dt[, L_start := L_start + query.start - 1] # -1 for 1 based
  results.dt[, L_width := length]
  results.dt[, L_end := L_start + L_width - 1]
  results.dt[, R_start := R_start + target.start - 1]
  results.dt[, R_width := length]
  results.dt[, R_end := R_start + R_width - 1]
  
  # Remove workings
  results.dt[, c("target.start", "query.start", "length", "L_sequence", "R_sequence", "fastaL", "fastaR") := NULL]
  
  return(results.dt)
  
}
