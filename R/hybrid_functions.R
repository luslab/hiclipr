# =============================================================================
# Load hybrid reads from BAM and identify valid reads (i.e both L and R arms aligned)
# =============================================================================

#' Loads hybrid BAM file
#'
#' Loads hybrid reads from BAM and identify valid reads (i.e both L and R arms aligned)
#'
#' @param bam.file Path to BAM file
#' @return data.table of valid hybrids
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @importFrom GenomicAlignments readGAlignments njunc
#' @export

LoadBAM <- function(bam.file) {
  
  processed.dir <- "results/processed_data"
  mapped.dir <- "results/mapped"
  
  # Load BAM file
  cat("Loading BAM...\n")
  
  bam <- readGAlignments(file.path(mapped.dir, bam.file), use.names = TRUE)
  bam <- bam[njunc(bam) == 0] # ensure no splice junctions (STAR)
  cat("There are", length(bam), "reads mapped.\n")
  # ? add in check to ensure no duplicate reads
  
  # Convert to data table, tidy tRNA and rRNA seqnames and ensure no MT seqnames
  bam.df <- as.data.frame(bam)
  bam.df$seqnames <- as.character(bam.df$seqnames)
  bam.df$qname <- row.names(bam.df)
  bam.dt <- data.table(bam.df)
  setkey(bam.dt, seqnames)
  bam.dt[grep("NR_023363.1", seqnames), seqnames := "rRNA_5S"]
  bam.dt[grep("NR_003285.2", seqnames), seqnames := "rRNA_5.8S"]
  bam.dt[grep("NR_003287.2", seqnames), seqnames := "rRNA_28S"]
  bam.dt[grep("NR_003286.2", seqnames), seqnames := "rRNA_18S"]
  bam.dt[grep("NR_046235.1", seqnames), seqnames := "rRNA_45S"]
  bam.dt[grep("tRNA", seqnames), seqnames := gsub("Homo_sapiens_", "", seqnames)]
  bam.dt <- bam.dt[seqnames != "MT"] # remove any aligning to MT
  bam.dt <- bam.dt[grep("GL", seqnames, invert = TRUE)]
  
  # Prepare for evaluating barcodes
  bam.dt[, read := sapply(strsplit(bam.dt$qname, "_"), "[[", 1)]
  bam.dt[, arm := sapply(strsplit(bam.dt$qname, "_"), "[[", 2)]
  left.dt <- bam.dt[arm == "L"] # split out arms
  right.dt <- bam.dt[arm == "R"]
  setkey(left.dt, read)
  setkey(right.dt, read)
  hybrids.dt <- merge(left.dt, right.dt) # merged back as hybrid reads (both L and R have to be mapped)
  setkey(hybrids.dt, read)
  cat("There are", nrow(hybrids.dt), "valid hybrid reads.\n")
  
  # ensure L arm comes before R arm
  hybrids.dt[, L_start := ifelse(start.x > start.y, start.y, start.x)] # ensure L arm comes before R arm
  hybrids.dt[, R_start := ifelse(start.x > start.y, start.x, start.y)]
  
  hybrids.dt[, L_seqnames := ifelse(L_start == start.x, seqnames.x, seqnames.y)] # swap the rest over if start swapped
  hybrids.dt[, L_strand := ifelse(L_start == start.x, as.character(strand.x), as.character(strand.y))]
  hybrids.dt[, L_width := ifelse(L_start == start.x, width.x, width.y)] # qwidth is query width, width is width that was actually aligned
  hybrids.dt[, L_qname := ifelse(L_start == start.x, qname.x, qname.y)]
  
  hybrids.dt[, R_seqnames := ifelse(L_start == start.x, seqnames.y, seqnames.x)]
  hybrids.dt[, R_strand := ifelse(L_start == start.x, as.character(strand.y), as.character(strand.x))]
  hybrids.dt[, R_width := ifelse(L_start == start.x, width.y, width.x)]
  hybrids.dt[, R_qname := ifelse(L_start == start.x, qname.y, qname.x)]
  
  results.dt <- hybrids.dt[, .(read, L_seqnames, L_start, L_strand, L_width, L_qname,
                               R_seqnames, R_start, R_strand, R_width, R_qname)]
  
  return(results.dt)
  
}

# =============================================================================
# Match arms from BAM and use random barcode to remove PCR duplicates
# =============================================================================

#' Evaluate Barcodes
#'
#' Matches arms from BAM and uses random barcode to remove PCR duplicates
#'
#' @param bam.dt BAM data table (from LoadBAM)
#' @param barcodeindex Path to barcode index table
#' @return data.table of unique hybrids
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @export

EvaluateBarcodes <- function(bam.dt, barcodeindex) {
  
  processed.dir <- "results/processed_data"
  
  # Load barcodes
  cat("Loading barcodes...\n")
  barcodes.dt <- fread(file.path(processed.dir, barcodeindex))
  barcodes.dt[, read := as.character(read)]
  setkey(barcodes.dt, read)
  setkey(bam.dt, read)
  
  # Add barcodes to BAM
  hybrids.dt <- barcodes.dt[bam.dt]
  
  # Add number of PCR duplicates
  # PCR duplicate defined as same barcode, chr, start, strand
  hybrids.dt[, duplicates := .N, by = .(random_barcode,
                                        L_seqnames, L_start, L_strand,
                                        R_seqnames, R_start, R_strand)]
  cat("This is the frequency table of PCR duplicates:")
  print(table(hybrids.dt$duplicates))
  
  # Remove PCR duplicates
  cat("Collapsing PCR duplicates.\n")
  hybrids.dt <- hybrids.dt[, .(read, random_barcode, duplicates, L_seqnames, L_start, L_strand, L_width, L_qname,
                               R_seqnames, R_start, R_strand, R_width, R_qname)]
  
  setkey(hybrids.dt, random_barcode, L_seqnames, L_start, R_seqnames, R_start, L_strand, R_strand)
  unique_hybrids.dt <- unique(hybrids.dt, by = c("random_barcode", "L_seqnames", "L_start", "R_seqnames", "R_start", "L_strand", "R_strand")) # == YS BED collapsed
  results.dt <- unique_hybrids.dt[, .(read, L_seqnames, L_start, L_strand, L_width, L_qname,
                                      R_seqnames, R_start, R_strand, R_width, R_qname)]
  cat("There are", nrow(results.dt), "unique hybrid reads.\n")
  
  # Return data.table of L and R unique hybrids - could change this to add more info if needed
  return(results.dt)
  
}

# =============================================================================
# Convert datatable to GRanges keeping metadata
# =============================================================================

#' Convert data.table to GRanges
#'
#' Convert data.table to GRanges keeping metadata
#'
#' @param hybrid.dt hybrid data.table to convert
#' @param arm arm to convert (all, left, right)
#' @return GRanges object
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @export


ConvertToGRanges <- function(hybrid.dt, arm = c("all", "left", "right")) {
  
  # Table has to have the following columns: "L_seqnames", "L_start", "L_strand", "L_width", "R_seqnames",  "R_start", "R_strand", "R_width"
  # This does not work for datatable with only 1 row.
  # First do L
  
  # Create dataframe
  L.gr <- GRanges(seqnames = Rle(hybrid.dt$L_seqnames),
                  ranges = IRanges(start = hybrid.dt$L_start, width = hybrid.dt$L_width),
                  strand = Rle(hybrid.dt$L_strand))
  # Add metadata
  L.names <- colnames(hybrid.dt)[!colnames(hybrid.dt) %in% c("L_seqnames", "L_start", "L_end", "L_strand", "L_width")] # filter out columns used to create GRanges
  L.names <- L.names[grep("R_", L.names, invert = TRUE)] # filter out R arm metadata
  df <- as.data.frame(hybrid.dt, stringsAsFactors = FALSE) # convert data.table to dataframe for lapply
  metadata <- sapply(L.names, function(x) df[, colnames(df) %in% x]) # get matrix of each metadata column
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE) # convert to dataframe without strings as factors
  mcols(L.gr) <- metadata # add to GRanges
  
  # Then do R
  
  # Create dataframe
  R.gr <- GRanges(seqnames = Rle(hybrid.dt$R_seqnames),
                  ranges = IRanges(start = hybrid.dt$R_start, width = hybrid.dt$R_width),
                  strand = Rle(hybrid.dt$R_strand))
  # Add metadata
  R.names <- colnames(hybrid.dt)[!colnames(hybrid.dt) %in% c("R_seqnames",  "R_start", "R_end", "R_strand", "R_width")] # filter out columns used to create GRanges
  R.names <- R.names[grep("L_", R.names, invert = TRUE)] # filter out L arm metadata
  df <- as.data.frame(hybrid.dt, stringsAsFactors = FALSE) # convert data.table to dataframe for lapply
  metadata <- sapply(R.names, function(x) df[, colnames(df) %in% x]) # get matrix of each metadata column
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE) # convert to dataframe without strings as factors
  mcols(R.gr) <- metadata # add to GRanges
  
  if(arm == "all") {
    combined_seqlevels <- unique(c(seqlevels(L.gr), seqlevels(R.gr)))
    seqlevels(L.gr) <- combined_seqlevels # need to set combined seqlevels for merge
    seqlevels(R.gr) <- combined_seqlevels
    names(mcols(L.gr)) <- gsub("L_", "", L.names) # remove L_ so has same colnames
    names(mcols(R.gr)) <- gsub("R_", "", R.names)
    L.gr$arm <- "L" # identify arm
    R.gr$arm <- "R"
    gr <- c(L.gr, R.gr)
    return(gr)
  } else if(arm == "left") {
    L.gr$arm <- "L"
    return(L.gr)
  } else {
    R.gr$arm <- "R"
    return(R.gr)
  }
}

# =============================================================================
# Annotate arm GRanges
# =============================================================================

#' Annotate GRanges arm
#'
#' Annotates GRanges arm (created using ConvertToGRanges) using rna.tc.gr
#'
#' @param arm.gr GRanges to annotate
#' @return Annotated GRanges object
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @export

AnnotateHybridArm <- function(arm.gr) {
  
  arm.gr$biotype <- as.character(NA) # need these otherwise will not fill later
  arm.gr$region <- as.character(NA)
  arm.gr$gene <- as.character(NA)
  
  # 1. Annotated those that mapped to transcriptome
  # Annotate based on start
  start.gr <- resize(arm.gr, width = 1, fix = "start")
  
  seqlevels.all <- unique(c(seqlevels(start.gr), seqlevels(rna.tc.gr))) # have not subset out as should not overlap with rna.tc.gr anyway
  seqlevels(start.gr) <- seqlevels.all
  seqlevels(rna.tc.gr) <- seqlevels.all
  
  overlaps <- findOverlaps(start.gr, rna.tc.gr)
  arm.gr$biotype[queryHits(overlaps)] <- rna.tc.gr$biotype[subjectHits(overlaps)] # can annotate arm directly
  arm.gr$region[queryHits(overlaps)] <- rna.tc.gr$region[subjectHits(overlaps)]
  arm.gr$gene[queryHits(overlaps)] <- as.character(seqnames(rna.tc.gr)[subjectHits(overlaps)]) # this should be the same as seqnames(arm.gr), but kept this way for consistency
  
  # Annotate negative strand of ENSG as intergenic - not needed if --noRC as alignment setting
  arm.gr$biotype[grepl("^ENSG*", as.character(seqnames(arm.gr))) & as.logical(strand(arm.gr) == "-")] <- "intergenic" # need as.logical as otherwise returns an Rle which you cannot use to subset
  arm.gr$region[grepl("^ENSG*", as.character(seqnames(arm.gr))) & as.logical(strand(arm.gr) == "-")] <- "INTERGENIC"
  arm.gr$gene[grepl("^ENSG*", as.character(seqnames(arm.gr))) & as.logical(strand(arm.gr) == "-")] <- "intergenic"
  
  # 2. Annotate those that mapped to genome
  # First get those that overlap with the ranges of known
  seqlevels.all <- unique(c(seqlevels(start.gr), seqlevels(rna.gc.gr)))
  seqlevels(start.gr) <- seqlevels.all
  seqlevels(rna.gc.gr) <- seqlevels.all
  
  overlaps <- findOverlaps(start.gr, rna.gc.gr)
  arm.gr$biotype[queryHits(overlaps)] <- rna.gc.gr$biotype[subjectHits(overlaps)] # can annotate arm.gr directly as same order as start.gr
  arm.gr$region[queryHits(overlaps)] <- "INTRON"
  arm.gr$gene[queryHits(overlaps)] <- rna.gc.gr$gene_id[subjectHits(overlaps)]
  
  # Anything else that mapped to genome is intergenic
  arm.gr$biotype[as.character(seqnames(arm.gr)) %in% chromosomes & is.na(arm.gr$biotype)] <- "intergenic"
  arm.gr$region[as.character(seqnames(arm.gr)) %in% chromosomes & is.na(arm.gr$region)] <- "INTERGENIC"
  arm.gr$gene[as.character(seqnames(arm.gr)) %in% chromosomes & arm.gr$region == "INTERGENIC"] <- "intergenic"
  
  # 3. Annotate those that map to rRNA and tRNA
  arm.gr$biotype[grep("rRNA", as.character(seqnames(arm.gr)))] <- "rRNA"
  arm.gr$region[grep("rRNA", as.character(seqnames(arm.gr)))] <- "rRNA"
  arm.gr$gene[grep("rRNA", as.character(seqnames(arm.gr)))] <- as.character(seqnames(arm.gr))[grep("rRNA", as.character(seqnames(arm.gr)))]
  arm.gr$biotype[grep("tRNA", as.character(seqnames(arm.gr)))] <- "tRNA"
  arm.gr$region[grep("tRNA", as.character(seqnames(arm.gr)))] <- "tRNA"
  arm.gr$gene[grep("tRNA", as.character(seqnames(arm.gr)))] <- as.character(seqnames(arm.gr))[grep("tRNA", as.character(seqnames(arm.gr)))]
  
  # Check everything annotated
  if(all(sum(table(arm.gr$biotype)) == length(arm.gr), sum(table(arm.gr$region)) == length(arm.gr))) {
    return(arm.gr)
  } else {
    cat("Not all arms annotated.")
  }
  
}

# =============================================================================
# Split GRanges to dt
# =============================================================================

#' Convert GRanges to data.table
#'
#' Converts from GRanges to data.table. If only one arm then returns single data.table, if both arms then returns hybrid data.table
#'
#' @param gr GRanges to convert
#' @return data.table or hybrid data.table
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @export

ConvertToDataTable <- function(gr) {
  
  left.gr <- gr[gr$arm == "L"]
  right.gr <- gr[gr$arm == "R"]
  
  left.dt <- data.table(as.data.frame(left.gr))
  setnames(left.dt, names(left.dt), paste0("L_", names(left.dt)))
  setnames(left.dt, "L_read", "read")
  setkey(left.dt, read)
  
  right.dt <- data.table(as.data.frame(right.gr))
  setnames(right.dt, names(right.dt), paste0("R_", names(right.dt)))
  setnames(right.dt, "R_read", "read")
  setkey(right.dt, read)
  
  if(nrow(left.dt) != 0 & nrow(right.dt) != 0) {
    results.dt <- merge(left.dt, right.dt)
  } else if(nrow(left.dt) == 0) {
    results.dt <- right.dt
  } else if(nrow(right.dt) == 0) {
    results.dt <- left.dt
  }
  
  return(results.dt)
  
}

# =============================================================================
# Combined function to go from unique.hybrids data.table to annotated unique hybrids data.table
# =============================================================================

#' Annotate hybrid data.table
#'
#' Annotates each arm of a hybrid data.table using rna.tc.gr
#'
#' @param hybrid.dt hybrid data.table to annotate
#' @return annotated data.table
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @export

AnnotateHybrid <- function(hybrid.dt) {
  
  # First convert to GRanges using "all" to get all the reads as a GRanges
  hybrid.gr <- ConvertToGRanges(hybrid.dt, arm = "all")
  
  # Then annotate all the reads
  hybrid.gr <- AnnotateHybridArm(hybrid.gr)
  
  # Then split arms back to hybrid data.table format
  results.dt <- ConvertToDataTable(hybrid.gr)
  
  return(results.dt)
  
}

# =============================================================================
# Select genic reads (only for those with both arms mapped to genes)
# =============================================================================

#' Select genic reads
#'
#' Select genic reads (only for those with both arms mapped to genes)
#'
#' @param hybrid.dt hybrid data.table
#' @param type type of genic read
#' @param class class of genic read (i.e. include introns or not)
#' @return hybrid data.table
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @export

SelectGenicHybrid <- function(hybrid.dt, type = c("all", "intra", "inter"), class = c("mature", "immature")) {
  
  if(class == "immature") {
    genic <- hybrid.dt[intersect(grep("^ENSG", L_gene), grep("^ENSG", R_gene))] # select those read for with both arms map to a gene, use gene column to include introns (i.e. immature tx)
  } else if(class == "mature") {
    genic <- hybrid.dt[intersect(grep("^ENSG", L_seqnames), grep("^ENSG", R_seqnames))] # use seqnames column for transcriptome
    genic <- genic[L_region != "INTERGENIC"] # for those reads that map to transcriptome but on negative strand
    genic <- genic[R_region != "INTERGENIC"]
  }
  
  if(type == "all") {
    return(genic)
  } else if(type == "intra") {
    intra <- genic[L_gene == R_gene] # map to same gene
    return(intra)
  } else if(type == "inter") {
    inter <- genic[L_gene != R_gene] # map to different genes
    return(inter)
  }
}

# =============================================================================
# Find confident islands that have > 1 read
# =============================================================================

#' Find confident islands
#'
#' Find confident hybrid islands
#'
#' @param hybrid.dt hybrid data.table
#' @param minimum.coverage minimum coverage for an island
#' @param minimum.width minimum width for an island
#' @return annotated hybrid data.table
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @importFrom IRanges slice
#' @importFrom GenomeInfoDb seqlevels
#' @export


FindIslands <- function(hybrid.dt, minimum.coverage, minimum.width) {
  
  ### This only works (currently) for intragenic hybrids
  
  # Create GRanges of L reads
  L.gr <- ConvertToGRanges(hybrid.dt, arm = "left")
  L.coverage <- coverage(L.gr) # calculate coverage per seqname as RleList
  L.island <- slice(L.coverage, minimum.coverage) # apply minimum coverage threshold and returns RleViewsList
  
  # Create GRanges of L islands
  L.island.gr <- GRanges(seqnames = Rle(rep(names(L.island), sapply(L.island, length))),
                         ranges = IRanges(start = unlist(start(L.island)), end = unlist(end(L.island))),
                         strand = Rle(rep("+", length(unlist(start(L.island))))),
                         peakcoverage = unlist(viewMaxs(L.island)))
  L.island.gr <- L.island.gr[width(L.island.gr) > minimum.width] # apply minimum island width threshold
  L.island.gr$islandid <- paste0(1:length(L.island.gr), rep("+", length(L.island.gr))) # Need to specify strand for later merging
  
  # Now get the matching R confident regions for each L confident region as a list with each item corresponding to one L confident region
  # FOR LOOP method - returns a GRangeslist, easier to manipulate than LAPPLY method (cannot unlist)
  result.grl <- GRangesList(L = GRanges(), R = GRanges())
  R.gr <- ConvertToGRanges(hybrid.dt, arm = "right") # taking this out of the loop will probably speed things up...
  
  for(x in 1:length(L.island.gr)) {
    
    # Get R reads that match the L confident region
    matching.L.gr <- subsetByOverlaps(L.gr, L.island.gr[x])
    matching.R.gr <- R.gr[R.gr$read %in% matching.L.gr$read] # select those R reads that match
    seqlevels(matching.R.gr) <- unique(as.character(seqnames(matching.R.gr))) # reduce R seqlevels to 1 to speed up splice operation in loop
    # Create GRanges of matching R reads
    
    R.coverage <- coverage(matching.R.gr) # calculate coverage per seqname as RleList
    R.island <- slice(R.coverage, minimum.coverage) # apply minimum coverage threshold and returns RleViewsList
    
    # Only create GRanges of R island if there is an island, otherwise return an empty GRanges
    if(length(unlist(start(R.island))) > 0) {
      R.island.gr <- GRanges(seqnames = Rle(rep(names(R.island), sapply(R.island, length))),
                             ranges = IRanges(start = unlist(start(R.island)), end = unlist(end(R.island))),
                             strand = Rle(rep("+", length(unlist(start(R.island))))),
                             peakcoverage = unlist(viewMaxs(R.island)),
                             islandid = L.island.gr[x]$islandid)
      R.island.gr <- R.island.gr[width(R.island.gr) > minimum.width] # apply minimum width threshold
      if(length(R.island.gr) > 0) {
        R.island.gr <- R.island.gr[R.island.gr$peakcoverage == max(R.island.gr$peakcoverage)]
        R.island.gr <- sort(R.island.gr)[1]
      }
    } else {
      R.island.gr <- GRanges()
    }
    
    seqlevels(R.island.gr) <- seqlevels(L.island.gr)
    result.grl$L <- c(result.grl$L, L.island.gr[x])
    result.grl$R <- c(result.grl$R, R.island.gr)
  }
  
  # Now need to get back to data.table format
  L.islands.gr <- result.grl$L
  names(L.islands.gr) <- NULL
  L.islands.gr <- AnnotateHybridArm(L.islands.gr) # annotate islands
  L.islands.dt <- data.table(as.data.frame(L.islands.gr))
  L.islands.dt[, arm := "L"]
  setkey(L.islands.dt, islandid)
  
  R.islands.gr <-result.grl$R
  names(R.islands.gr) <- NULL
  R.islands.gr <- AnnotateHybridArm(R.islands.gr)
  R.islands.dt <- data.table(as.data.frame(R.islands.gr))
  R.islands.dt[, arm := "R"]
  setkey(R.islands.dt, islandid)
  
  islands.dt <- merge(L.islands.dt, R.islands.dt)
  setnames(islands.dt,
           c("islandid", "seqnames.x", "start.x", "end.x", "width.x", "strand.x", "peakcoverage.x", "biotype.x", "region.x", "gene.x", "arm.x",
             "seqnames.y", "start.y", "end.y", "width.y", "strand.y", "peakcoverage.y", "biotype.y", "region.y", "gene.y", "arm.y"),
           c("islandid", "L_seqnames", "L_start", "L_end", "L_width", "L_strand", "L_peakcoverage", "L_biotype", "L_region", "L_gene", "L_arm",
             "R_seqnames", "R_start", "R_end", "R_width", "R_strand", "R_peakcoverage", "R_biotype", "R_region", "R_gene", "R_arm"))
  
  return(islands.dt)
}

# =============================================================================
# Find hybrid reads that do not overlap an island
# =============================================================================

#' Find non-island hybrid reads
#'
#' Find hybrid reads that do not overlap an island
#'
#' @param reads hybrid data.table of reads
#' @param islands hybrid data.table of islands
#' @return hybrid data.table
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @export

FindNonIslandHybrids <- function(reads, islands) {
  
  arms <- c("left", "right") # create a list of L reads that overlap L islands and R reads that overlap R islands (as L and R have previously been sorted so L is 5' of R shouldn't miss alternate hybrids where 1 arm overlaps an island)
  overlapping.reads.grl <- lapply(arms, function(x) {
    
    # First convert data.tables to GRanges
    reads.gr <- ConvertToGRanges(reads, arm = x)
    islands.gr <- ConvertToGRanges(islands, arm = x)
    
    # Then find reads overlapping the islands
    overlaps <- findOverlaps(reads.gr, islands.gr)
    overlapping.reads.gr <- reads.gr[queryHits(overlaps)]
    return(overlapping.reads.gr)
  })
  
  
  intersecting.reads <- intersect(overlapping.reads.grl[[1]]$read, overlapping.reads.grl[[2]]$read) # gets read names for those reads that overlap both L and R islands (as L and R have previously been sorted so L is 5' of R shouldn't miss alternate hybrids where only 1 arm overlaps an island)
  nonisland.reads <- reads[!(read %in% intersecting.reads)] # filter out the overlapping reads from all reads
  nonisland.reads[, L_peakcoverage := 1]
  nonisland.reads[, R_peakcoverage := 1]
  
  return(nonisland.reads)
}

# =============================================================================
# Join island and non-island hybrids to create non-redundant set of hybrids
# =============================================================================

#' Create non-redundant set of hybrids
#'
#' Join island and non-island hybrids to create non-redundant set of hybrids
#'
#' @param islands data.table of islands
#' @param nonislands data.table of non-island hybrid reads
#' @return annotated data.table
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @export

CreateNonRedundantHybrids <- function(islands, nonislands) {
  
  # Remove conflicted/extraneous columns
  setnames(islands, "islandid", "id")
  setnames(nonislands, "read", "id")
  nonislands[, c("L_qname", "R_qname") := NULL]
  
  # Merge
  nonredundant <- list(islands, nonislands)
  nonredundant <- rbindlist(nonredundant, use.names = TRUE, fill = TRUE)
  
  return(nonredundant)
  
}

# =============================================================================
# Convert from transcriptomic to genomic coordinates
# =============================================================================

#' Convert to genomic coordinates
#'
#' Convert from transcriptomic to genomic coordinates
#'
#' @param hybrid.dt hybrid data.table
#' @return GRangesList of length two (L and R arms)
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import GenomicRanges
#' @export

# Based on F. Agostini's convertCoordinates function

ConvertToGenomicCoordinates <- function(hybrid.dt) {
  # For this (esp. in terms of negative strand genes, "start" in genomic coordinates means the left-most coordinate, i.e. what would return if did start(gr) [this is actually the end as we would understand it for negative genes])
  
  results.grl <- GRangesList(L = GRanges(), R = GRanges())
  
  for(i in 1:nrow(hybrid.dt)) {
    
    # Get details to add to GRanges at end - read/id and coverage/score
    if("read" %in% colnames(hybrid.dt)) {
      id <- hybrid.dt[i]$read # if reads
    } else {
      id <- hybrid.dt[i]$id # if after islands/duplexes
    }
    
    if("L_peakcoverage" %in% colnames(hybrid.dt)) {
      score <- min(hybrid.dt[i]$L_peakcoverage, hybrid.dt[i]$R_peakcoverage)
    } else {
      score <- 1
    }
    
    # First do left arm
    # Seqnames are ENSG, get corresponding ENST
    tx_id <- longest.rna.dt$ensembl_transcript_id[match(hybrid.dt[i]$L_seqnames, longest.rna.dt$ensembl_gene_id)]
    tx_length <- longest.rna.dt$tx_len[match(hybrid.dt[i]$L_seqnames, longest.rna.dt$ensembl_gene_id)]
    
    # Get corresponding exons.gr for that tx
    exons.gr <- unlist(exons_tx.grl[names(exons_tx.grl) == tx_id])
    exons.gr$cumwidths <- cumsum(width(exons.gr))
    
    # Define start and end (and switch if tx/gene is on negative strand)
    if(unique(as.character(strand(exons.gr))) == "+") {
      start <- hybrid.dt[i]$L_start
      end <- hybrid.dt[i]$L_end
    } else if(unique(as.character(strand(exons.gr))) == "-") {
      start <- hybrid.dt[i]$L_end
      end <- hybrid.dt[i]$L_start
    }
    
    # Find the exon in which the start position is located
    startexon <- length(which(exons.gr$cumwidths - start < 0)) + 1
    exon.gr <- exons.gr[startexon] # first identifies the exons before the start and then +1 to get the exon the start is located in, because using GRanges accounts for strand in order
    
    # Calculate genomic coordinates
    if(unique(as.character(strand(exons.gr)) == "+")) {
      if(exon.gr$exon_rank == 1) {
        start.g <- start + start(exon.gr) - 1
      } else if(exon.gr$exon_rank > 1) {
        start.g <- start + start(exon.gr) - exons.gr$cumwidths[startexon - 1] - 1 # this is the left most coord, - 1 for 1-based, - the widths of the preceding exons
      }
    } else if(unique(as.character(strand(exons.gr)) == "-")) {
      start.g <- start(exon.gr) + (exons.gr$cumwidths[startexon] - start) # no need to +/- 1
    }
    
    # Find the exon in which the end position is located
    endexon <- length(which(exons.gr$cumwidths - end < 0)) + 1
    exon.gr <- exons.gr[endexon] # first identifies the exons before the start and then +1 to get the exon the start is located in, because using GRanges accounts for strand in order
    
    # Calculate genomic coordinates
    if(unique(as.character(strand(exons.gr)) == "+")) {
      if(exon.gr$exon_rank == 1) {
        end.g <- end + start(exon.gr) - 1
      } else if(exon.gr$exon_rank > 1) {
        end.g <- end + start(exon.gr) - exons.gr$cumwidths[endexon - 1] - 1 # this is the left most coord, - 1 for 1-based, - the widths of the preceding exons
      }
    } else if(unique(as.character(strand(exons.gr)) == "-")) {
      end.g <- start(exon.gr) + (exons.gr$cumwidths[endexon] - end) # no need to +/- 1
    }
    
    # Account for any reads that cross a splice junction
    if(endexon == startexon) {
      L_result.gr <- GRanges(seqnames = seqnames(exon.gr),
                             ranges = IRanges(start = start.g, width = hybrid.dt[i]$L_width),
                             strand = strand(exon.gr),
                             id = id,
                             score = score)
      results.grl$L <- c(results.grl$L, L_result.gr)
    } else {
      # first do start
      L_start.gr <- GRanges(seqnames = seqnames(exon.gr),
                            ranges = IRanges(start = start.g, end = end(exons.gr[startexon])),
                            strand = strand(exon.gr),
                            id = id,
                            score = score)
      results.grl$L <- c(results.grl$L, L_start.gr)
      # then do end
      L_end.gr <- GRanges(seqnames = seqnames(exon.gr),
                          ranges = IRanges(start = start(exons.gr[endexon]), end = end.g),
                          strand = strand(exon.gr),
                          id = id,
                          score = score)
      results.grl$L <- c(results.grl$L, L_end.gr)
    }
    
    # Next do right arm
    # Seqnames are ENSG, get corresponding ENST
    tx_id <- longest.rna.dt$ensembl_transcript_id[match(hybrid.dt[i]$L_seqnames, longest.rna.dt$ensembl_gene_id)]
    tx_length <- longest.rna.dt$tx_len[match(hybrid.dt[i]$L_seqnames, longest.rna.dt$ensembl_gene_id)]
    
    # Get corresponding exons.gr for that tx
    exons.gr <- unlist(exons_tx.grl[names(exons_tx.grl) == tx_id])
    exons.gr$cumwidths <- cumsum(width(exons.gr))
    
    # Define start and end (and switch if tx/gene is on negative strand)
    if(unique(as.character(strand(exons.gr))) == "+") {
      start <- hybrid.dt[i]$R_start
      end <- hybrid.dt[i]$R_end
    } else if(unique(as.character(strand(exons.gr))) == "-") {
      start <- hybrid.dt[i]$R_end
      end <- hybrid.dt[i]$R_start
    }
    
    # Find the exon in which the start position is located
    startexon <- length(which(exons.gr$cumwidths - start < 0)) + 1
    exon.gr <- exons.gr[startexon] # first identifies the exons before the start and then +1 to get the exon the start is located in, because using GRanges accounts for strand in order
    
    # Calculate genomic coordinates
    if(unique(as.character(strand(exons.gr)) == "+")) {
      if(exon.gr$exon_rank == 1) {
        start.g <- start + start(exon.gr) - 1
      } else if(exon.gr$exon_rank > 1) {
        start.g <- start + start(exon.gr) - exons.gr$cumwidths[startexon - 1] - 1 # this is the left most coord, - 1 for 1-based, - the widths of the preceding exons
      }
    } else if(unique(as.character(strand(exons.gr)) == "-")) {
      start.g <- start(exon.gr) + (exons.gr$cumwidths[startexon] - start) # no need to +/- 1
    }
    
    # Find the exon in which the end position is located
    endexon <- length(which(exons.gr$cumwidths - end < 0)) + 1
    exon.gr <- exons.gr[endexon] # first identifies the exons before the start and then +1 to get the exon the start is located in, because using GRanges accounts for strand in order
    
    # Calculate genomic coordinates
    if(unique(as.character(strand(exons.gr)) == "+")) {
      if(exon.gr$exon_rank == 1) {
        end.g <- end + start(exon.gr) - 1
      } else if(exon.gr$exon_rank > 1) {
        end.g <- end + start(exon.gr) - exons.gr$cumwidths[endexon - 1] - 1 # this is the left most coord, - 1 for 1-based, - the widths of the preceding exons
      }
    } else if(unique(as.character(strand(exons.gr)) == "-")) {
      end.g <- start(exon.gr) + (exons.gr$cumwidths[endexon] - end) # no need to +/- 1
    }
    
    # Account for reads that cross a splice junction
    if(endexon == startexon) {
      R_result.gr <- GRanges(seqnames = seqnames(exon.gr),
                             ranges = IRanges(start = start.g, width = hybrid.dt[i]$R_width),
                             strand = strand(exon.gr),
                             id = id,
                             score = score)
      results.grl$R <- c(results.grl$R, R_result.gr)
    } else {
      # first do start
      R_start.gr <- GRanges(seqnames = seqnames(exon.gr),
                            ranges = IRanges(start = start.g, end = end(exons.gr[startexon])),
                            strand = strand(exon.gr),
                            id = id,
                            score = score)
      results.grl$R <- c(results.grl$R, R_start.gr)
      # then do end
      R_end.gr <- GRanges(seqnames = seqnames(exon.gr),
                          ranges = IRanges(start = start(exons.gr[endexon]), end = end.g),
                          strand = strand(exon.gr),
                          id = id,
                          score = score)
      results.grl$R <- c(results.grl$R, R_end.gr)
    }
  }
  
  return(results.grl)
  
}


# =============================================================================
# Function to export genomic coordinate BED file from converted coords output
# =============================================================================

#' Create BED file
#'
#' Export genomic coordinate BED file from converted coords output
#'
#' @param gc.grl GRangesList of hybrids in genomic coordinates
#' @param filename Filename to save
#' @return BED file
#'
#' @import GenomicRanges
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @importFrom rtracklayer export.bed
#' @export

CreateBED <- function(gc.grl, filename) {
  
  results.dir <- "results"
  
  seqlevelsStyle(gc.grl) <- "UCSC"
  L <- gc.grl$L
  # L$id <- 1:length(gc.grl$L) # give each hybrid an ID
  R <- gc.grl$R
  # R$id <- 1:length(gc.grl$R)
  combined <- c(L, R)
  combined.grl <- split(combined, combined$id) # split into grl where each entry is a GRanges of 2 - L and R arm
  export.bed(combined.grl, con = file.path(results.dir, filename))
  
  # Add score
  score <- data.table(id = combined$id, score = combined$score)
  setkey(score, id)
  score <- unique(score, by = "id") # need to get unique score by id
  
  bed <- read.delim(file.path(results.dir, filename), header = FALSE)
  bed[, 5] <- score$score
  write.table(bed, file = file.path(results.dir, filename), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}

# =============================================================================
# Wrapper function to export genomic coordinate BED file from hybrid data.table
# =============================================================================

#' Export BED file
#'
#' Export genomic coordinate BED file from hybrid data.table
#'
#' @param hybrid.dt hybrid data.table
#' @param filename Filename to save
#' @return BED file
#'
#' @export

ExportBED <- function(hybrid.dt, filename) {
  
  hybrid.gc <- ConvertToGenomicCoordinates(hybrid.dt)
  CreateBED(hybrid.gc, filename = filename)
  
  print("BED file exported")
  
}