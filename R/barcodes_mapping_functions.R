# =============================================================================
# Extract random barcodes and leave experiment barcode at front
# =============================================================================

#' Extract random barcodes and leave experiment barcode at front
#'
#' Extract random barcodes and leave experiment barcode at front
#'
#' @param fq.file Path to FASTQ file
#' @param barcode.filename File for barcode index table
#' @return FASTQ file and random barcode tsv
#'
#' @importFrom ShortRead readFastq writeFastq ShortReadQ sread id
#' @importFrom Biostrings BStringSet quality narrow xscat BString DNAStringSet DNAString
#' @export

ExtractRandomBarcodes <- function(fq.file, barcode.filename) {
  
  data.dir <- "data/"
  processed.dir <- "results/processed_data/"
  
  # Step 1 - Adds read number to header for later identification (have kept quality data, unlike original paper)
  cat("Assigning read number to FASTQ header\n")
  fq <- readFastq(paste0(data.dir, fq.file))
 
  # Ensure reads are long enough, otherwise narrow won't work
  if(min(width(fq)) < 10) {
  print("Some reads are shorter than random + experimental barcode. Aborting")
  stop()
  }
 
  id <- BStringSet(1:length(fq)) # renames by read number
  fq <- ShortReadQ(sread = sread(fq), id = id, quality = quality(fq))
  
  # Step 2 - reconfigure barcode layout (replicates swapBarcode.py)
  cat("Extracting random barcode...")
  
  # 123456789
  # RRRXXXXRR
  # to
  # 4567NNNN	12389
  # XXXXNNNN	RRRRR
  
  # Extract reads and qualities from ShortRead object
  reads.ds <- sread(fq)
  quality.bs <- quality(quality(fq)) # to get to BString
  
  # Extract 1st and 2nd random barcodes and concatenate
  r1.barcode.ds <- narrow(reads.ds, start = 1, width = 3)
  r1.quality.bs <- narrow(quality.bs, start = 1, width = 3)
  r2.barcode.ds <- narrow(reads.ds, start = 8, width = 2)
  r2.quality.bs <- narrow(quality.bs, start = 8, width = 2)
  r.barcode.ds <- xscat(r1.barcode.ds, r2.barcode.ds)
  r.quality.bs <- xscat(r1.quality.bs, r2.quality.bs)
  
  # Create data frame of random barcodes
  r.barcode <- data.frame(read = as.character(id), random_barcode = as.character(r.barcode.ds))
  write.table(r.barcode, file = paste0(processed.dir, barcode.filename), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  # Extract experimental barcode
  x.barcode.ds <- narrow(reads.ds, start = 4, width = 4)
  x.quality.bs <- narrow(quality.bs, start = 4, width = 4)
  
  # Remove barcodes and add back experimental barcode only
  reads.ds <- narrow(reads.ds, start = 10) # Trim 1st random barcode and experimental barcode
  reads.ds <- xscat(x.barcode.ds, reads.ds) # Then add experimental barcode back
  
  # Do same for corresponding quality scores
  quality.bs <- narrow(quality.bs, start = 10) # Trim 1st random barcode and experimental barcode
  quality.bs <- xscat(x.quality.bs, quality.bs)
  
  # Recreate ShortReadQ and write reordered FASTQ
  fq <- ShortReadQ(sread = reads.ds, id = id(fq), quality = quality.bs)
  filename <- gsub(".fq.gz", "_reordered.fq.gz", fq.file)
  writeFastq(fq, file = paste0(processed.dir, filename)) # This is very slow... might be an idea to only write once
  # return(r.barcode)
  
}

# =============================================================================
# Demultiplex reads using experimental barcode
# =============================================================================

#' Demultiplex reads using experimental barcode
#'
#' Demultiplex reads using experimental barcode
#'
#' @param fq.file Path to FASTQ file
#' @param ... list of experiments and barcodes e.g. experiment_1 = "NNNN"
#' @return FASTQ files as per experimental barcodes
#'
#' @export

DemultiplexReads <- function(fq.file, ...) {
  
  processed.dir <- "results/processed_data/"
  
  # Get experimental barcodes into correct format
  args <- list(...)
  exp <- character()
  for(n in 1:length(args)){
    exp[n] <- paste0(names(args)[n], "=^", args[[n]])
  }
  exp <- paste0(exp, collapse = " -g ")
  
  # Cutadapt command
  cmd <- paste0("cutadapt -e 0 --no-indels -g ", exp, " -o ", processed.dir, "{name}.fq.gz ", processed.dir, fq.file)
  print(cmd)
  system(cmd)
  
  # zcat LUs27_196_reordered.fq.gz | fastx_barcode_splitter.pl --bcfile LUs27_196.barcode --bol --exact --prefix fastx_ --suffix .fq
  # cutadapt is slower than fastx_toolkit, but gets same results. I have kept with it to reduce number of dependencies
  
}

# =============================================================================
# Remove adapters to select hybrid reads (version 2 - more like YS)
# =============================================================================

#' Remove adapters to select hybrid reads
#'
#' Remove A and B adapters to select hybrid reads expanding on YS's strategy
#'
#' @param fq.file Path to FASTQ file
#' @param adapterA Adapter A sequence
#' @param adapterB Adapter B sequence
#'
#' @return FASTQ file of L and R arms (unzipped)
#'
#' @importFrom ShortRead readFastq writeFastq ShortReadQ sread id
#' @importFrom Biostrings BStringSet quality DNAStringSet BString DNAString vmatchPattern vcountPattern
#' @importFrom BiocGenerics append
#' @importFrom IRanges IRanges IRangesList
#' @export

ExtractHybrids <- function(fq.file, adapterA, adapterB) {
  
  # Set filename base for outputs
  processed.dir <- "results/processed_data/"
  temp.dir <- "results/processed_data/temp/"
  filename <- basename(substr(fq.file, 1, regexpr(".fq.gz", fq.file, fixed = TRUE) -1))
  
  # First remove B-A at 3' end (...BA)
  cat("Removing BA from ...BA\n")
  input <- paste0(processed.dir, fq.file)
  trimmed.output <- paste0(temp.dir, filename, "_BA.fq.gz")
  untrimmed.output <- paste0(temp.dir, filename, "_noBA.fq.gz")
  cmd <- paste0("cutadapt -a ", adapterB, adapterA, " -e 0.06 -n 10 -o ", trimmed.output, " --untrimmed-output ", untrimmed.output, " ", input)
  system(cmd)
  
  # Then remove A from 3' end of those without B-A (...A, not ...BA)
  cat("Removing A from ...A\n")
  input <- paste0(temp.dir, filename, "_noBA.fq.gz")
  trimmed.output <- paste0(temp.dir, filename, "_A_noBA.fq.gz")
  untrimmed.output <- paste0(temp.dir, filename, "_noA_noBA.fq.gz")
  cmd <- paste0("cutadapt -a ", adapterA, " -e 0.06 -n 10 -O 10 -o ", trimmed.output, " --untrimmed-output ", untrimmed.output, " ", input)
  system(cmd)
  
  # Then select 1. B and A, but not BA to identify (...B...A) and 2. B and BA to identify ...B...BA - input will be ...B...; should not be possible to get ...BA...A, or ...A...A
  cat("Removing B from ...B...(A/BA)\n")
  
  # Loads both FASTQ files, combine and select only those with an adapter B for further processing
  fq1 <- readFastq(paste0(temp.dir, filename, "_A_noBA.fq.gz"))
  fq2 <- readFastq(paste0(temp.dir, filename, "_BA.fq.gz"))
  fq <- append(fq1, fq2) # combines the two
  B <- vcountPattern(adapterB, subject = sread(fq), max.mismatch = 2, with.indels = FALSE) # counts instance of B
  fq <- fq[B > 0] # selects only those with 1 (or more) B
  
  # Extract components of ShortRead
  reads.ds <- sread(fq)
  quality.bs <- quality(fq)@quality # extract quality as a BString
  names(reads.ds) <- id(fq)
  names(quality.bs) <- id(fq)
  
  # Find locations of adapter B and reduce instances of BB/BBB
  hybrid.ir <- unlist(vmatchPattern(adapterB, subject = reads.ds, max.mismatch = 2, with.indels = FALSE)) # can unlist and still know which read it refers to as has name
  hybrid.irl <- split(hybrid.ir, names(hybrid.ir)) # split into list by read
  hybrid.ir <- unlist(reduce(hybrid.irl)) # reduce polyB and unlist back
  
  # Remove those that are ...B...B...
  duplicated <- names(hybrid.ir)[duplicated(names(hybrid.ir))] # identify those that are ...B...B... and remove from ongoing analysis (may remove a few hybrids that are ...B...B, but these are very few ~150)
  hybrid.ir <- hybrid.ir[!(names(hybrid.ir) %in% duplicated), ] # so these are the ranges of B for reads that are ...B...A or ...B...BA
  reads.ds <- reads.ds[names(reads.ds) %in% names(hybrid.ir), ] # and these are the reads for the above hybrids
  quality.bs <- quality.bs[names(quality.bs) %in% names(hybrid.ir)] # and there are the qualities
  
  # Order reads
  hybrid.ir <- hybrid.ir[order(names(hybrid.ir))] # order both by read name for trimming step
  reads.ds <- reads.ds[order(names(reads.ds))]
  quality.bs <- quality.bs[order(names(quality.bs))]
  
  #Remove reads that are too short for trimming
  tooshort <- end(hybrid.ir) > width(reads.ds) # need to remove reads that are too short
  hybrid.ir <- hybrid.ir[!tooshort]
  reads.ds <- reads.ds[!(tooshort)]
  quality.bs <- quality.bs[!(tooshort)]
  
  # Remove reads that have a negative start (i.e. adapter B overhang)
  overhang <- start(hybrid.ir) < 1
  hybrid.ir <- hybrid.ir[!overhang]
  reads.ds <- reads.ds[!(overhang)]
  quality.bs <- quality.bs[!(overhang)]

  # Create R and L arms of hybrid reads
  cat("Splitting right and left arms\n")
  R.reads.ds <- DNAStringSet(reads.ds, start = 1, end = (start(hybrid.ir)-1)) # R arm is from start of read to start of adapter B
  R.quality.bs <- BStringSet(quality.bs, start = 1, end = (start(hybrid.ir)-1))
  L.reads.ds <- DNAStringSet(reads.ds, start = (end(hybrid.ir)+1)) # L arm is from end of adapter B to end of read
  L.quality.bs <- BStringSet(quality.bs, start = (end(hybrid.ir)+1))
  
  # Only keep those that have both R and L arms > 16 nt
  R.reads.ds <- R.reads.ds[width(R.reads.ds) > 16]
  L.reads.ds <- L.reads.ds[width(L.reads.ds) > 16]
  R.reads.ds <- R.reads.ds[names(R.reads.ds) %in% names(L.reads.ds), ]
  L.reads.ds <- L.reads.ds[names(L.reads.ds) %in% names(R.reads.ds), ]
  
  # Get corresponding qualities
  R.quality.bs <- R.quality.bs[names(R.quality.bs) %in% names(R.reads.ds)]
  L.quality.bs <- L.quality.bs[names(L.quality.bs) %in% names(L.reads.ds)]
  
  # Create ShortReads
  R.hybrid <- ShortReadQ(sread = R.reads.ds[order(names(R.reads.ds))],
                         id = BStringSet(paste0(names(R.reads.ds)[order(names(R.reads.ds))], "_R")),
                         quality = R.quality.bs[order(names(R.reads.ds))])
  L.hybrid <- ShortReadQ(sread = L.reads.ds[order(names(L.reads.ds))],
                         id = BStringSet(paste0(names(L.reads.ds)[order(names(L.reads.ds))], "_L")),
                         quality = L.quality.bs[order(names(L.reads.ds))])
  hybrid <- append(R.hybrid, L.hybrid)
  
  # Write out FASTQ files
  cat("Writing FASTQ files\n")
  writeFastq(hybrid, compress = FALSE, file = paste0(processed.dir, filename, "_hybrid.fq")) # Do not compress as aligning using Bowtie
  
}

# =============================================================================
# Align and merge with Bowtie == Step 6 & 7 of runHybridAnalysis_transcripts-2.0.0.py & mergeSamHybrid.R
# =============================================================================

#' Align and merge reads using Bowtie (as in original paper, except no RC for step 3)
#'
#' Align reads and merge as directed
#'
#' @param fq.file FASTQ file (unzipped)
#' @param bam.file output BAM filename
#' @param processors Number of processors to use
#'
#' @return SAM file of alignments
#'
#' @export

# Also removes unmapped reads as in Step 7 (-F 4 in samtools view)

AlignReadsBowtie <- function(fq.file, bam.file, processors) {
  
  # Set paths to directories
  ref.dir <- file.path("ref/bowtie")
  mapped.dir <- file.path("results/mapped")
  temp.dir <- file.path("results/mapped/temp")
  processed.dir <- file.path("results/processed_data")
  
  # First map hybrid reads to rRNAs and tRNAs
  # Allow multiple hits and select best
  cmd <- paste0("bowtie -p ", processors, " -v 2 -k 1 --best --sam --un ", temp.dir, "/unmapped1.fq ", ref.dir, "/rRNA_tRNA ", file.path(processed.dir, fq.file), " | samtools view -hu -F 4 - | sambamba sort -t ", processors, " -o ", temp.dir, "/mapped1.bam /dev/stdin")
  print(paste0("Running ", cmd))
  system(cmd)
  
  # Second map those hybrid reads not mapped, rRNAs or tRNAs to mtDNA and pre-rRNA
  # Allow multiple hits and select best
  cmd <- paste0("bowtie -p ", processors, " -v 2 -k 1 --best --sam --un ", temp.dir, "/unmapped2.fq ", ref.dir, "/mtDNA_prerRNA ", temp.dir, "/unmapped1.fq | samtools view -hu -F 4 - | sambamba sort -t ", processors, " -o ", temp.dir, "/mapped2.bam /dev/stdin")
  print(paste0("Running ", cmd))
  system(cmd)
  
  # Third map those hybrid reads not mapped to pre-rRNAs, rRNAs and tRNAs or mtDNA to transcripts, do not allow reverse complement alignment for this setp
  # Allow multiple hits and select best - filtered out on loading for barcode analysis
  cmd <- paste0("bowtie -p ", processors, " -v 2 -m 1 --best --strata --sam --norc --un ", temp.dir, "/unmapped3.fq ", ref.dir, "/longestRNA ", temp.dir, "/unmapped2.fq | samtools view -hu -F 4 - | sambamba sort -t ", processors, " -o ", temp.dir, "/mapped3.bam /dev/stdin")
  print(paste0("Running ", cmd))
  system(cmd)
  
  # Fourth and finally map those hybrid reads not mapped to pre-rRNAs, rRNAs and tRNAs, mtDNA or transcripts to genome
  # If wanted to, could align unmapped2.fq to both transcripts and genome, then load in with MAPQ filter of 40 and only select those reads mapped to genome that were not mapped to transcriptome
  # However, multimappers to transcriptome would also be multimappers to genome (some might map unique as bowtie2 is not splice aware, but these should be excluded as actually they are multimappers really)
  cmd <- paste0("bowtie -p ", processors, " -v 2 -m 1 --best --strata --sam --un ", temp.dir, "/unmapped4.fq ", ref.dir, "/GRCh37.75 ", temp.dir, "/unmapped3.fq | samtools view -hu -F 4 - | sambamba sort -t ", processors, " -o ", temp.dir, "/mapped4.bam /dev/stdin")
  print(paste0("Running ", cmd))
  system(cmd)
  
  # Then merge mapped reads
  cmd <- paste0("sambamba merge -t ", processors, " ", mapped.dir, "/", bam.file, " ", temp.dir, "/mapped1.bam ", temp.dir, "/mapped3.bam ", temp.dir, "/mapped4.bam")
  print(paste0("Running ", cmd))
  system(cmd)
  
}