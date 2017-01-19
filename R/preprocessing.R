# =============================================================================
# Create directory structure
# =============================================================================

#' Create directory structure
#'
#' Create directory structure
#'
#' @param folder Folder in which to create directory structure
#' @return New directories
#'
#' @export

CreateDirectories <- function(folder) {
  
  print(paste("The directory structure will be created in", folder))
  
  if(dir.exists(folder) == FALSE) {
    dir.create(folder)
  }
  
  setwd(folder)
  
  # Create directories if they do not exist
  directories <- c("ref",
                   "data",
                   "results",
                   "plots",
                   "doc")
  
  for(dir in directories) {
    if(dir.exists(dir) == FALSE) {
      dir.create(dir)
    }
  }
  
  # Create subdirectories if they do not exist
  subdirectories <- c("results/processed_data",
                      "results/mapped",
                      "results/rnahybrid",
                      "results/mapped/temp",
                      "results/processed_data/temp",
                      "ref/bowtie")
  
  for(subdir in subdirectories) {
    if(dir.exists(subdir) == FALSE) {
      dir.create(subdir)
    }
  }
  
  print(paste("Directory structure created."))
}

# =============================================================================
# Create Bowtie indices
# =============================================================================

#' Create Bowtie indices
#'
#' Create Bowtie indices
#'
#' @return Bowtie indices
#' @importFrom rentrez entrez_fetch
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom Biostrings DNAStringSet getSeq
#' @importFrom ShortRead writeFasta readFasta
#'
#' @export

CreateBowtieIndices <- function() {
  
  ref.dir <- file.path("ref/bowtie")
  
  # =============================================================================
  # rRNA from NCBI
  # =============================================================================
  
  if(!file.exists(file.path(ref.dir, "rRNA.fa"))) {
    print("Getting rRNA FASTA")
    rRNA.ids <- c("NR_023363.1","NR_003285.2","NR_003287.2","NR_003286.2") # NCBI rRNA 5S, 5.8S, 28S, 18S
    rRNA.fa <- entrez_fetch(db="nucleotide", id=rRNA.ids, rettype="fasta")
    write(rRNA.fa, file=file.path(ref.dir, "rRNA.fa"))
  }
  
  # =============================================================================
  # tRNA from GtRNAdb
  # =============================================================================
  
  if(!file.exists(file.path(ref.dir, "rRNA_tRNA.fa"))) {
    print("Getting tRNA FASTA")
    tRNA.fa <- file.path(ref.dir, "tRNA.fa")
    download.file("http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.fa", destfile = file.path(ref.dir, "tRNA.fa"))
    hs.tRNA <- readFasta(file.path(ref.dir, "tRNA.fa"))
    hs.tRNA.CCA <- DNAStringSet(paste(sread(hs.tRNA), "CCA", sep = ""))
    names(hs.tRNA.CCA) <- id(hs.tRNA)
    writeFasta(hs.tRNA.CCA, file.path(ref.dir, "tRNA.fa"))
    cmd <- paste0("cat ", file.path(ref.dir, "rRNA.fa"), " ", file.path(ref.dir, "tRNA.fa"), " > ", file.path(ref.dir, "rRNA_tRNA.fa")) # combine rRNA and tRNA files
    system(cmd)
  }
  
  # =============================================================================
  # mtRNA from NCBI (rCRS version 2013) and pre-rRNA
  # =============================================================================
  
  if(!file.exists(file.path(ref.dir, "mtDNA_prerRNA.fa"))) {
    print("Getting mtDNA and pre-rRNA FASTA")
    mtDNA.prerRNA.fa <- entrez_fetch(db = "nucleotide", id = c("NC_012920", "NR_046235.1"), rettype = "fasta") # mtDNA and pre-rRNA 45S
    write(mtDNA.prerRNA.fa, file = file.path(ref.dir, "mtDNA_prerRNA.fa"))
  }
  
  # =============================================================================
  # longest RNA (pc and nc)
  # =============================================================================
  
  if(!file.exists(file.path(ref.dir, "longestRNA.fa"))) {
    print("Getting longest RNA FASTA")
    # Get sequences for these RNA transcripts (i.e. introns)
    seqlevelsStyle(exons_tx.grl) <- "UCSC" # switch style for BSgenome
    exons_tx.seq <- getSeq(Hsapiens, exons_tx.grl, as.character = TRUE) # gets sequences for each exon in each transcript
    
    # Write FASTA
    exons_tx.seq <- sapply(exons_tx.seq, function(x) paste(x, collapse = "")) # collapse the exons for each transcript down into one character vector
    exons_tx.ds <- DNAStringSet(exons_tx.seq)
    names(exons_tx.ds) <- longest.rna.dt$ensembl_gene_id[match(names(exons_tx.ds), longest.rna.dt$ensembl_transcript_id)] # Add matching gene_id as seqname
    writeFasta(exons_tx.ds, file = file.path(ref.dir, "longestRNA.fa"))
  }
  
  # =============================================================================
  # Genome from Ensembl (GRCh37.75 to match TxDb)
  # =============================================================================
  
  if(!file.exists(file.path(ref.dir, "GRCh37.75.fa"))) {
    print("Getting genome primary assembly FASTA")
    download.file("ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz", destfile = file.path(ref.dir, "GRCh37.75.fa.gz"))
    cmd <- paste0("gunzip ", file.path(ref.dir, "GRCh37.75.fa.gz")) # unzip for bowtie2-build
    system(cmd)
  }
  
  # =============================================================================
  # Build Bowtie indices
  # =============================================================================
  
  print("Building Bowtie indices")
  indices <- c(file.path(ref.dir, "rRNA_tRNA.fa"), file.path(ref.dir, "mtDNA_prerRNA.fa"), file.path(ref.dir, "longestRNA.fa"), file.path(ref.dir, "GRCh37.75.fa"))
  for(index in indices) {
    index.name <- substr(index, 1, regexpr(".fa", index, fixed = TRUE) -1)
    cmd <- paste0("bowtie-build ", index, " ", index.name)
    system(cmd)
  }
  
}