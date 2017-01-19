# =============================================================================
# Plot hybrid source
# =============================================================================

#' Plot hybrid source
#'
#' Plot the hybrid region for all arms
#'
#' @param hybrid.dt hybrid data.table
#' @param filename filename for plot with extension
#' @return plot and table to command line
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import ggplot2
#' @export

# Requires data.table, ggplot2, scales, RColorBrewer

PlotHybridRegion <- function(hybrid.dt, filename) {
  
  plot.dir <- "plots"
  
  # Extract data
  data <- data.frame(Region = c(hybrid.dt$L_region, hybrid.dt$R_region))
  data$Region <- factor(data$Region, levels = c("UTR5", "CDS", "UTR3", "INTRON", "rRNA", "tRNA", "lincRNA", "miRNA", "other", "INTERGENIC"))
  
  # Print summary table
  cat("Hybrid region [number]:")
  print(table(data$Region))
  cat("Hybrid region [%]:")
  print(round(prop.table(table(data$Region)) * 100, digits = 0))
  
  # Plot
  p <- ggplot(data, aes(x = factor(1), fill = Region)) +
    geom_bar(width = 1) +
    coord_polar(theta = "y") +
    labs(title = "Hybrid read region",
         x = "",
         y = "") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    scale_fill_manual(values = c("#eff3ff", "#bdd7e7", "#6baed6", "#2171b5", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#969696"))
  ggsave(filename = file.path(plot.dir, filename), plot = p)
  return(p)
}

# =============================================================================
# Plot genic region
# =============================================================================

#' Plot genic region
#'
#' Plot the hybrid region for all arms for genic reads and facet by intra/inter
#'
#' @param genic.dt hybrid data.table
#' @param type type of plot, either stacked or bar
#' @param filename filename for plot with extension
#' @return plot
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import ggplot2
#' @import scales
#' @export

# Requires data.table, ggplot2, scales, RColorBrewer

PlotGenicRegion <- function(genic.dt, type = c("stacked", "bar"), filename) {
  
  plot.dir <- "plots"
  
  # Prepare genic.dt for plotting
  data <- genic.dt
  data[, genic := ifelse(L_gene == R_gene, "intragenic", "intergenic")] # annotate with type of read
  L <- data[, .(genic, L_region)] # split out L and R arms
  setnames(L, "L_region", "region")
  R <- data[, .(genic, R_region)]
  setnames(R, "R_region", "region")
  data <- rbind(L, R)
  
  # Now get proportional table to get relative values
  prop <- melt(prop.table(table(data), 1))
  prop$region <- factor(prop$region, levels = c("UTR5", "CDS", "UTR3", "lincRNA", "miRNA", "other")) # re-order factor for plotting
  prop <- prop[order(prop$region), ] # reorder data.frame by factor for plotting
  
  if(type == "stacked") {
    p <- ggplot(prop, aes(x = genic, y = value, fill = region)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_y_continuous(labels = percent, limits = c(0, 1)) +
      labs(title = "Region of origin of hybrid arm",
           y = "Percentage",
           x = "Type of genic hybrid") +
      theme_classic() +
      scale_fill_manual(values = c("#eff3ff", "#bdd7e7", "#6baed6", "#fb6a4a", "#ef3b2c", "#cb181d"))
    
  } else if(type == "bar") {
    p <- ggplot(prop, aes(x = region, y = value, fill = region)) +
      geom_bar(stat = "identity") +
      facet_grid(genic ~ .) +
      scale_y_continuous(labels = percent, limits = c(0, 1)) +
      labs(title = "Region of origin of hybrid arm",
           x = "Region",
           y = "Percentage") +
      theme_classic() +
      theme(legend.position = "none") +
      scale_fill_manual(values = c("#eff3ff", "#bdd7e7", "#6baed6", "#fb6a4a", "#ef3b2c", "#cb181d"))
  }
  
  ggsave(filename = file.path(plot.dir, filename), plot = p)
  return(p)
  
}

# =============================================================================
# Plot intra-inter ratio (for genic reads)
# =============================================================================

#' Plot intra-inter ratio (for genic reads)
#'
#' Plot the intra-inter ratio (for genic reads)
#'
#' @param genic.dt hybrid data.table of genic reads
#' @param filename filename for plot with extension
#' @return plot
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import ggplot2
#' @export

PlotIntraInterRatio <- function(genic.dt, filename) {
  
  ratio <- data.frame(Type = ifelse(genic.dt$L_gene == genic.dt$R_gene, "Intrahybrid", "Interhybrid"))
  ratio$Type <- factor(ratio$Type, levels = c("Intrahybrid", "Interhybrid"))
  
  plot.dir <- "plots"
  
  # Plot
  p <- ggplot(ratio, aes(x = factor(1), fill = Type)) +
    geom_bar(width = 1) +
    coord_polar(theta = "y") +
    labs(title = "Hybrid read type",
         x = "",
         y = "") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    scale_fill_manual(values = c("#08619c", "#deebf7"))
  ggsave(filename = file.path(plot.dir, filename), plot = p)
  return(p)
}

# =============================================================================
# Plot region mate
# =============================================================================

#' Plot region mate
#'
#' Plot the mate for a given region
#'
#' @param hybrid.dt hybrid data.table
#' @param region region for which to find mates
#' @param filename filename for plot with extension
#' @return plot
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import ggplot2
#' @export

PlotRegionMate <- function(hybrid.dt, region = c("miRNA", "lincRNA", "ncRNA", "CDS", "UTR3", "UTR5", "INTRON", "INTERGENIC", "tRNA", "rRNA", "other"), filename) {
  
  plot.dir <- "plots"
  
  # Manipulate data
  data <- hybrid.dt[, .(L_region, R_region)]
  data[L_region == region, mate := R_region] # get mate when L_region == region
  data[R_region == region, mate := L_region]
  data[!is.na(mate)] # remove others
  
  # Plot
  p <- ggplot(data, aes(x = mate)) +
    geom_bar() +
    labs(title = paste0("Mate for hybrid arm mapped to ", region),
         x = "Mate region",
         y = "Count") +
    theme_classic()
  ggsave(filename = file.path(plot.dir, filename), plot = p)
  return(p)
}

# =============================================================================
# Plot loop length
# =============================================================================

#' Plot loop length
#'
#' Plot the loop length for genic reads (CDS-CDS and UTR3-UTR3)
#'
#' @param genic.dt hybrid data.table of genic reads
#' @param weighting weight by coverage of islands or not
#' @param filename filename for plot with extension
#' @return plot
#'
#' @importFrom data.table ":=" data.table fread melt rbindlist setkey setnames ".N"
#' @import ggplot2
#' @export

PlotLoopLength <- function(genic.dt, weighting = c(TRUE, FALSE), ratio = c(TRUE, FALSE), filename) {
  
  plot.dir <- "plots"
  
  # Calculate loop length and add weighting
  genic.dt[, loop := R_start - L_end - 1] # -1 for 1 based
  setkey(genic.dt, L_gene)
  
  # Add gene region information for relative loop length calculation
  info <- longest.rna.dt[, .(ensembl_gene_id, cds_len, utr3_len)]
  setkey(info, ensembl_gene_id)
  genic.dt <- info[genic.dt]
  
  if(weighting == TRUE) {
    # Expand according to weighting
    genic.dt[, weighting := min(L_peakcoverage, R_peakcoverage), by = id] # weights result by peak coverage
    genic.dt <- genic.dt[rep(1:nrow(genic.dt), genic.dt$weighting, each = TRUE)]
  }
  
  # Select out UTR3-UTR3 and CDS-CDS hybrids and calulate loop ratios
  utr3utr3 <- genic.dt[L_region == "UTR3" & R_region == "UTR3"]
  utr3utr3[, hybrid := "UTR3-UTR3"]
  utr3utr3[, loop_ratio := loop/utr3_len]
  
  cdscds <- genic.dt[L_region == "CDS" & R_region == "CDS"]
  cdscds[, hybrid := "CDS-CDS"]
  cdscds[, loop_ratio := loop/cds_len]
  
  data <- list(utr3utr3 = utr3utr3[, .(hybrid, loop, loop_ratio)],
               cdscds = cdscds[, .(hybrid, loop, loop_ratio)])
  data <- rbindlist(data)
  data <- data[loop > 0]
  
  if(ratio == FALSE) {
    
    print(paste("There are", length(utr3utr3$loop), "UTR3-UTR3 loops"))
    print(paste("There are", length(cdscds$loop), "CDS-CDS loops"))
    print(wilcox.test(utr3utr3$loop, cdscds$loop, conf.int = TRUE))
    p <- ggplot(data, aes(x = log10(loop), fill = hybrid)) +
      geom_density(alpha = 0.3) +
      labs(title = "Distribution of loop lengths",
           x = "log10(loop length)",
           y = "Density") +
      scale_fill_manual(name = "Hybrid", values = c("#67a9cf", "#ef8a62")) +
      theme_classic()
    
  } else if(ratio == TRUE) {
    
    print(paste("There are", length(utr3utr3$loop_ratio), "UTR3-UTR3 loops"))
    print(paste("There are", length(cdscds$loop_ratio), "CDS-CDS loops"))
    print(wilcox.test(utr3utr3$loop_ratio, cdscds$loop_ratio, conf.int = TRUE))
    p <- ggplot(data, aes(x = log10(loop_ratio), fill = hybrid)) +
      geom_density(alpha = 0.3) +
      labs(title = "Distribution of loop length ratios",
           x = "log10(loop length ratio)",
           y = "Density") +
      scale_fill_manual(name = "Hybrid", values = c("#67a9cf", "#ef8a62")) +
      theme_classic()
  }
  
  ggsave(filename = file.path(plot.dir, filename), plot = p)
  return(p)
}
