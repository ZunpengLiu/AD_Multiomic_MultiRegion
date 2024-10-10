# Load required libraries
library(devtools)
library(data.table)
library(GenomicRanges)
library(HMMt)
library(stringr)

# HMM function to run Hidden Markov Model on genomic regions
HMM <- function(normalized, na_solution, file) {
  
  # Validate the 'na_solution' input
  if (!na_solution %in% c("NA", "keep", "-")) {
    stop("Unknown na_solution")
  }
  
  # Add index to retain the original order
  normalized.tmp <- normalized
  normalized.tmp$idx <- 1:nrow(normalized)
  
  # Sort by chromosome and start position
  normalized.tmp <- normalized.tmp[order(normalized.tmp$chr, normalized.tmp$start), ]
  
  # Use the bridge function to compute a sequence over genomic regions
  br <- bridge(normalized.tmp[, 1:4])
  
  # Run Baum-Welch algorithm for Hidden Markov Model training
  sink(file=file)
  fit <- BaumWelchT(x=br$x, series.length=br$series.length)
  sink()
  
  # Identify the most likely hidden state (active vs. repressive)
  boundstate <- which.max(fit$mu)
  model <- 0 + (fit$ViterbiPath[br$nonvirtuals] == boundstate)
  model <- ifelse(model == 1, "Active", "Repressive")
  
  # Create a DataFrame containing the results
  df.hmm <- cbind(normalized.tmp[, c(1:3, 5)], model)
  
  # Handle NA values based on the 'na_solution' parameter
  if (na_solution == "NA") {
    df.hmm$model[is.na(normalized.tmp$score)] <- NA
  } else if (na_solution == "-") {
    df.hmm$model[is.na(normalized.tmp$score)] <- "-"
  }
  
  # Restore the original order of the data
  df.hmm <- df.hmm[order(df.hmm$idx), ]
  
  return(df.hmm[, c(1:3, 5)])  # Return only the required columns
}

# Function to get the active or repressive state ranges from GRanges object
getHMMRanges <- function(gr, i = 1, score = 1) {
  gr <- gr[!is.na(mcols(gr)[, i])]  # Filter out NAs
  gr <- gr[mcols(gr)[, i] == score]  # Keep only regions with the desired score
  return(reduce(gr))  # Reduce overlapping or adjacent regions
}

# Define chromosome list
chr.list <- paste0("chr", c(seq(1, 22, 1), "X"))

# Read all the samples (cell types) from files ending in ".bed"
samples <- list.files(pattern = ".bed")
cell_types <- samples
ratio.df <- data.frame()

# Process each sample (cell type)
for (cell in cell_types) {
  
  tmp <- fread(cell)  # Load the bed file
  tmp <- tmp[tmp$V1 %in% chr.list, ]  # Filter for chromosomes in chr.list
  
  # Remove ".100k.bed" from the sample name
  cell2 <- str_remove(cell, pattern = ".100k.bed")
  
  # Scale the score column and update the chromosome with sample name
  tmp$V1 <- paste0(cell2, ";", tmp$V1)
  tmp$V4 <- scale(tmp$V4)
  
  # Append to the ratio DataFrame
  ratio.df <- rbind(ratio.df, tmp)
  
  # Store the DataFrame in the environment (optional)
  assign(cell2, tmp)
}

# Rename columns of the ratio DataFrame
colnames(ratio.df) <- c("chr", "start", "end", "score")

# Convert ratio.df to GRanges object
ratio.range <- GRanges(
  seqnames = ratio.df$chr,
  ranges = IRanges(start = ratio.df$start + 1, end = ratio.df$end),
  score = ratio.df$score
)

# Read the high-quality 100kb bins
bin100K <- read.csv("hg38_100kb.HighQuality.bed", sep = "\t", header = F)

# Prepare DataFrame of 100kb bins, filter by chromosome list
bin100K.df <- bin100K[, c(1, 2, 3)]
bin100K.df <- bin100K.df[bin100K.df$V1 %in% chr.list, ]

# Replicate the 100kb bins for all cell types
bin100K.df2 <- data.frame(V1 = rep(bin100K.df$V1, length(cell_types)),
                          V2 = rep(bin100K.df$V2, length(cell_types)),
                          V3 = rep(bin100K.df$V3, length(cell_types)))

# Create a combined sample/chromosome name for each bin
sample_chr <- c()
for (cell in cell_types) {
  cell2 <- str_remove(cell, pattern = ".100k.bed")
  sample_chr <- c(sample_chr, paste0(cell2, ";", bin100K.df$V1))
}

# Add combined sample/chromosome names to bin100K DataFrame
bin100K.df2$V1 <- sample_chr
colnames(bin100K.df2) <- c('chrom', 'start', 'end')
bin100K.df2$start <- bin100K.df2$start + 1
bin100K.df2$queryHits <- seq(1, nrow(bin100K.df2), 1)

# Convert bin100K.df2 to GRanges object
bin100K.range <- GRanges(
  seqnames = bin100K.df2$chrom,
  ranges = IRanges(start = bin100K.df2$start, end = bin100K.df2$end)
)

# Find overlaps between 100kb bins and the ratio ranges
overlap.df <- data.frame(findOverlaps(bin100K.range, ratio.range))
overlap.df$score <- ratio.df[overlap.df$subjectHits, ]$score

# Merge overlap scores with 100kb bins
bin100K.df <- merge(bin100K.df2,
                    overlap.df[, c('queryHits', 'score')],
                    by = 'queryHits', all.x = TRUE)

# Write the combined data to file
write.table(bin100K.df,
            "./All.ClassBrainRegion_Pathology.Zscore.tsv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Prepare bin100K ratio DataFrame for HMM
bin100K.ratio.df <- bin100K.df[, 2:5]

# Run HMM on the ratio data
hmm_calls.1 <- HMM(bin100K.ratio.df, na_solution = "NA", file = "./Active_compartment.stat")

# Extract active compartment ranges
hmm_ranges.1 <- getHMMRanges(as(hmm_calls.1, "GRanges"), score = "Active")

# Write the active compartments to a bed file
write.table(data.frame(hmm_ranges.1)[, 1:3],
            "./All.Active_compartment.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Extract active compartment data and split into chromosome/sample name
Active <- data.frame(hmm_ranges.1)[, 1:3]
library(stringr)
Active$chr <- unlist(str_split(Active$seqnames, ";"))[seq(2, nrow(Active) * 2, 2)]
Active$sample <- unlist(str_split(Active$seqnames, ";"))[seq(1, nrow(Active) * 2, 2)]

# Write active compartments for each cell type
for (sa in cell_types) {
  cell2 <- str_remove(sa, pattern = ".100k.bed")
  write.table(Active[Active$sample == cell2, c(4, 2, 3)],
              paste0(cell2, ".Active.100k.bed"),
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
}

# Save the R session workspace
save.image("Epigenomic_Compartments.RData")

# System commands to find repressive compartments
# (Run in terminal)
'''
for i in `ls *.Active.100k.bed`; do
    bedtools intersect \
        -a hg38_100kb.HighQuality.bed \
        -b $i -v | bedtools merge > ${i%.Active.100k.bed}.Repressive.100k.bed
done
'''
