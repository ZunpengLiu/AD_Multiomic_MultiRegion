#library(devtools, lib.loc="/net/bmc-lab6/data/lab/kellis/users/ssz/software/Miniconda3/envs/devtools/lib/R/library")

Run_SCAVENGE <- function(SE_Data, trait_file){
  # trait_file: <chr> <snp_loc> <snp_loc> <snp_id> <posterior probability> <filename>
  require(SCAVENGE)
  require(chromVAR)
  require(gchromVAR)
  require(BuenColors)
  require(SummarizedExperiment)
  require(data.table)
  require(BiocParallel)
  require(BSgenome.Hsapiens.UCSC.hg38)
  require(dplyr)
  require(igraph)
  require(parallel)
  get_sigcell_simple <- function (knn_sparse_mat = mutualknn30, seed_idx = seed_p0.05,
      topseed_npscore = topseed_npscore, permutation_times = 1000,
      true_cell_significance = 0.05, rda_output = F, out_rda = "true_cell_df.rda",
      mycores = 4, rw_gamma = 0.05){
      if (permutation_times < 100) {
          warning("Permutation times less than 100")
      }
      stopifnot(`true_cell_significance must be a numeric value between 0, 1` = (true_cell_significance >
          0 & true_cell_significance < 1))
      message("Get started!")
      cell_mat <- data.frame(cell = 1:nrow(knn_sparse_mat), degree = colSums(knn_sparse_mat))
      cell_table <- data.frame(table(cell_mat$degree))
      seed_mat_top <- data.frame(seed = which(seed_idx), degree = colSums(knn_sparse_mat[, which(seed_idx==TRUE)]))
      seed_table_top <- data.frame(table(seed_mat_top$degree))
      xx_top <- tapply(cell_mat[, 1], cell_mat[, 2], list)
      xx2_top <- xx_top[names(xx_top) %in% seed_table_top$Var1]
      permutation_score_top <- mclapply(1:permutation_times, mc.cores = mycores,
          function(i) {
              sampled_cellid <- xx2_top %>% mapply(sample, ., seed_table_top$Freq) %>%
                  unlist %>% sort
              xx <- randomWalk_sparse(intM = knn_sparse_mat, queryCells = rownames(knn_sparse_mat)[as.numeric(sampled_cellid)],
                  gamma = rw_gamma)
              if (i%%100 == 0) {
                  message(i)
              }
              return(xx)
          })
      names(permutation_score_top) <- paste0("permutation_", 1:permutation_times)
      permutation_score_top <- data.frame(sapply(permutation_score_top,
          c))
      permutation_df_top <- data.frame(matrix(nrow = nrow(knn_sparse_mat),
          ncol = permutation_times))
      permutation_df_top <- apply(permutation_score_top, 2, function(x) {
          temp <- x > topseed_npscore
          return(temp)
      })
      message("cells passed 0.001 threshold: ", round(sum(rowSums(permutation_df_top) <=
          0.001 * permutation_times) * 100/nrow(permutation_df_top),
          2), "%")
      message("cells passed 0.01 threshold: ", round(sum(rowSums(permutation_df_top) <=
          0.01 * permutation_times) * 100/nrow(permutation_df_top),
          2), "%")
      message("cells passed 0.05 threshold: ", round(sum(rowSums(permutation_df_top) <=
          0.05 * permutation_times) * 100/nrow(permutation_df_top),
          2), "%")
      true_cell_top_idx <- rowSums(permutation_df_top) <= true_cell_significance *
          permutation_times
      message("your emprical P value threshold: ", true_cell_significance)
      message("what propertion of enriched cells over all cells: ",
          round((sum(true_cell_top_idx) * 100)/nrow(permutation_df_top),
              2), "%")
      message("what propertion of seed cells that are enriched cells: ",
          round(sum(true_cell_top_idx & seed_idx) * 100/sum(seed_idx),
              2), "%")
      message("fold of true cell over seed: ", round(sum(true_cell_top_idx)/sum(seed_idx),
          2))
      true_cell_top_filter_idx <- ture_cell_df <- data.frame(seed_idx,
          true_cell_top_idx)
      if (rda_output) {
          save(ture_cell_df, permutation_score_top, file = out_rda)
      }
      return(ture_cell_df)
  }
  SE_Data <- addGCBias(SE_Data, genome = BSgenome.Hsapiens.UCSC.hg38)
  trait_import <- importBedScore(rowRanges(SE_Data), trait_file, colidx = 5)
  SE_Data_bg <- getBackgroundPeaks(SE_Data, niterations=200)
  SE_Data_DEV <- computeWeightedDeviations(SE_Data, trait_import, background_peaks = SE_Data_bg)
  z_score_mat <- data.frame(colData(SE_Data), z_score=t(assays(SE_Data_DEV)[["z"]]) %>% c)
  z_score_mat <- z_score_mat[complete.cases(z_score_mat),]
#  print(paste("counts:",ncol(counts)))
#  print(paste("z_score_mat:",nrow(z_score_mat)))
  seed_idx <- seedindex(z_score_mat$z_score, 0.05)
  scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)
  peak_by_cell_mat <- assay(SE_Data)
  tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)
  lsi_mat <- do_lsi(tfidf_mat, dims=30)
  mutualknn30 <- getmutualknn(lsi_mat, 30)
  np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)
  omit_idx <- np_score==0
  sum(omit_idx)
  mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
  np_score <- np_score[!omit_idx]
  TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
  TRS <- TRS * scale_factor
  mono_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)
  # please set @mycores >= 1 and @permutation_times >= 1,000 in the real setting
  mono_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, seed_idx=mono_mat$seed_idx, topseed_npscore=mono_mat$np_score, permutation_times=1000, true_cell_significance=0.05, rda_output=F, mycores=8, rw_gamma=0.05)
  mono_mat2 <- data.frame(mono_mat, mono_permu)
  mono_mat2$rev_true_cell_top_idx <- !mono_mat2$true_cell_top_idx
  return(mono_mat2)
}
library(gridExtra)
library(Signac)
library(Seurat)
library(dplyr)
library(rlist)
library(reshape2)
