#'to run NS Forest program from python 
#' @param sce instance of SingleCellExperiment
#' @param numtrees numeric(1) number of trees built in random forest step
#' @param numthread numeric(1) number of threads for python execution
#' @param expressionLevel numeric(1) median expression level
#' @param num_inf_genes numeric(1) number of genes to regard as informative for each cluster
#' @param num_top_ranked numeric(1) number of genes to retain in reporting
#' @param beta_value numeric(1) set values for fbeta weighting
#' @examples 
#' require("pcmp")
#' data(Ab10k)
#' r1 = runNSFromPython(Ab10k) # running with default parameters
#' @export
runNSFromPython <- function(sce,numtrees=100, numthread=1, expressionLevel=0, 
                            num_inf_genes=5, num_top_ranked=3, beta_value=0.5){
  workdir = tempdir()
  setwd(workdir)
  tsvtarg = tempfile(fileext=".tsv")
  chk = sceTotsv(sce, tsvtarg, "km9_assay")
  file = system.file("python/coreAnalysis.py", package="NSForestR")
  cmd = paste0("python ",file," -fileName ",tsvtarg," -rfTrees ",numtrees," -threads ",numthread," -Median_Expression_Level ",expressionLevel,
               " -InformativeGenes ",num_inf_genes," -Genes_to_testing ",num_top_ranked)
  system(cmd)
  output_csv = list.files(pattern="*.csv")
  lapply(output_csv, read.csv)
}