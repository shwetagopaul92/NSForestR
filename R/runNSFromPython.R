# Code copied and adapted from https://github.com/JCVenterInstitute/NSForest
#'to run NS Forest program from python 
#' @require reticulate 
#' @param tsvfile 
#' @param rfTrees numeric(1) number of trees built in random forest step
#' @param threads numeric(1) number of threads for python execution
#' @param Median_Expression_Level numeric(1) median expression level
#' @param InformativeGenes numeric(1) number of genes to regard as informative for each cluster
#' @param Genes_to_testing numeric(1) number of genes to retain in reporting
#' @param betaValue numeric(1) set values for fbeta weighting
#' @examples 
#' tsvfile = system.file("data/Ab10k.tsv", package="NSForestR")
#' runNSFromPython(tsvfile, rfTrees=100L, threads=1L, Median_Expression_Level=0L, InformativeGenes=5L,
#'                 Genes_to_testing=3L, betaValue=0.5)
#' @export
runNSFromPython <- function(tsvfile, rfTrees, threads, Median_Expression_Level, 
                            InformativeGenes, Genes_to_testing, betaValue) {
  allFunctions_file = system.file("python/allFunctions.py", package="NSForestR")
  coreAnalysis_file = system.file("python/coreAnalysis.py", package="NSForestR")
  source_python(allFunctions_file)
  source_python(coreAnalysis_file)
  runNSForest(tsvfile, rfTrees, threads, Median_Expression_Level, InformativeGenes, Genes_to_testing, betaValue)
}