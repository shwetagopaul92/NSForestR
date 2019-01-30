#' transform fixed template into program using user parameter settings
#' @param tsvname character(1)
#' @export
doSubs = function(tsvname, clustVblName="Clusters",
                  numtrees=100, numthread=1, num_inf_genes = 5,
                  num_top_ranked = 3) {
  #
  # grab template file
  templ = readLines(system.file("python/NSforestTemplate.py", package="NSForestR"))
  prog = gsub("%%TSVNAME%%", tsvname, templ)
  prog = gsub("%%CLUSTVBLNAME%%", clustVblName, prog)
  prog = gsub("%%NUMTREES%%", numtrees, prog)
  prog = gsub("%%NUMTHREADS%%", numthread, prog)
  prog = gsub("%%NUM_INF_GENES%%", num_inf_genes, prog)
  prog = gsub("%%NUM_TOP_RANKED%%", num_top_ranked, prog)
  chk = grep("%%", prog)
  stopifnot(length(chk)==0)
  prog
}

# generate the TSV File from an SCE

# invoke program having written output of doSubs to a tempfile
# where you source it  -- run it in a temp folder and get
# all the csv read in and put in a list

#'run doSubs to generate python program and source
#'from reticulate
#'@export
runNSForest <- function(){
  writeLines(doSubs("Ab10k.tsv"), "temp/NSForestprogram.py")
  setwd("temp")
  source_python("NSForestprogram.py")
  output_csv = list.files(pattern="*.csv")
  myfiles = lapply(output_csv, read.delim)
}

#readFiles <- function(){
 # f1 <- read.delim(file='temp/Binary_scores_Supplmental_results.csv')
#  f2 <- read.delim(file='temp/NS-Forest_v2_results.csv')
#  f3 <- read.delim(file='temp/NSForest_v2_topResults.csv')
#  f4 <- read.delim(file='temp/Function_medianValues.csv')
#  f5 <- read.delim(file='temp/NSForest_v2_maxF-scores.csv')
#}

# version 2 -- minimize the use of files
