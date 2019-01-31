#' transform fixed template into program using user parameter settings
#' @param assaytsv character(1) arbitrary name used for image of assay data in a tsv format as assumed by the NSforest python code
#' @export
doSubs = function(assaytsv, clustVblName="Clusters",
                  numtrees=100, numthread=1, num_inf_genes = 5,
                  num_top_ranked = 3) {
  # grab template file
  templ = readLines(system.file("python/NSforestTemplate.py", package="NSForestR"))
  prog = gsub("%%TSVNAME%%", assaytsv, templ)
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
#'@param tsvfile character filename
#'@export
runNSForest <- function(sce, clustVblName="Clusters",
numtrees=100, numthread=1, num_inf_genes = 5,
 num_top_ranked = 3) {
  tsvtarg = tempfile(fileext=".tsv")
  chk = sceTotsv(sce, tsvtarg)
  pytarg = tempfile(fileext=".py")
  mydir = tempdir()
  prog = doSubs(tsvtarg, clustVblName=clustVblName, 
     numtrees=numtrees, numthread=numthread, num_inf_genes=
       num_inf_genes, num_top_ranked=num_top_ranked)
  writeLines(prog, pytarg)
  workdir = tempdir()
  curd = getwd()
  on.exit(setwd(curd))
  setwd(workdir)
  source_python(pytarg)  # check errors?
  output_csv = list.files(pattern="*.csv")
  lapply(output_csv, read.delim)
}


# version 2 -- minimize the use of files
