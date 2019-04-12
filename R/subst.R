#Code copied and adapted from https://github.com/JCVenterInstitute/NSForest.
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

#' run NSForest python code on data emitted by this function, and collect the results
#' use source_python from reticulate in this version
#' @param sce instance of SingleCellExperiment
#' @param clustVblName character(1) component of colData(sce) to use for cluster label
#' @param numtrees numeric(1) number of trees built in random forest step
#' @param numthread numeric(1) number of threads for python execution
#' @param num_inf_genes numeric(1) number of genes to regard as informative for each cluster
#' @param num_top_ranked numeric(1) number of genes to retain in reporting
#' @examples
#' requireNamespace("pcmp") # from github vjcitn/pcmp
#' sce = pcmp::sce300xx
#' sce$Clusters = as.numeric(sce$PIMMquart) # grouped by percent reads aligned to introns
#' sds = matrixStats::rowSds(assay(sce))
#' ok = which(sds > quantile(sds, .95))
#' r1 = runNSForest(sce[ok,], num_inf_genes=8, num_top_ranked=5)
#' str(r1)
#'@export
runNSForest <- function(sce, clustVblName="Clusters",
numtrees=100, numthread=1, num_inf_genes = 5,
 num_top_ranked = 3) {
  if (clustVblName != "Clusters") {
    if ("Clusters" %in% names(colData(sce)))
        warning(paste("replacing content of 'Clusters' in colData with contents of '", clustVblName,"'"))
        sce$Clusters = sce[[clustVblName]]
    }
  tsvtarg = tempfile(fileext=".tsv")
  chk = sceTotsv(sce, tsvtarg, "Clusters")
  pytarg = tempfile(fileext=".py")
  mydir = tempdir()
  prog = doSubs(tsvtarg, clustVblName="Clusters",
     numtrees=numtrees, numthread=numthread, num_inf_genes=
       num_inf_genes, num_top_ranked=num_top_ranked)
  writeLines(prog, pytarg)
  workdir = tempdir()
  curd = getwd()
  on.exit(setwd(curd))
  setwd(workdir)
  source_python(pytarg)  # check errors?
  output_csv = list.files(pattern="*.csv")
  lapply(output_csv, read.csv)
}


# version 2 -- minimize the use of files, and allow proper handling of flexible selection of cluster table
