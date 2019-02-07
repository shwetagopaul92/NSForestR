#' transform sce to tsv file
#' @param sce SingleCellExperiment object
#' @param destination character filepath to tsv
#' @param clustvbl character(1) variable in colData to use as cluster label
#' @examples 
#' requireNamespace("pcmp")
#' data(pcmp::Ab10k)
#' tf = tempfile(fileext=".tsv")
#' sceTotsv(Ab10k,tf)
#' head(read.delim(tf,sep="\t"))
#' @export
sceTotsv <- function(sce, destination, clustvbl){
  asy = t(assay(sce))
  assay.df = as.data.frame(asy)
  clusters = colData(sce)[[clustvbl]]
  assay.df$Clusters = clusters
  assay.df$Clusters <- sub("^","CL",assay.df$Clusters)
  write.table(assay.df,file=destination,quote=FALSE,sep="\t",col.names = NA)
}
