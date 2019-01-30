#' transform sce to tsv file
#' @param sce SingleCellExperiment object
#' @param destination character filepath to tsv
#' @examples 
#' require(pcmp)
#' data(Ab10k)
#' sceTotsv(Ab10k,"temp/Ab10k.tsv")
#' @export
sceTotsv <- function(sce, destination){
  asy = t(assay(sce))
  assay.df = as.data.frame(asy)
  clusters = colData(sce)$km9_assay
  assay.df$Clusters = clusters
  assay.df$Clusters <- sub("","CL",assay.df$Clusters)
  write.table(assay.df,file=destination,quote=FALSE,sep="\t",col.names = NA)
}