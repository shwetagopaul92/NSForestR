#'import python modules with reticulate
#'@export
impMod = function() {
  np <- import("numpy", delay_load=TRUE, convert=FALSE)
  pd <- import("pandas", delay_load=TRUE)
  graphviz <- import("graphviz", delay_load=TRUE)
  numexpr <- import("numexpr", delay_load=TRUE)
  itertools <- import("itertools", delay_load=TRUE)
  list(np=np, pd=pd, graphviz=graphviz, numexpr=numexpr, itertools=itertools)
}