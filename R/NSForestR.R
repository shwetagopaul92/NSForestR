#'run RandomForest
#'@param column numeric(1) column number 
#'@param dataFull data.frame tsvfile for from a SingleCellExperiment
#'@param dataDummy matrix dummy columns for one vs all Random Forest modeling
#'@param rfTrees numeric(1) number of trees
#'@param cind integer index of "Clusters" column
#'@export
runRandomForest<-function(column, dataFull, dataDummy, rfTrees, cind=cind){
  require(randomForest)
  x_train=dataFull[,-cind]
  y_train=as.factor(dataDummy[,column])
  rf=randomForest(x=x_train,y=y_train, ntree=rfTrees, importance=TRUE)
  imp=data.frame(importance(rf))
  imp=imp[order(imp["MeanDecreaseGini"], decreasing=TRUE),]
  if(nrow(imp)>30){return(imp[1:30,])}else{return(imp)}
}

#'run negativeOut on output from randomForest()
#'@param impDF data.frame 
#'@param column numeric(1) column number 
#'@param medianValues data.frame  cluster median expression values
#'@param MedianExpressionLevel numeric(1) parameter for filtering and ranking genes
#'@export
negativeOut<-function(impDF, column, medianValues, MedianExpressionLevel){
  keep=sapply(rownames(impDF),function(x){
    medianValues[column,x]>MedianExpressionLevel
  })
  return(impDF[keep,])
}

#'run BinaryScore
#'@param negOutDF data.frame output from negativeOut()
#'@param informativeGenes numeric(1) number of genes to regard as informative for each cluster
#'@param medianValues data.frame cluster median expression values
#'@param column numeric(1) column number
#'@param nclust numeric(1) number of clusters
binaryScore<-function(negOutDF, informativeGenes, medianValues, column, nclust){
  negOutDF=negOutDF[1:informativeGenes,]
  keep=intersect(rownames(negOutDF), names(medianValues))
  subMedians=medianValues[,keep]
  subMediansScaled=data.frame(do.call("rbind",apply(subMedians,1, function(x){
    x/subMedians[column,]
  })))
  subMediansScaled=1-subMediansScaled
  subMediansScaledMod=apply(subMediansScaled,2, function(x){ifelse(x>0,x,0)})
  myColSums=colSums(subMediansScaledMod)
  rescaled=myColSums/nclust
  res=data.frame(symbol=names(rescaled), scaled=rescaled)
  negOutDF$symbol=rownames(negOutDF)
  res=merge(res,negOutDF)
  res$scaledGini=res[["MeanDecreaseGini"]]/sum(res[["MeanDecreaseGini"]])
  res=res[order(res$scaled, res[["MeanDecreaseGini"]], decreasing=TRUE),]
  if(nrow(res)>informativeGenes){return(res[1:informativeGenes,])}else{return(res)}
}

#'to get expression cutoffs for f-beta testing
#'@param binaryDF data.frame 
#'@param clusterCol numeric(1) cluster column 
#'@param dataDummy matrix dummy columns for one vs all Random Forest modeling
#'@param dataFull data.frame tsvfile for from a SingleCellExperiment
dtCutoff<-function(binaryDF, clusterCol, dataDummy, dataFull ){
  require(tree)
  res=
    sapply(binaryDF$symbol,function(x){
      temp=tree(dataDummy[,clusterCol]~dataFull[[as.character(x)]])$frame$splits[1,2]
      temp=sub(as.character(temp),pattern=">",replacement="")
      temp=as.numeric(temp)
    }
    )
  res2=data.frame(symbol=binaryDF$symbol, cutoff=res, stringsAsFactors=FALSE)
  return(res2)
}

#'to make the prediction columns from trees
#'@param cutoffs data.frame holds cutoff values for each symbol
#'@param dataFull data.frame tsvfile for from a SingleCellExperiment
genPredDF<-function(cutoffs, dataFull){
  res=data.frame(matrix(nrow=nrow(dataFull)))
  for(i in 1:nrow(cutoffs)){
    res[,i]=dataFull[[as.character(cutoffs[i,"symbol"])]]>=cutoffs[i,"cutoff"]
  }
  names(res)=cutoffs$symbol
  return(res)
}

#'run fBeta test on prediction results
#'@param predDF data.frame prediction results for each symbol
#'@param clusterName numeric(1) cluster number
#'@param betaValue numeric(1) value for beta weighting
#'@param maxSymbol numeric(1) number of symbols to retain
#'@export
fBetaTest<-function(predDF, clusterNum, betaValue,maxSymbol){
  require(MLmetrics)
  predDF=predDF[,1:maxSymbol]
  mycols=1:ncol(predDF)
  mynams=names(predDF)
  res=lapply(1:ncol(predDF), function(x){
    combs=combn(mycols,x)
    if(x>1){
      combs=data.frame(combs)
      names(combs)=apply(combs,2,function(y){paste0(mynams[y],collapse="_")})
      temp=apply(combs,2,function(y){
        as.numeric(apply(predDF[,y],1,all))
      })
    }else{
      names(combs)=mynams[combs]
      temp=sapply(combs, function(y){as.numeric(predDF[,y])})
    }
    apply(temp,2,function(z){
      FBeta_Score(y_pred=z, positive=1,y_true=dataDummy[,clusterNum],beta=betaValue)
    })
  })
  res
}

#'run NSForestCL 
#'@param file csv file
#'@param clusterColumn numeric(1) cluster column 
#'@param cluster_value numeric(1) value for cluster
#'@param rfTrees numeric(1) number of trees
#'@param medianValues data.frame  cluster median expression values
#'@param informativeGenes numeric(1) number of genes to regard as informative for each cluster
#'@param num_top_ranked numeric(1) number of symbols to retain
#'@param betaValue numeric(1) value for beta weighting
#'@param column numeric(1) column number
#'@export
NSForestCL <- function(file, clusterColumn, clusterValue, rfTrees, informativeGenes, num_top_ranked, betaValue, column){
  require(matrixStats)
  dataFull = read.csv(file, sep="\t",header=TRUE,row.names=1)
  cind=which(names(dataFull)=="Clusters")
  mylev1=levels(dataFull$Clusters)[1]
  haslev1=as.numeric(as.character(dataFull$Clusters)==mylev1)
  dataDummy = cbind(haslev1,model.matrix(~Clusters, data=dataFull)[,-1])
  colnames(dataDummy)=as.character(levels(dataFull$Clusters))
  #Creates matrix of cluster median expression values
  splitdat=split(dataFull,dataFull$Clusters)
  mymeds=lapply(splitdat,function(x){colMedians(as.matrix(x[,-cind]))})
  mynames=names(mymeds)
  mymeds=do.call("rbind",mymeds)
  medianValues=data.frame(mymeds)
  names(medianValues)=names(dataFull)[-cind]
  medianValues$cluster=mynames
  
  res=runRandomForest(column=2, dataFull, dataDummy,  rfTrees=100, cind)
  
  binres=binaryScore(res,informativeGenes=5, medianValues,column=2,nclust=8)
  
  cuts=dtCutoff(binres,clusterCol = 2, dataDummy, dataFull)
  
  preds=genPredDF(cuts, dataFull)
  
  fBetaTest(preds,2, .5, 3)
  
  
}


