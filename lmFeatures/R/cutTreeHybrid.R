#' Get clusters in m/z domain
#'
#' Function calculates hierarchical clustering with 'average' method by hclust.
#' The number of clusters and their composition is detected by cutreeDynamic function.
#'
#' @param p. data.frame with column 'mz' containing m/z values to cluster
#'
#' @return result of cutreeDynamic call
#' @export
#' @import dynamicTreeCut
#'
#' @examples
#'
#' p<-data.frame(mz=c(1.1*runif(100),1.5*runif(100),3.1*runif(100),5.1*runif(100),))
#' p$time<-seq_along(p$mz)
#' cth<-cutTreeHybrid(p)
#' .pal <- brewer.pal(12,'Paired')
#' p$color<-.pal[cth+1]
#' qplot(time,mz,data = p.,color=color)
#'
cutTreeHybrid<-function(p.){
  dissim1 = dist(p.$mz)

  dendro1 <- hclust(d = dissim1, method = 'average')
  ct1 <- cutreeDynamic(
    dendro1,
    cutHeight = NULL,
    minClusterSize = 30,
    method = "hybrid",
    deepSplit = 0,
    pamStage = TRUE,
    distM = as.matrix(dissim1),
    maxPamDist = 0,
    verbose = 0
  )

  return(ct1)
}
