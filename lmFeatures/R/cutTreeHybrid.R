cutTreeHybrid <-
function(p.){
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
