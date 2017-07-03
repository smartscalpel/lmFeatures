
#' Calculate Rle representation of feature in time domain.
#'
#' Function creates vector of zeros and for each nonzero scan
#' sets the value to 1. The vector obtained converted into Rle
#' representation and returned
#'
#' @param scan vector indices of nonzero intensity values
#' @param length total number of scans in the spectrum
#'
#' @return Rle representation of the feature in time domain
#' @export
#'
get_Rle<-function(scan,length){
  require(IRanges)
  xf <- rep(0, length)
  xf[scan] <- 1
  Rle(xf) -> r
  return(r)
}

#' Calculate Rle representation of feature in time domain with gaps.
#'
#' Function takes Rle representation created by get_Rle and fill gaps
#' up to gap in length.
#'
#' @param scan vector indices of nonzero intensity values (see get_Rle)
#' @param length total number of scans in the spectrum (see get_Rle)
#' @param gap max length of gap between two nonzero runs
#' @param factor
#'
#' @return with gaps filled
#' @import S4Vectors
get_Rle_nogaps<-function(scan,length,gap=1,factor=10){
  r<-get_Rle(scan,length)
  rll<-runLength(r)
  rl1<-which(rll<=gap)
  rv0<-which(runValue(r)==0)
  if(any(rl1%in%rv0)){
    rl1v0<-rl1[rl1%in%rv0]
    seq<-cumsum(rll)
    check<-function(.x){
      l<-FALSE
      r<-l
      if(.x==1){
        l<-TRUE
      }else if(rll[.x-1]>=factor*gap){
        l<-TRUE
      }
      if(.x==length(rll)){
        r<-TRUE
      }else if(rll[.x+1]>=factor*gap){
        r<-TRUE
      }
      if(l&r){
        return(1)
      }else{
        return(0)
      }
    }
    subs<-unlist(sapply(rl1v0,check))
    r[seq[rl1v0]]<-subs
  }
  return(r)
}

#' Vector of length of feature nonzero runs
#'
#' @param scan vector indices of nonzero intensity values (see get_Rle)
#' @param length total number of scans in the spectrum (see get_Rle)
#'
#' @return vector of length of continious nonzero runs
#' @export
#' @import S4Vectors
get_runLen <- function(scan, length) {
  r<-get_Rle(scan,length)
  return(runLength(r)[which(runValue(r)==1)])
}



freqTable<-function(x,breaks,labels){
  if(length(breaks)-length(labels)!=1)
    stop("labels should be 1 element shorter then breaks \n")

  hist(x,breaks,plot=FALSE)->hrl
  res<-t(data.frame(h=hrl$counts,row.names = labels))

  return(res)
}

runLenTable<-function(scans,scanLen){
  rl <- get_runLen(scans, scanLen)
  return(prepRunLenTable(rl))
}

gapRunLenTable<-function(scans,scanLen){
  r<-get_Rle_nogaps(scans,scanLen)
  rl<-runLength(r)[which(runValue(r)==1)]
  return(prepRunLenTable(rl,prefix = 'gl'))
}

scanOverlapTable<-function(scan){
  breaks<-c(0:5,Inf)
  labels<-c(paste0('ns',1:5),'ns>5')
  return(freqTable(table(scan),breaks = breaks,labels = labels))
}

prepRunLenTable<-function(rl,scanLen,prefix='rl'){
  breaksD<-c(0,10,50,100,500,1000,5000,10000,50000)
  breaks<-c(breaksD[breaksD<scanLen],scanLen)
  labels<-paste0(prefix,breaks[-1])
  return(freqTable(rl,breaks,labels))
}


rlTableFun <- function(.x,scanLen) {
  return(cbind(
    scanOverlapTable(.x$scan),
    gapRunLenTable(.x$scan, scanLen = scanLen),
    runLenTable(.x$scan, scanLen)
  ))
}
