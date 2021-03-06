#' Converts XCMS data into MZ table
#'
#'Creates data.table with five columns: 'mz', 'intensity', 'ion', 'time', 'scan'
#' @param xcms xcmzRaw data to convert
#'
#' @return ata.table with five columns: 'mz', 'intensity', 'ion', 'time', 'scan'
#' @export
#' @import xcms data.table
#'
makePeaks <- function(xcms){
    sel <- profRange(xcms)
    peaks<-data.frame(mz=100.0,
                      intensity=1.0e6,
                      ion=NA,
                      time=0.0,
                      scan=1)[FALSE,]
    l<-list()
    l<-lapply(sel$scanidx,function(i){
      scan <- as.data.table(getScan(xraw, sel$scanidx[i], sel$mzrange))
      p<-scan
      p$ion<-NA
      p$time<-xraw@scantime[i]
      p$scan<-i
      return(p)
    })
    peaks<-as.data.table(ldply(l,rbind))
    return(peaks)
}
