makePeaks <-
function(.xcms){
    sel <- profRange(.xcms)
    peaks<-data.frame(mz=100.0,
                      intensity=1.0e6,
                      ion=NA,
                      MP=NA,
                      time=0.0,
                      scan=1,
                      delta=NA,
                      adduct=NA)[FALSE,]
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
