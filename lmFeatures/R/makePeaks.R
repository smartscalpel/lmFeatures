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
    
    for (i in seq(along = sel$scanidx)) {
        scan <- as.data.table(getScan(xraw, sel$scanidx[i], sel$mzrange))
        p<-scan
        p$ion<-NA
        p$MP<-NA
        p$time<-xraw@scantime[i]
        p$scan<-i
        p$delta<-NA
        p$adduct<-NA
        peaks<-rbind(peaks,p)
    }
    return(peaks)
}
