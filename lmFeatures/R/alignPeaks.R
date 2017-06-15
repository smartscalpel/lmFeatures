alignPeaks <-
function(.peaks,.epsilon=5e-3,.epsC13=5e-4){
    #Align topmost peak
    idx<-which(is.na(.peaks$ion))
    cP<-which.max(.peaks$intensity[idx])
    ionI<-idx[cP]
    dmz<-.peaks$mz[idx]-.peaks$mz[ionI]
    candI<-which(abs(dmz)<.epsilon)
    tbl<-which(table(.peaks$time[idx[candI]])>1)
    if(length(tbl)>0){
        for(t in unique(as.numeric(names(tbl)))){
            dtidx<-which(.peaks$time[idx[candI]]==t)
            candI<-candI[-(order(dmz[candI[dtidx]])[-1])]
        }
    }
    .peaks$ion[idx[candI]]<-.peaks$mz[ionI]
    #find isotopes for topmost peak
    # for(i in candI){
    #   tidx<-which(.peaks$time[idx]==.peaks$time[idx[i]])
    #   dmz<-.peaks$mz[idx[tidx]]-.peaks$mz[idx[i]]
    #   for(c13 in 1:3){
    #     cidx<-which(abs(dmz-1.003*c13)<.epsC13)
    #     if(length(cidx)==1){
    #       .peaks$ion[idx[tidx[cidx]]]<-.peaks$ion[idx[i]]
    #       .peaks$delta[idx[tidx[cidx]]]<-c13
    #     }else if(length(cidx)>1){
    #       cidx<-which.min(abs(dmz-1.003*c13))
    #       .peaks$ion[idx[tidx[cidx]]]<-.peaks$ion[idx[i]]
    #       .peaks$delta[idx[tidx[cidx]]]<-c13
    #     }
    #   }
    #   #cat(paste0(i,'\n'))
    # }
    return(.peaks)
}
