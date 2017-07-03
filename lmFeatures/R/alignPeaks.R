#' Simplest possible feature detection method
#'
#' Gets the highest peak in the dataset and collect all other peaks in its
#' vicinity. Finally filter multiple peaks in one scan by choosing the highest one.
#'
#' @param peaks data.frame with columns 'mz' and 'time' for m/z values and scan ids.
#' @param epsilon width of feature allowed region
#'
#' @return peaks data.frame with column 'ion', which contains feature m/z value
#' for the feature peaks and NA otherwise
#' @export
#' @family feature detection functions
#' @import stats
alignPeaks <-
function(peaks,epsilon=5e-3){
    #Align topmost peak
    idx<-which(is.na(peaks$ion))
    cP<-which.max(peaks$intensity[idx])
    ionI<-idx[cP]
    dmz<-peaks$mz[idx]-peaks$mz[ionI]
    candI<-which(abs(dmz)<epsilon)
    tbl<-which(table(peaks$time[idx[candI]])>1)
    if(length(tbl)>0){
        for(t in unique(as.numeric(names(tbl)))){
            dtidx<-which(peaks$time[idx[candI]]==t)
            candI<-candI[-(order(dmz[candI[dtidx]])[-1])]
        }
    }
    peaks$ion[idx[candI]]<-peaks$mz[ionI]
    return(peaks)
}

#' Feature detection by linear regression
#'
#' Check the precondition for the continious feature and if met fit linear
#' regression to the data and assign to the feature all peaks for which
#' standartized residual is less then \code{\link{thRes}}.
#'
#' The preconditions are:
#' \itemise
#'
#' @param peaks data.frame with columns \code{mz} and \code{time} for m/z values and scan ids.
#' @param epsilon width of feature allowed region. Default 5e-3.
#' @param minAbsLen minimal length for the continious feature in number of scans.
#' Default 300.
#' @param minRelLen minimal length for the continious feature relative to the
#' total length of the registration. Default 5e-1 (50\%).
#' @param minRunLen minimal length of uninterapted nonzero run for the
#' continious feature. Default 100.
#' @param thRes threshold for standartized linear regression value. Default 2.0.
#'
#' @return peaks data.frame with column \code{ion}, which contains feature m/z value
#' for the feature peaks and NA otherwise
#' @family feature detection functions
#' @export
#' @import stats
alignPeaksDT<-function(peaks,
                       epsilon=5e-3,
                       minAbsLen=300,
                       minRelLen=0.5,
                       minRunLen=100,
                       thRes=2.0){
  #Align topmost peak
  idx<-which(is.na(peaks$ion))
  cP<-which.max(peaks$intensity[idx])
  ionI<-idx[cP]
  maxMZ<-peaks$mz[ionI]
  idx<-which(abs(peaks$mz-maxMZ)<(epsilon)&is.na(peaks$ion))
  scanLen<-max(peaks$scan)
  rl<-get_runLen(peaks$scan[idx],scanLen)
  if(length(idx)>=min(minAbsLen,minRelLen*scanLen)|max(rl)>minRunLen){
    lm(mz~time,data=peaks[idx,])->md
    if(any(abs(md$residuals) >thRes * sd(md$residuals))){
      idx<-idx[-which(abs(md$residuals) >thRes * sd(md$residuals))]
      lm(mz~time,data=peaks[idx,])->md1
    }else{
      md1<-md
    }
    peaks$lm[idx]<-md1$coefficients[2]
    peaks$ion[idx]<-md1$coefficients[1]
  }else{
    # peaks$ion[idx]<-maxMZ
    peaks$ion[ionI]<- -Inf
    #peaks<-alignPeaks(peaks,epsilon*0.1)
  }
  return(peaks)
}

#' Feature detection by clusterisation and linear regression
#'
#' @param peaks data.frame with columns 'mz' for the m/z values,
#' 'time' for the scan time value and
#' 'scan' for the scan ids.
#' @param ppm width of feature allowed region
#' @param minAbsLen minimal length for the continious feature in number of scans.
#' See also alignPeaksDT.
#' @param minRelLen minimal length for the continious feature relative to the
#' total length of the registration. See also alignPeaksDT.
#' @param minRunLen minimal length of uninterapted nonzero run for the
#' continious feature. See also alignPeaksDT.
#' @param thRes threshold for standartized linear regression value. See also alignPeaksDT.
#'
#' @return peaks data.frame with column 'ion', which contains feature m/z value
#' for the feature peaks and NA otherwise
#'
#' @export
#' @family feature detection functions
#'
#' @import stats dynamicTreeCut
alignPeaksDTrle <- function(peaks,
                            ppm = 20,
                            minAbsLen=300,
                            minRelLen=0.5,
                            minRunLen=100,
                            thRes=2.0) {
  idx <- which(is.na(peaks$ion))
  cP <- which.max(peaks$intensity[idx])
  ionI <- idx[cP]
  maxMZ <- peaks$mz[ionI]
  idx <- which(abs(peaks$mz - maxMZ) < (ppm * 1e-6 * maxMZ) &
                 is.na(peaks$ion))
  scanLen <- max(peaks$scan)
  p. <- peaks[idx, ]
  rl<-get_runLen(peaks$scan[idx],scanLen)
  if(length(idx) >= min(minAbsLen,minRelLen*scanLen)|max(rl)>minRunLen){
    ct1 <- cutTreeHybrid(p.)
    cct1 <- ct1[which.max(p.$intensity)]
    id.1 <- idx[ct1 == cct1]
    mz1 <- median(peaks$mz[id.1])
    rl <- get_runLen(peaks$scan[id.1], scanLen)
    if (length(id.1) >= min(minAbsLen,minRelLen*scanLen) || max(rl) > minRunLen) {
      peaks$clIon1[id.1] <- mz1
#      peaks$ion[id.1] <- mz1
      lm(mz ~ time, data = peaks[id.1, ]) -> md
      if (any(abs(md$residuals) > thRes * sd(md$residuals))) {
        id.2 <- id.1[-which(abs(md$residuals) > thRes * sd(md$residuals))]
        lm(mz ~ time, data = peaks[id.2, ]) -> md1
      } else{
        md1 <- md
        id.2 <- id.1
      }
      peaks$lm[id.2] <- md1$coefficients[2]
      peaks$clIon2[id.2] <- md1$coefficients[1]
      peaks$ion[id.2] <- md1$coefficients[1]
    } else{
      peaks$ion[ionI] <- -Inf
    }
  } else{
    peaks$ion[ionI] <- -Inf
  }
  return(peaks)
}


#' Feature detection by linear regression with region width in ppm
#'
#' The same functionality as in 'alignPeaksDT' but the region width vary with
#' m/z value
#'
#' @param peaks data.frame with columns 'mz' for the m/z values,
#' 'time' for the scan time value and
#' 'scan' for the scan ids.
#' @param ppm width of feature allowed region in ppm.
#' @inheritParams alignPeaksDT
#' @param minAbsLen minimal length for the continious feature in number of scans.
#' See also alignPeaksDT.
#' @param minRelLen minimal length for the continious feature relative to the
#' total length of the registration. See also alignPeaksDT.
#' @param minRunLen minimal length of uninterapted nonzero run for the
#' continious feature. See also alignPeaksDT.
#' @param thRes threshold for standartized linear regression value. See also alignPeaksDT.
#'
#' @return peaks data.frame with column 'ion', which contains feature m/z value
#' for the feature peaks and NA otherwise
#'
#' @seealso \code{\link{alignPeaksDT}} which this function wraps
#' @export
#' @family feature detection functions
alignPeaksDTppm<-function(peaks,
                          ppm=20,
                          minAbsLen=300,
                          minRelLen=0.5,
                          minRunLen=100,
                          thRes=2.0){
  #Align topmost peak
  idx<-which(is.na(peaks$ion))
  cP<-which.max(peaks$intensity[idx])
  ionI<-idx[cP]
  maxMZ<-peaks$mz[ionI]
  epsilon<-ppm*1e-6*maxMZ
  return(alignPeaksDT(peaks=peaks,
                      epsilon=epsilon,
                      minAbsLen=minAbsLen,
                      minRelLen=minRelLen,
                      minRunLen=minRunLen,
                      thRes=thRes))
}

