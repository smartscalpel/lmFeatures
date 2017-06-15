\name{makePeaks}
\alias{makePeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makePeaks(.xcms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.xcms}{
%%     ~~Describe \code{.xcms} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (.xcms) 
{
    sel <- profRange(.xcms)
    peaks <- data.frame(mz = 100, intensity = 1e+06, ion = NA, 
        MP = NA, time = 0, scan = 1, delta = NA, adduct = NA)[FALSE, 
        ]
    for (i in seq(along = sel$scanidx)) {
        scan <- as.data.table(getScan(xraw, sel$scanidx[i], sel$mzrange))
        p <- scan
        p$ion <- NA
        p$MP <- NA
        p$time <- xraw@scantime[i]
        p$scan <- i
        p$delta <- NA
        p$adduct <- NA
        peaks <- rbind(peaks, p)
    }
    return(peaks)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
