\name{alignPeaks}
\alias{alignPeaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
alignPeaks(.peaks, .epsilon = 0.005, .epsC13 = 5e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.peaks}{
%%     ~~Describe \code{.peaks} here~~
}
  \item{.epsilon}{
%%     ~~Describe \code{.epsilon} here~~
}
  \item{.epsC13}{
%%     ~~Describe \code{.epsC13} here~~
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
function (.peaks, .epsilon = 0.005, .epsC13 = 5e-04) 
{
    idx <- which(is.na(.peaks$ion))
    cP <- which.max(.peaks$intensity[idx])
    ionI <- idx[cP]
    dmz <- .peaks$mz[idx] - .peaks$mz[ionI]
    candI <- which(abs(dmz) < .epsilon)
    tbl <- which(table(.peaks$time[idx[candI]]) > 1)
    if (length(tbl) > 0) {
        for (t in unique(as.numeric(names(tbl)))) {
            dtidx <- which(.peaks$time[idx[candI]] == t)
            candI <- candI[-(order(dmz[candI[dtidx]])[-1])]
        }
    }
    .peaks$ion[idx[candI]] <- .peaks$mz[ionI]
    return(.peaks)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
