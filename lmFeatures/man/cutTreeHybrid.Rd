\name{cutTreeHybrid}
\alias{cutTreeHybrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cutTreeHybrid(p.)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p.}{
%%     ~~Describe \code{p.} here~~
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
function (p.) 
{
    dissim1 = dist(p.$mz)
    dendro1 <- hclust(d = dissim1, method = "average")
    ct1 <- cutreeDynamic(dendro1, cutHeight = NULL, minClusterSize = 30, 
        method = "hybrid", deepSplit = 0, pamStage = TRUE, distM = as.matrix(dissim1), 
        maxPamDist = 0, verbose = 0)
    return(ct1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
