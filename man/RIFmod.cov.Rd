\name{RIFmod.cov}
\alias{RIFmod.cov}
\title{
RIF Model Covariance
}
\description{
Extracts covariance matrix on RIF models over quantiles
}
\usage{
RIFmod.cov(object, pred)
}

\arguments{
  \item{object}{An object of class multivariate lm representing lm models fit using RIF transformed response variable at each quantile.}
  \item{pred}{Name (character or string) of the genotype variable (e.g. SNP)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A covariance matrix of tau dimension}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Akram Alyass
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
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
