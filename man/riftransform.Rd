\name{riftransform}
\alias{riftransform}

\title{
RIF Transformation
}
\description{
Performs recentered influence function (RIF) transformation on continuous response variable at each specified quantile
}
\usage{
riftransform(y,taus)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Numeric vector containing continuous response variable of interest}
  \item{taus}{Vector of quantiles at which to perform transformation}
}
\details{
}
\value{
Matrix of RIF transformed response in rows and quantiles in columns.
}
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
