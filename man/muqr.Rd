\name{muqr}
\alias{muqr}
\title{
CQR Meta-Regression
}
\description{
This function fits a meta-regression model of CQR estimates over quantiles. It takes in quantile estimates from UQR models across tau quantiles as well a corresponding covariance matrix said estimates across tau and fits a metaregression model to estimate linear change in quantile estimates (dependent variable) over tau (independent variable).
}
\usage{
muqr(Beta,COV,taus,test="z")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Beta}{Numerical vector of unconditional quantile estimates across response quantiles (tau).}
  \item{COV}{Covariance matrix of Betas over taus}
  \item{taus}{vector of percentiles at which unconditional quantile estimates was fitted}
  \item{test}{Type of test. Default is "z". "t" is also available but not recomended as it has not been tested.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of containing the results of the meta-regression model estimating the effects of tau unconditional quantile estimates.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Akram Alyass and Arkan Abadi
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
