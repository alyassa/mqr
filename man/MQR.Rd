\name{MQR}
\alias{MQR}
\title{
Meta-Quantile Regression
}
\description{
Returns estimates, standard errors, test statistic and p-value from a fitted Meta-Quantile Regression model (specified either conditional or unconditional).
}
\usage{
MQR(datatable,y,g,covariates=NULL,tau=seq(0.05, 0.95, by=0.05),
                mqr.method="UQR",boot.m="mcmb",boot.R=200,seed=31371,enable.dither=TRUE,
                Univariable=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datatable}{An object of class `data.table` containing the data to be analysed with samples in rows and the dependent variable and all independent variables in columns. Ensure all variables are of the correct class (e.g. `numeric`, `factor`)
  \item{y}{Name (character or string) of the continuous response variable}
  \item{g}{Name (character or string) of the genotype variable (e.g. SNP)}
  \item{covariates}{Name (character or string) of any additional variables (covariates) to fit the model. NOTE : Additional covariates will impact the genotype results only in CQR. UQR additional covariates will have little impact. Default is NULL}
  \item{tau}{quantiles percentiles at which to fit UQR and CQR models. Should be between 0 and 1. Should be at least 3 quantiles. Selecting too many quantiles may result in overfitting. Default is every 0.05 quantile from 0.05 to 0.95 (i.e. `seq(0.05,0.95,by=0.05)`).}
    \item{mqr.method}{Specify either "UQR" or "CQR".}
    \item{boot.m}{Bootstrapping method for CQR. Default is "mcmb". See `summary.rq.boot` details for additional options but use these with caution as they may increase false positive rate.}
    \item{boot.R}{Number of bootstrap replicates for CQR. Default is 200. This number has been found to empirically adequate and larger numbers increase compute time exponentially.}
    \item{seed}{Set seed to ensure reproducibility and proper bootsrap sampling.}
    \item{enable.dither}{UQR and CQR can sometimes fail to resolve quantiles when the response contains many numerically identical values. Dithering introduces small amount of noise to the response. This is similar to `gitter` except the noise is all right-sided. During testing dithering had minimal impact on estimates, improved computation speed, and decreased failiure rate. Default is `TRUE`}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector with 4 elements containg the MQR estimate, standard error, t-value and p-value for the effects of genotype g on response y.
}
\references{
[1] Koenker, R. W. (2005). Quantile Regression, Cambridge U. Press.
[2] Firpo, S., Fortin, N. M., & Lemieux, T. (2009). Unconditional quantile regressions. Econometrica. http://doi.org/10.3982/ECTA6822
}
\author{
Arkan Abadi and Akram Alyass
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
