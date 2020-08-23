\name{Levene}
\alias{Levene}
\title{Levene Test of Variance Heterogeneity}
\description{Returns test statistics for Levene's test of variance heterogeneity across genotype categories.}
\usage{Levene(datatable, y, g, covariates=NULL, method="Brown-Forsythe",
                   log.p=FALSE)}
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{datatable}{An object of class `data.table` containing the data to be analysed with samples in rows and the dependent variable and all independent variables in columns. Ensure all variables are of the correct class (e.g. `numeric`, `factor`)}
  \item{y}{Name (character or string) of the continuous response variable}
  \item{g}{Name (character or string) of the genotype variable (e.g. SNP)}
  \item{covariates}{Name (character or string) of any additional variables (covariates) to fit the model. NOTE : here linear model will be fit with covariates and levene test will be performed on residuals from the linear model}
  \item{method}{Two options are implemented, "Classic-Levene"  and "Brown-Forsythe" modification. The former is known to be prone to elevated type 1 error when response is skewed while the latter is robust to skew.
  \item{log.p}{return log of p-value?}}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A vector with 4 elements containg the F statistic, degrees of freedom 1 and 2 and the p-value.}
\references{
[1] Paré, G., Cook, N. R., Ridker, P. M., & Chasman, D. I. (2010). On the use of variance per genotype as a tool to identify quantitative trait interaction effects: a report from the Women's Genome Health Study. PLoS Genetics, 6(6), e1000981. http://doi.org/10.1371/journal.pgen.1000981
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
