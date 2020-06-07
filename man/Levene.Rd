\name{levene}
\alias{levene}
\title{Levene Function}
\description{This function performs levene test of variance heterogeneity across genotype categories.}
\usage{levene(datatable,y,g,covariates=NULL,method="Brown-Forsythe",
                   log.p=FALSE)}
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{datatable}{An object of class `data.table` containing the data to be analysed with samples in rows and the dependent variable and all independent variables in columns. Ensure all variables are of the correct class (e.g. `numeric`, `factor`)}
  \item{y}{Name (character or string) of the continuous response variable}
  \item{g}{Name (character or string) of the genotype variable (e.g. SNP)}
  \item{covariates}{Name (character or string) of any additional variables (covariates) to fit the model. NOTE : Additional covariates will impact the genotype results only in CQR. UQR additional covariates will have little impact. Default is NULL}
  \item{method}{Specify 'Classic-Levene' or 'Brown-Forsythe'}
  \item{log.p}{Specify whether to return log p-value.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A vector of containing the results of the levene test of variance heterogeneity of outcome across genotype categories.}
\references{
[1] Paré, G., Cook, N. R., Ridker, P. M., & Chasman, D. I. (2010). On the use of variance per genotype as a tool to identify quantitative trait interaction effects: a report from the Women's Genome Health Study. PLoS Genetics, 6(6), e1000981. http://doi.org/10.1371/journal.pgen.1000981
[2] Gauderman, W. J., Zhang, P., Morrison, J. L., & Lewinger, J. P. (2013). Finding novel genes by testing G × E interactions in a genome-wide association study. Genetic Epidemiology, 37(6), 603–613. http://doi.org/10.1002/gepi.21748
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
