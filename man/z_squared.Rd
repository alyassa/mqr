\name{z_squared}
\alias{z_squared}
\title{Z-Squared of Assessing Variance Heterogeneity}
\description{Returns test statistics for z-squared testing of variance heterogeneity.}
\usage{Levene(datatable, y, g, covariates=NULL)}
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{datatable}{An object of class `data.table` containing the data to be analysed with samples in rows and the dependent variable and all independent variables in columns. Ensure all variables are of the correct class (e.g. `numeric`, `factor`)}
  \item{y}{Name (character or string) of the continuous response variable}
  \item{g}{Name (character or string) of the genotype variable (e.g. SNP)}
  \item{covariates}{Name (character or string) of any additional variables (covariates) to fit the model. NOTE : here linear model will be fit with covariates and levene test will be performed on residuals from the linear model}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A vector with 4 elements containg the effect size estimate, standard error, Z statistic and p-value.}
\references{
[1] Yang, J., Loos, R. J. F., Powell, J. E., Medland, S. E., Speliotes, E. K., Chasman, D. I., et al. (2012). FTO genotype is associated with phenotypic variability of body mass index. Nature, 490(7419), 267–272. http://doi.org/10.1038/nature11401
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
