\name{geno_generator}
\alias{geno_generator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotype Generator
}
\description{
Function to generate genotypes with set paramenters usefull in simulations. Importantly this functions allows generation of genotypes that deviate from hardy-weinberg equilibrium.
}
\usage{
geno_generator(N=10000,MAF=0.05,hwe.p=1,min.grp.size=5,Pare.Encoding=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Sample size, must be positive integer}
  \item{MAF}{Minor allele frequency}
  \item{hwe.p}{p-value of hardy weinberg equilibrium, distirbution of homozygotes (pp,qq) and heterozygotes (pq) works backwards from here according to the chi test of hardy weinberg equilibrium.}
  \item{Pare.Encoding}{Codes genotypes according to Pare 2010 (see references) to have an approximate mean of zero (-2p,1-2p,2-2p). If FALSE, genotypes are encoded by number of copies of minor allele (0,1,2).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of length N, containing sample genotypes.
}
\references{
Paré, G., Cook, N. R., Ridker, P. M., & Chasman, D. I. (2010). On the use of variance per genotype as a tool to identify quantitative trait interaction effects: a report from the Women's Genome Health Study. PLoS Genetics, 6(6), e1000981. http://doi.org/10.1371/journal.pgen.1000981
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
function (x, ...)
{
    mapply(function(...) rbind(..., fill = TRUE), x, ..., SIMPLIFY = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
