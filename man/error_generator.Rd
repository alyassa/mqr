\name{error_generator}
\alias{error_generator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error Generator
}
\description{
Function to generate error distribution. Useful for generator various shpes (e.g. skewed).
}
\usage{
error_generator(N=10000,Type="Normal",a,b,Skew.Dir="Right")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Sample size, must be positive integer}
  \item{Type}{Specify the type of distribution as follows
  METHOD 0 : Normal (not skewed)
             Type="Normal"
  METHOD 1 : uses the Skew-N Distribution. This is pretty good but may not produce
             ideal shape of skewing. Increasing the term "alpha" term increases
             skewing. Dont bother with alpha < 15 because skewing is not
             substantial enough.
             Type="Skew-N"
  METHOD 2 : uses the chi distribution which is really effective, increasing "df" term
             brings it closer and closer to normal distribution. 2 < df < 10 is ideal
             for testing.
             Type="Chi"
  METHOD 3 : uses the beta distribution which allows independent control of right and
             left side of the distribution, "a" term specifies steepness on the left
             side of the distribution (increasing the term decreases the steepness)
             and "b" is the same except on the right side of the distribution.
             a=10,b=10 is pretty close-ish to normal (but not really). To skew,
             set a and b to lower or higher values. ie left skewed distribution a<<b
             try a=2,b=10 for only minor skewing ideally a=1 and 10<b<20
             SET a=2, b=2 FOR FAT DISTRIUBTION WITH SHORT TAILS
             Type="Beta"
  METHOD 4 : log-normal distribution. Not much to this, no room to modify skewness
             throwing it in there because BMI is often log transformed
             Type="LogNormal"
  METHOD 5 : Laplace distribution. Used to Specify SKINY DISTRIBUTION AND LONG TAILS
             Type="Laplace"
  METHOD 6 : Uniform distribution. Used to Specify FAT DISTRIBUTION AND SHORT TAILS
             Type="Uniform"
  METHOD 7 : Adding outliers to distribution. Specify left or right side or both. Specify
           : percent of samples to be made into outliers. a=0.02
             Type="Left-Outliers", "Right-Outliers", "Both-Outliers"
  METHOD 6 : Exponential distribution, needed for more extreme skewing. a <=1
             Type="Exp"}
  \item{a}{alpha parameter for many distribution types.}
  \item{b}{beta parameter for many distribution types.}
  \item{Skew.Dir}{Direction of skew. Can be "Right" or "Left".}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of length N, containing sample genotypes.
}
\references{
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
