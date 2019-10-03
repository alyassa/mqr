\name{Scaling.Adjustor}
\alias{Scaling.Adjustor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Scaling.Adjustor(FML, data, weights, model, response, predictor, method = "Iteratively_ReWeighted")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FML}{
%%     ~~Describe \code{FML} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{response}{
%%     ~~Describe \code{response} here~~
}
  \item{predictor}{
%%     ~~Describe \code{predictor} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (FML, data, weights, model, response, predictor, method = "Iteratively_ReWeighted") 
{
    DT <- copy(data)
    if (missing(model)) {
        DT[, `:=`(w, 1/(DT[[weights]]^2))]
        if (method == "Iteratively_ReWeighted") {
            Model <- rlm(FML, data = DT, weights = w, method = "M", 
                wt.method = "inv.var")
        }
        else if (method == "Median_rq") {
            Model <- rq(FML, tau = 0.5, data = DT, weights = w, 
                method = "fn")
        }
        else if (method == "ols") {
            Model <- lm(FML, data = DT)
        }
        y.name <- colnames(Model$model)[1]
        y.se <- gsub("Beta", "SE", y.name)
        x.name <- colnames(Model$model)[2]
        x.se <- gsub("Beta", "SE", x.name)
        Pred <- Model$model[, x.name] * Model$coefficients[[2]] + 
            Model$coefficients[[1]]
        Pred.se <- DT[[x.se]] * Model$coefficients[[2]] + Model$coefficients[[1]]
        Beta <- Model$residuals
        se <- sqrt(DT[[y.se]]^2 + Pred.se^2)
        tval <- Beta/se
        p.value <- 2 * pnorm(abs(tval), lower = FALSE)
        Results <- cbind(Pred, Pred.se, Beta, se, tval, p.value)
        r.name <- gsub("Beta", "", y.name)
        colnames(Results) <- c(paste0("Pred.", r.name, "Beta"), 
            paste0("Pred.", r.name, "SE"), paste0("Adj.", r.name, 
                "Beta"), paste0("Adj.", r.name, "SE"), paste0("Adj.", 
                r.name, "tval"), paste0("Adj.", r.name, "p.value"))
        Results <- data.table(Results)
        return(list(Model, Results))
    }
    else {
        Model <- copy(model)
        x.se <- gsub("Beta", "SE", predictor)
        y.se <- gsub("Beta", "SE", response)
        Pred <- DT[[predictor]] * Model$coefficients[[2]] + Model$coefficients[[1]]
        Pred.se <- DT[[x.se]] * Model$coefficients[[2]] + Model$coefficients[[1]]
        Beta <- DT[[response]] - Pred
        se <- sqrt(DT[[y.se]]^2 + Pred.se^2)
        tval <- Beta/se
        p.value <- 2 * pnorm(abs(tval), lower = FALSE)
        Results <- cbind(Pred, Pred.se, Beta, se, tval, p.value)
        r.name <- gsub("Beta", "", response)
        colnames(Results) <- c(paste0("Pred.", r.name, "Beta"), 
            paste0("Pred.", r.name, "SE"), paste0("Adj.", r.name, 
                "Beta"), paste0("Adj.", r.name, "SE"), paste0("Adj.", 
                r.name, "tval"), paste0("Adj.", r.name, "p.value"))
        Results <- data.table(Results)
        return(Results)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
