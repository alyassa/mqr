\name{Scaling.Adjustor}
\alias{Scaling.Adjustor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjusting for Scaling
}
\description{
This functions fits a model to estimate scaling (i.e. correlation between variant quantile-wideining effects and marginal effects). It then computes the estimated quantile-widening effects from the marginal effect of each variant and finally deducts this estimated quantile-widening effect from the actual quantile-widening effect to produce scaling-adjusting quantile widening effects for each variant.
}
\usage{
Scaling.Adjustor(FML, data, weights, model, response, predictor, method = "Median_rq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FML}{Univariable model formula, with quantile-widening effect response (dependent variable) and marginal effects as the predictor (independent variable).}
  \item{data}{Model frame of class `data.table` with quantile-widening effects and marginal effects for columns and variants for rows.}
  \item{weights}{Optional name  (character or string) of column containing weights. e.g. standard errors of marginal effect estimates can be used for inverse-variance weighting.}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{response}{
%%     ~~Describe \code{response} here~~
}
  \item{predictor}{
%%     ~~Describe \code{predictor} here~~
}
\item{method}{Type of scaling model to fit. Default is 'Median_rq', which implements median regression. Other options include 'Iteratively_ReWeighted', which implements robust regression; and 'ols', which implements least squares regression.}
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
