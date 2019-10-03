comb <- function(x, ...) {
  #rbind comb for data.table
  mapply(function(...) rbind(...,fill=TRUE),x,...,SIMPLIFY=FALSE)
}
