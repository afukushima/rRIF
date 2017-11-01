##' @export
print.rRIF <- function(x, ...) {
  print.data.frame(make.rRIFtable(x))
}
##' @export
summary.rRIF <- function(object, ...) {
  res <- make.rRIFtable(object)
  
  cat(paste("\n\n===== rRIF Summary =====\n"))
  cat(paste("Top 5 RIF1:\n"))
  utils::head(res[base::order(abs(res$RIF1), decreasing = TRUE),], 5)
}

make.rRIFtable <- function(x) {
  tgt <- names(x$RIF1)
  Ai <- x$Ai[tgt]
  dEi <- x$dEi[tgt]
  PIFi <- x$PIFi[tgt]
  RIFtable <- data.frame(cbind(Ai, dEi, PIFi, x$RIF1, x$RIF2))
  colnames(RIFtable) <- c("A", "DE", "PIF", "RIF1", "RIF2")
  return(RIFtable)
}
