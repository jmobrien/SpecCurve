
# Reporting function ------------------------------------------------------



#' Title Display s.curve results to console
#'
#' @param s.curve a fitted specification curve object
#'
#' @return
#' @export
#'
#' @examples
s.curve.report <- function(s.curve){
  cat(s.curve$report)
  if(!is.null(s.curve$resample.report.significance)){
    cat(s.curve$resample.report.significance)
    cat(s.curve$resample.report.median)

  }
}
