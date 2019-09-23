
# Tablemaking tool for s.curve output -------------------------------------



#' Creates a summary table of specification curve results (for export to spreadsheet)
#'
#' @param s.curve.mod fitted specification curve model
#' @param ordering character, choice of orderings from the model.  Accepts
#'   "estimate", "p.value", "perm.p.value", "AIC", "BIC"
#'
#' @return dataframe containing the treatment parameter estimates, statistics,
#'   and permutation test results for each specification
#' @export
#'
#' @examples
s.curve.table <- 
  function(
    s.curve.mod, 
    ordering = NULL
    ){
  
  table.displaylist <-
    intersect(
      c("formulas", "subset.used", "weights.used", 
        "estimate", "std.error", "p.value",
        "perm.lower", "perm.median", "perm.upper"),
      names(s.curve.mod$results)
    )
  
  
  order.touse <- 
    if(!is.null(ordering)){
      order(s.curve.mod$results[[ordering]], decreasing = TRUE)
    } else { 
      TRUE 
    }
  
  s.table <- 
    data.frame(model = 1:nrow(s.curve.mod$results))
  
  s.table[table.displaylist] <-
    s.curve.mod$results[order.touse, table.displaylist]
  
  return(s.table)
}


