# Subsetting tool ---------------------------------------------------------



#' Utility function for subsetting specification curves based on specification characteristics
#'
#' @param s.curve.mod fitted specification curve model object
#' @param ... character, logical statements that can be used to filter the specification based on the results
#'
#' @return calls s.curve.update and returns results and output from that function, including permutation tests if appropriate
#' @export
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @examples
s.curve.subset <-
  function(s.curve.mod,
           ...
  ){


  # (STEP 6 from readme.md) FILTERING --------------------------
  # .(6a) ------------------------------------------------------
  # Filters a subset of the data based on criteria w/in the results:
  s.curve.mod$results <-
    eval(substitute(
      s.curve.mod$results %>%
        filter(...)
      ))

  # Filters the specification list equivalently, using the index:
  s.curve.mod$spec.list <-
    s.curve.mod$spec.list %>%
    filter(specification.index %in%
             s.curve.mod$results$specification.no)

  # Filters the permutation test equivalently:
  if(!is.null(s.curve.mod$perm.test)){
    s.curve.mod$perm.test <-
      s.curve.mod$perm.test %>%
      map(~.x %>%
            filter(specification.no %in%
                     s.curve.mod$results$specification.no)
      )
  }

  # .(6b) ------------------------------------------------------
  # Passes this back to the s.curve.update for reporting:
  s.curve.update(s.curve.mod)
}


