
## cl - linear model clustering function ----------------------------------------

## Included for use so that the specification curve can cluster
## standard errors if needed obtained from:
## https://www.ne.su.se/polopoly_fs/1.216115.1426234213!/menu/standard/file/clustering1.pdf





#' cl - linear model clustering function
#'
#' @param dat data frame from which model was fitted
#' @param fm fitted model including varible to be clustered
#' @param cluster string, variable by which to cluster standard errors
#'
#' @return
#' @export
#' @importFrom sandwich sandwich
#' @importFrom lmtest coeftest
#'
#' @examples
cl <- function(
  dat,
  fm,
  cluster
) {
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  coeftest(fm, vcovCL)
}

