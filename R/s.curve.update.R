
# s.curve.update - Calculates output & incorporates perm test data ---------------

#' Calculates output & incorporates overall results, including permutation test data if present
#'
#' @param s.curve.mod estimated s.curve model from curveRunner
#' @param resample.dat permutation test data to add (if run externally)
#' @param resample.pvals logical, calculate p-values against permutation test data? (if NULL, uses default in model object)
#'
#' @return updated s-curve model with results summaries (also printed to console)
#' @export
#'
#' @examples
s.curve.update <- function(s.curve.mod,
                           resample.dat = NULL, # For combining permutation data in (if run externally)
                           resample.pvals = NULL # Calculate model-derived p-values
){
  # (STEP 5 from readme.md) -------------


  # .(5a) Import perm test data ---------------------
  # Used to incorporate externally run permutation tests (e.g. on batching server)
  if(!is.null(resample.dat)){
    s.curve.mod$resample.test <- resample.dat
    s.curve.mod$iterations <- length(resample.dat)
  }
  if(is.null(resample.pvals)){
    resample.pvals <- s.curve.mod$resample.pvals
  }



  # .(5b) Results ---------------------

  # .. Significance rates ----
  # Calculate number and percent significant across types:
  s.curve.mod$n.specifications <- nrow(s.curve.mod$results)
  s.curve.mod$n.sig <- sum(s.curve.mod$results$significant)
  s.curve.mod$prop.sig <- mean(s.curve.mod$results$significant)

  # Calculate mean estimate:
  s.curve.mod$mean.estimate <- mean(s.curve.mod$results$estimate)
  # Calculate median estimate:
  s.curve.mod$median.estimate <- median(s.curve.mod$results$estimate)


  # .. Permutation test rates ----

  # Significance rates for permutation test code
  if(!is.null(s.curve.mod$resample.test)){

    # Get proportion significance across all permuted models:
    s.curve.mod$resample.prop.sig <-
      vapply(
        s.curve.mod$resample.test,
        function(x){mean(x$significant)},
        1
      )

    # Mean significance rate across all permutations:
    s.curve.mod$resample.prop.sig.mean <-
      mean(s.curve.mod$resample.prop.sig)
    s.curve.mod$resample.prop.sig.median <-
      median(s.curve.mod$resample.prop.sig)


    # How many are more significant than the main set:
    s.curve.mod$resample.sig.compare <-
      s.curve.mod$prop.sig >= s.curve.mod$resample.prop.sig
    s.curve.mod$resample.sig.compare.proportion <-
      mean(s.curve.mod$resample.sig.compare)


    # MEDIAN EFFECT SIZE COMPARISONS ----
    s.curve.mod$resample.median.estimate <-
      vapply(
        s.curve.mod$resample.test,
        function(x){median(x$estimate)},
        1
      )

    # Mean significance rate across all models:
    s.curve.mod$resample.median.estimate.mean <-
      mean(s.curve.mod$resample.median.estimate)
    s.curve.mod$resample.median.estimate.median <-
      median(s.curve.mod$resample.median.estimate)

    s.curve.mod$resample.median.compare <-
      s.curve.mod$median.estimate >= s.curve.mod$resample.median.estimate
    s.curve.mod$resample.median.compare.proportion <-
      mean(s.curve.mod$resample.median.compare)

    # WITHIN SPECIFICATION COMPARISONS (not used in this paper) ----

    if(resample.pvals){
      # Obtain matrix of estimates by specification
      s.curve.mod$estimate.matrix <-
        sapply(
          s.curve.mod$resample.test,
          function(x){x$estimate},
          simplify = "matrix"
        )

      # Names of quantiles:
      resample.quantile.names <-
        c("resample.lower", "resample.median", "resample.upper")

      # Default quantiles:
      resample.quantiles <-
        c(s.curve.mod$alpha/2, .5, 1-(s.curve.mod$alpha/2))

      if(!is.null(s.curve.mod$tail)){

        if(s.curve.mod$tail == "upper"){

          # Median and upper quantile only, taking full alpha:
          resample.quantiles <-
            c(.5, 1-s.curve.mod$alpha)

          resample.quantile.names <-
            c("resample.median", "resample.upper")

        }

        if (s.curve.mod$tail == "lower") {

          # Median and lower quantile only, taking full alpha:
          resample.quantiles <-
            c(s.curve.mod$alpha, .5)

          resample.quantile.names <-
            c("resample.lower", "resample.median")
        }
      }

      # Obtain quantiles of estimate from each speficiation:
      s.curve.mod$results[resample.quantile.names] <-
        t(
          apply(
            s.curve.mod$estimate.matrix,
            1,
            quantile,
            probs = resample.quantiles
          ))


      # Obtain P-value from exact test approximation:
      # (comparing absolute values for two-tailed based on CLT,
      # so as to make it direction-agnostic.  Could do better, probably).

      if(is.null(s.curve.mod$tail)){

        s.curve.mod$results$resample.p.value <-
          vapply(
            seq_len(nrow(s.curve.mod$results)),
            function(n){

              # Get p-value from empirical cumulative distribution function:
              distfunc <- ecdf(s.curve.mod$estimate.matrix[n,])

              # Calculate two-tailed p-value
              2*(.5-(abs(distfunc(s.curve.mod$results$estimate[n])-.5)))

            }, 1)

      }

      if (s.curve.mod$tail == "upper") {

        s.curve.mod$results$resample.p.value <-
          vapply(
            seq_len(nrow(s.curve.mod$results)),
            function(n){

              # Get p-value from empirical cumulative distribution function:
              distfunc <- ecdf(s.curve.mod$estimate.matrix[n,])

              # Calculate upper p-value
              1-distfunc(s.curve.mod$results$estimate[n])

            }, 1)

      }

      if (s.curve.mod$tail == "lower") {

        s.curve.mod$results$resample.p.value <-
          vapply(
            seq_len(nrow(s.curve.mod$results)),
            function(n){

              # Get p-value from empirical cumulative distribution function:
              distfunc <- ecdf(s.curve.mod$estimate.matrix[n,])

              # Calculate lower p-value
              distfunc(s.curve.mod$results$estimate[n])

            }, 1)

      }

      # Critical value test, boolean:
      s.curve.mod$results$resample.significant <-
        s.curve.mod$results$resample.p.value <= s.curve.mod$alpha

    }
  }  # End permutation test code


  # ORDERING ----

  # Ordering of value for use in plotting:
  s.curve.mod$ordering <-
    data.frame(
      estimate = order(s.curve.mod$results$estimate),
      p.value = order(s.curve.mod$results$p.value),
      AIC = order(s.curve.mod$results$AIC),
      BIC = order(s.curve.mod$results$BIC)
    )

  # Ordering by permutation test p-values if available:
  if(!is.null(s.curve.mod$resample.test) && resample.pvals) {
    s.curve.mod$ordering$resample.p.value <-
      order(s.curve.mod$results$resample.p.value)
  }


  # REPORTING ----

  # print summary to screen for easy reading:
  s.curve.mod$report <-
    paste0(
      "Percentage of significant p-values in treatment term", "\n",
      "across all models is ",
      round(100*s.curve.mod$prop.sig, 2), "% (",
      s.curve.mod$n.sig, " of ", s.curve.mod$n.specifications, " specifications)\n")

  cat(s.curve.mod$report)


  # Print permutation summary to screen for easy reading:
  if(!is.null(s.curve.mod$resample.test)){
    s.curve.mod$resample.report.significance  <-   paste0(
      "\n",
      "Across ", s.curve.mod$iterations, " iterations, ",
      "the main s-curve analysis had the same or more significant\n",
      "results than ",
      round(100*s.curve.mod$resample.sig.compare.proportion, 2), "% ",
      "of s-curves with resampled treatment variables.\n",
      "Mean significance rate: ",
      round(100*s.curve.mod$resample.prop.sig.mean, 2), "%\n",
      "Median significance rate:",
      round(100*s.curve.mod$resample.prop.sig.median, 2), "%\n"
    )

    s.curve.mod$resample.report.median  <-   paste0(
      "\n",
      "Across ", s.curve.mod$iterations, " iterations, ",
      "the main s-curve analysis had an equivalent or larger median \n",
      "effect size than ",
      round(100*s.curve.mod$resample.median.compare.proportion, 2), "% ",
      "of s-curves with resampled treatment variables.\n",
      "Mean of median estimate: ",
      round(s.curve.mod$resample.median.estimate.mean, 3), "\n",
      "Median of median estimates: ",
      round(s.curve.mod$resample.median.estimate.median, 3), "\n"
    )

    cat(s.curve.mod$resample.report.significance)
    cat(s.curve.mod$resample.report.median)

  }

  # FINAL OUTPUT ----
  return(s.curve.mod)

}


