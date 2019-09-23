
# s.curve.update - Calculates output & incorporates perm test data ---------------

#' Calculates output & incorporates overall results, including permutation test data if present
#'
#' @param s.curve.mod estimated s.curve model from curveRunner
#' @param perm.dat permutation test data to add (if run externally)
#' @param perm.pvalues logical, calculate p-values against permutation test data? (if NULL, uses default in model object)
#'
#' @return updated s-curve model with results summaries (also printed to console)
#' @export
#'
#' @examples
s.curve.update <- function(s.curve.mod, 
                           perm.dat = NULL, # For combining permutation data in (if run externally)
                           perm.pvalues = NULL # Calculate model-derived p-values
){
  # (STEP 5 from readme.md) -------------
  
  
  # .(5a) Import perm test data ---------------------
  # Used to incorporate externally run permutation tests (e.g. on batching server) 
  if(!is.null(perm.dat)){
    s.curve.mod$perm.test <- perm.dat
    s.curve.mod$permutations <- length(perm.dat)
  }
  if(is.null(perm.pvalues)){
    perm.pvalues <- s.curve.mod$perm.pvalues
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
  if(!is.null(s.curve.mod$perm.test)){
    
    # Get proportion significance across all permuted models:
    s.curve.mod$perm.prop.sig <-
      vapply(
        s.curve.mod$perm.test,
        function(x){mean(x$significant)},
        1
      )
    
    # Mean significance rate across all permutations:
    s.curve.mod$perm.prop.sig.mean <-
      mean(s.curve.mod$perm.prop.sig)
    s.curve.mod$perm.prop.sig.median <-
      median(s.curve.mod$perm.prop.sig)
    
    
    # How many are more significant than the main set:
    s.curve.mod$perm.sig.compare <-
      s.curve.mod$prop.sig >= s.curve.mod$perm.prop.sig
    s.curve.mod$perm.sig.compare.proportion <-
      mean(s.curve.mod$perm.sig.compare)
    
    
    # MEDIAN EFFECT SIZE COMPARISONS ----
    s.curve.mod$perm.median.estimate <-
      vapply(
        s.curve.mod$perm.test,
        function(x){median(x$estimate)},
        1
      )
    
    # Mean significance rate across all models:
    s.curve.mod$perm.median.estimate.mean <-
      mean(s.curve.mod$perm.median.estimate)
    s.curve.mod$perm.median.estimate.median <-
      median(s.curve.mod$perm.median.estimate)
    
    s.curve.mod$perm.median.compare <-
      s.curve.mod$median.estimate >= s.curve.mod$perm.median.estimate
    s.curve.mod$perm.median.compare.proportion <-
      mean(s.curve.mod$perm.median.compare)
    
    # WITHIN SPECIFICATION COMPARISONS (not used in this paper) ----
    
    if(perm.pvalues){
      # Obtain matrix of estimates by specification
      s.curve.mod$estimate.matrix <-
        sapply(
          s.curve.mod$perm.test,
          function(x){x$estimate},
          simplify = "matrix"
        )
      
      # Names of quantiles:
      perm.quantile.names <-
        c("perm.lower", "perm.median", "perm.upper")
      
      # Default quantiles:
      perm.quantiles <-
        c(s.curve.mod$alpha/2, .5, 1-(s.curve.mod$alpha/2))
      
      if(!is.null(s.curve.mod$tail)){
        
        if(s.curve.mod$tail == "upper"){
          
          # Median and upper quantile only, taking full alpha:
          perm.quantiles <-
            c(.5, 1-s.curve.mod$alpha)
          
          perm.quantile.names <-
            c("perm.median", "perm.upper")
          
        }  
        
        if (s.curve.mod$tail == "lower") {
          
          # Median and lower quantile only, taking full alpha:
          perm.quantiles <-
            c(s.curve.mod$alpha, .5)
          
          perm.quantile.names <-
            c("perm.lower", "perm.median")
        }
      }
      
      # Obtain quantiles of estimate from each speficiation:
      s.curve.mod$results[perm.quantile.names] <-
        t(
          apply(
            s.curve.mod$estimate.matrix,
            1,
            quantile,
            probs = perm.quantiles
          ))
      
      
      # Obtain P-value from exact test approximation:
      # (comparing absolute values for two-tailed based on CLT,
      # so as to make it direction-agnostic.  Could do better, probably).
      
      if(is.null(s.curve.mod$tail)){
        
        s.curve.mod$results$perm.p.value <-
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
        
        s.curve.mod$results$perm.p.value <-
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
        
        s.curve.mod$results$perm.p.value <-
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
      s.curve.mod$results$perm.significant <-
        s.curve.mod$results$perm.p.value <= s.curve.mod$alpha
      
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
  if(!is.null(s.curve.mod$perm.test) && perm.pvalues) {
    s.curve.mod$ordering$perm.p.value <-
      order(s.curve.mod$results$perm.p.value)
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
  if(!is.null(s.curve.mod$perm.test)){
    s.curve.mod$perm.report.significance  <-   paste0(
      "\n",
      "Across ", s.curve.mod$permutations, " permutations, ",
      "the main p-curve analysis had the same or more significant\n",
      "results than ",
      round(100*s.curve.mod$perm.sig.compare.proportion, 2), "% ",
      "of permuted p-curve sets.\n",
      "Mean significance rate: ",
      round(100*s.curve.mod$perm.prop.sig.mean, 2), "%\n",
      "Median significance rate:",
      round(100*s.curve.mod$perm.prop.sig.median, 2), "%\n"
    )
    
    s.curve.mod$perm.report.median  <-   paste0(
      "\n",
      "Across ", s.curve.mod$permutations, " permutations, ",
      "the main p-curve analysis had an equivalent or larger median \n",
      "effect size than ",
      round(100*s.curve.mod$perm.median.compare.proportion, 2), "% ",
      "of permuted p-curve sets.\n",
      "Mean of median estimate: ",
      round(s.curve.mod$perm.median.estimate.mean, 3), "\n",
      "Median of median estimates: ",
      round(s.curve.mod$perm.median.estimate.median, 3), "\n"
    )
    
    cat(s.curve.mod$perm.report.significance)
    cat(s.curve.mod$perm.report.median)
    
  }
  
  # FINAL OUTPUT ----
  return(s.curve.mod)
  
}


