# Adjusted p-values:
adjusted_pvals <- function(results, resample.test){

  b_boot_ests <-
    resample.test %>%
    map("estimate") # Pull out estimates

  mean_b_boot <-
    b_boot_ests %>%
    map_dbl(mean, na.rm=T) %>%  # Take the mean
    mean(na.rm=T) # Take the mean of means

  b_boot_biasadj_abs <-
    b_boot_ests %>%
    map(~abs(.x - mean_b_boot))

  adjusted_pvals <-
    b_boot_biasadj_abs %>%
    pmap(c) %>%
    map2_dbl(results$estimate, ~mean(.x >= abs(.y), na.rm=T))

  return(adjusted_pvals)
}

# Stouffer aggregator:

stouffer <- function(p) {
  stouffers_Z <-
    qnorm(p, lower.tail = FALSE) %>% # Gets normal deviate associated w/P-value
    pmin(4) %>% # If Z > 4, take 4 (upper bounding)
    pmax(-4) %>% # If Z < 4 take -4 (lower bounding)
    sum(na.rm=T) %>% # Sum of Z-values
    magrittr::divide_by(sqrt(sum(!is.na(p))))# Divide by square root of number of valid cases in resample

  return(stouffers_Z)
}
