

# Density Plots -----------------------------------------------------------


#

#' Significance density plot for permutation testing:
#'
#' @param s.curve.mod fitted s-curve model
#' @param subtitle.add subtitle to add to main title
#'
#' @return Outputs a ggplot object that plots a density curve of significance on the
#' @export
#'
#' @examples
s.curve.sig.density.plot <-
  function(
    s.curve.mod,
    subtitle.add = NULL
  ){

    env <- new.env(parent = globalenv())

    env$perm.prop.sig <- s.curve.mod$perm.prop.sig
    env$prop.sig <- s.curve.mod$prop.sig
    env$subtitle.add <- subtitle.add
    ds.plot <-
      with(env, {
        ggplot(mapping = aes(x = perm.prop.sig)) +
          stat_density(geom = "line") +
          stat_density(geom = "area", fill = "blue", alpha = .5) +
          geom_vline(xintercept = prop.sig, color = "red", linetype = 4) +
          theme_minimal() +
          xlab("Proportion of Significant Models") +
          ggtitle(label = "Density of Simulated-Null Significance Rates (Permutation Test)",
                  subtitle = paste0(subtitle.add, " (Actual Data Significance Rate ",
                                    round(100*prop.sig, 1), "%)")
          )})
    return(ds.plot)
  }

# Treatment estimate significance density plot for perm testing:

s.curve.est.density.plot <-
  function(
    s.curve.mod,
    subtitle.add = NULL
  ){
    env <- new.env(parent = globalenv())

    env$median.estimate <- s.curve.mod$median.estimate
    env$perm.median.estimate <- s.curve.mod$perm.median.estimate
    env$subtitle.add <- subtitle.add

    de.plot <-
      with(env, {
        ggplot(mapping = aes(x = perm.median.estimate)) +
          stat_density(geom = "line") +
          stat_density(geom = "area", fill = "green", alpha = .4) +
          geom_vline(xintercept = median.estimate, color = "red", linetype = 4) +
          theme_minimal() +
          xlab("Median Effect Size") +
          ggtitle(label = "Density of Simulated-Null Median Effect Sizes (Permutation Test)",
                  subtitle = paste0(subtitle.add,
                                    " (Actual Data Median Effect Size ",
                                    round(median.estimate, 3), ")")
          )
      })
    return(de.plot)
  }



