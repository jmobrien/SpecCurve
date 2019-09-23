
# Plotting tool for s.curve output ----------------------------------------



#' Plotting tool for specification curve output from s.curve
#'
#' @param s.curve.mod fitted model object from function s.curve() or
#'   s.curve.update()
#' @param what character, for each specification, what result is to be plotted?
#'   Accepts "estimate" (treatment estimate, default), "p.value", "perm.p.value"
#'   (permutation-test-estimated p-value), "AIC", "BIC"
#' @param title character, what should the title of the plot be?
#' @param plot.order character - what, if anything, should the values be sorted
#'   on? Accepts same values as estimate "parameter"
#' @param decreasing logical, when sorting, use decreasing order (default TRUE)
#' @param show.sig logical, should the points be colorized by significance?
#'   (default TRUE)
#' @param include.permutations logical, if permutations are present, include CI
#'   boundaries of permutations? (default FALSE)
#' @param perm.sig logical, use permutation p-values for significance markings?
#'   (default FALSE)
#' @param pointsize numeric, size of points to use in plotting (default 2)
#' @param x.alt numeric vector, with same length as # specifications.  If
#'   included, manually positions specifications
#' @param plot.theme character, name of theme available to ggplot to override
#'   default (e.g. "minimal", "bw", "classic")
#'
#' @return a ggplot plot summarizing the results of the fitted s.curve object
#'   passed to it.  Can be saved, or will render plot in viewer if printed to
#'   console.
#' @export
#'
#' @examples
s.curve.plot <-
  function(
    s.curve.mod,
    what = "estimate",
    title = NULL,
    plot.order = NULL,
    decreasing = TRUE,
    show.sig = TRUE,
    include.permutations = FALSE,
    perm.sig = FALSE,
    pointsize = 2,
    x.alt = NULL,
    plot.theme = NULL
  ){

    order.name <- setNames(c("Estimated Unstandardized Effect Size",
                             "P-value",
                             "Permutation Test Calculated P-value",
                             "AIC", "BIC"),
                           c("estimate", "p.value", "perm.p.value",
                             "AIC", "BIC"))
    ## Create a new plotting environment (needed so the plot objects don't
    ## hold all the extra perm test data)
    env <- new.env(parent = globalenv())

    ## Pass the function parameters to the new environment
    ## (avoids file problems when ggplot object defined in global env against large s-curve):
    env$what <- what
    env$order.name <- order.name
    env$title <- title
    env$plot.order <- plot.order
    env$decreasing <- decreasing
    env$show.sig <- show.sig
    env$include.permutations <- include.permutations
    env$perm.sig <- perm.sig
    env$tail <- s.curve.mod$tail
    env$pointsize <- pointsize
    env$plot.theme <- plot.theme
    env$y.label <- order.name[what]

    ## Sort by the provided metric, if needed:
    if(!is.null(plot.order)){
      new.order <- order(s.curve.mod$results[plot.order], decreasing = decreasing)

      env$plot.dat <-  s.curve.mod$results[new.order,]
      env$x.label <- paste0("Specification, Sorted by ", order.name[plot.order])
    } else {

      env$plot.dat <- s.curve.mod$results
      env$x.label <- "Specification Number"

    }

    ## Apply any alternative x-values (used for visibility positional
    ## adjustments):

    if(is.null(x.alt) || all(is.na(x.alt))){

      env$plot.dat$xval <- seq_len(nrow(env$plot.dat))

    } else {

      env$plot.dat$xval <- x.alt

    }


    # Make plot labels/params  -------------------------------------------------------

    ## Significance legend title, adjusted if perm-test p-values used::
    env$sig.title <-
      "Significance:"

    ## Using permutation-test based p-values?
    if(perm.sig && !is.null(s.curve.mod$permutations)){
      env$sig.title <-
        "Permutation-test\nestimated\nsignificance:"
    }

    ## P-value marker
    env$sig.labels <- c(
      paste0("p ≥ ", substring(s.curve.mod$alpha, 2)),
      paste0("p < ", substring(s.curve.mod$alpha, 2))
      )

    ## One tailed, if necessary:
    if(s.curve.mod$tail %in% c("upper", "lower")){
      env$sig.labels <-
        paste0(env$sig.labels, ", one-tailed")
    }


    ## Add in permutation testing results labels:
    if (!is.null(s.curve.mod$permutations) &&
        what == "estimate" &&
        include.permutations){
      env$CI.bound <- ifelse(is.null(s.curve.mod$tail),
                             paste0(100*(1 - s.curve.mod$alpha/2), "% CI"),
                             ifelse(s.curve.mod$tail == "upper",
                                    paste0("Upper ", 100*(1 - s.curve.mod$alpha), "% CI"),
                                    paste0("Lower ", 100*(1 - s.curve.mod$alpha), "% CI")
                             ))
    }

    ## Apply selected theme or select minimal:
    if (!is.null(plot.theme)){
      env$theme_choice <-
        match.fun(paste0("theme_", plot.theme))
    } else {
      env$theme_choice <- theme_minimal
    }


    ## Create Plot -------------------------------------------------------------

    newplot <-
      with(env, {

        ## Get the data for the y-axis
        plot.dat$yval <-
          plot.dat[[what]]

        ## Create a plot:
        new.plot <-
          ggplot(plot.dat,
                 aes(x=xval,
                     y=yval)) +
          theme_choice() +
          ylab(y.label) +
          xlab(x.label) +
          ggtitle(title)

        ## Colorize if significance colors desired:
        if(show.sig){

          if (perm.sig  && !is.null(permutations)){
            ## Using permutation-test based p-values:
            new.plot <-
              new.plot +
              geom_point(aes(color=perm.significant), size = pointsize) +
              scale_color_manual(
                values = c(`FALSE` = "black", `TRUE` = "red"),
                labels = sig.labels,
                name = sig.title)

          } else {
            ## Using standard model-estimated p-values:
            new.plot <-
              new.plot +
              geom_point(aes(color=significant), size = pointsize) +
              scale_color_manual(
                values = c(`FALSE` = "black", `TRUE` = "red"),
                labels = sig.labels,
                name = sig.title)
          }

        } else {

          ## Just make the plot w/o colors
          new.plot <-
            new.plot +
            geom_point(size = pointsize)

        }

        ## Draw lines for permutation test significance cutoffs:
        if ("perm.median" %in% names(plot.dat) &&
            what == "estimate" &&
            include.permutations){

          ## Median line:
          new.plot <-
            new.plot +
            geom_line(aes(y=perm.median, linetype = "Median")) +
            scale_linetype_manual(
              name = "Permutation\nTest Estimates",
              values =
                setNames(
                  c("solid", rep("dotted", length(CI.bound))),
                  c("Median", CI.bound)
                )
            )

          if(!is.na(tail)){
            ## Upper tail:
            if (tail == "upper"){
              new.plot <-
                new.plot +
                geom_line(aes(y=perm.upper, linetype = CI.bound))
            }

            ## Lower tail:
            if(tail == "lower"){
              new.plot <-
                new.plot +
                geom_line(aes(y=perm.lower, linetype = CI.bound))

            }

          } else {

            new.plot <-
              new.plot +
              geom_line(aes(y=perm.lower, linetype = CI.bound))
            geom_line(aes(y=perm.upper, linetype = CI.bound))
          }
        }
        ## output from the environment:
        return(new.plot)
      })

    ## Ouput the plot:
    return(newplot)


  }


