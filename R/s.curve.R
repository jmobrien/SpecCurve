
#' Specification curve maker function
#'
#' @param dat a dataframe containing the variables
#' @param outcomes vector of outcome variables
#' @param treatment vector of treatment variables
#' @param cov.list named list, sets of covariates / moderators - write
#'   moderators as "var1:var2"; eg: cov.list = list(gender = c("gender.roster",
#'   "gender.selfreport")
#' @param no.cov.exclude Will add a model where each item of the list above is
#'   missing, unless specified here
#' @param ranef.list Random effects to cross, if model uses them (always included in such models)
#' @param extra.models Models written literally; will be appended to the set
#'   (still crossed against subsets and weights)
#' @param extra.treatment Treatment variable in extra models (length 1 or same
#'   as other, needed if doing inverse probability weighting
#' @param mod.type bare function or string of functions--takes lm, glm, lmer, glmer, glmmTMB
#' @param mod.family if family needed
#' @param alpha nominal alpha value (if one tailed, will test against p <= .1
#'   on that side)
#' @param tail One tailed test?  Takes "upper" and "lower"
#' @param subsets runs subsets of data based on the specifications listed here
#'   (vector of conditions)
#' @param subsets.exclude if subsets added, includes an un-subsetted version
#' @param weights vector of weight variables names to add to runs, re-named if
#'   desired or, a weighting method (currenly accepts "ipw.calc" for inverse
#'   probability weighting)
#' @param weights.exclude Includes an unweighted version where weights added
#' @param weights.ipw.vars in addition to treatment variable, any other (factor
#'   or categorical) variables to consider when (re)weighting subsets
#' @param permutations optional, number of permutations for p-curve
#' @param perm.pvalues logical, calculates permutation test-based pvalues (if
#'   permutation test active)
#' @param perm.seed optional, RNG seed to use for permutation test (integer,
#'   will create one if not present)
#' @param cluster optional, Cluster robust standard errors
#' @param cluster.var optional, Variable on which to cluster, string
#' @param robust.se optional, heteroskedasticity-consistent standard error
#'   adjustment provided by package "sandwich".  See details in ?vcovHC.
#' @param cat.percent logical, displays summary output of % significant at the
#'   end of the data for convenience in interactive mode. Set to false if using
#'   as a part of an Rmarkdown file.
#' @param keep.tidy.models logical, set to false for large model samples
#' @param keep.full.models logical, set to false for large model samples
#' @param model.only logical, outputs the model for later use, rather than
#'   running (default is FALSE)
#'
#' @return depending on paramters, either a fitted s-curve object or a template for fitting an s-curve object
#' @export
#'
#' @examples
s.curve <- function(dat, outcomes, treatment,
                    cov.list, no.cov.exclude = NULL,
                    ranef.list = NULL,
                    extra.models = NULL,
                    extra.treatment = NULL,
                    mod.type = lm, mod.family = NULL,
                    alpha = .05, tail = NULL,
                    subsets = NULL, subsets.exclude = TRUE,
                    weights = NULL, weights.exclude = TRUE,
                    weights.ipw.vars = NULL,
                    permutations = NULL,
                    perm.pvalues = FALSE, perm.seed = NULL,
                    cluster = FALSE, cluster.var = NULL,
                    robust.se = NULL,
                    cat.percent = TRUE,
                    keep.tidy.models = TRUE,
                    keep.full.models = FALSE,
                    model.only = FALSE) {

  ## Type of model:
  if(class(mod.type) == "function"){
    model.type <- deparse(substitute(mod.type))
  }  else {
    # If it's already a string, then copy it over:
    model.type <- mod.type
  }

  # Load packages, if necessary:
  if("lmer" %in% model.type) stopifnot(require(lmerTest))
  if("glmer" %in% model.type) stopifnot(require(lme4))
  if("glmmTMB" %in% model.type) stopifnot(require(glmmTMB))


  # Set variable indicating whether there are single-level models:
  has.ols.mods <-
  ifelse(any(c("lm", "glm") %in% model.type),
         TRUE,
         FALSE
  )


  # Check whether random effects are present:
  if(any(c("lmer", "glmer", "glmmTMB") %in% model.type)){

    if(is.null(ranef.list)){
      stop("multi-level models require random effects terms to be specified")
    }

    # Set variable indicating that there are models needing random effects:
    has.mlm.mods <- TRUE

  } else {
    # Set to indicate absence:
    has.mlm.mods <- FALSE

  }


  # Check whether generalized models are present (needs "family" argument(s))
  if(any(c("glm", "glmer", "glmmTMB") %in% model.type)){

    if(is.null(mod.family)){
      stop('generalized models require "mod.family" argument with one or more elements')
    }

    # Set variable indicating that there are models needing random effects:
    has.gen.mods <- TRUE

  }

  # If tail is null, set to "2tail"
  if(is.null(tail)) {
    tail <- "2tail"
    }

  ## INITIALIZATION ----
  s.curve.mod <-
    list(
      dat = dat,
      treatment = treatment,
      outcomes = outcomes,
      alpha = alpha,
      tail = tail,
      mod.type = model.type,
      mod.family = mod.family,
      cluster = cluster,
      cluster.var = cluster.var,
      robust.se = robust.se,
      permutations = permutations,
      perm.pvalues = perm.pvalues,
      perm.seed = perm.seed,
      weights.ipw.vars = weights.ipw.vars,
      keep.tidy.models = keep.tidy.models,
      keep.full.models = keep.full.models
    )

  ## (STEP 2 from readme.md) FORMULA BUILDING ----

  ## Create a vector of all included variables,
  ## outcomes, treatment variables, and covariates:
  vars.list <-
    unique(c(outcomes, treatment, unlist(cov.list), unlist(ranef.list)))

  ## .(2a & 2b) cross treatment, outcomes, & covariates ----

  ## Add NA's to beginning of each category to model exclusion of that
  ## set from specification, excepting any categories marked for
  ## constant inclusion:
  cov.list.na <-
    sapply(
      names(cov.list),
      function(x){
        if(x %in% no.cov.exclude){
          cov.list[[x]]
        } else {
          c(NA, cov.list[[x]])
        }
      },
      USE.NAMES = TRUE,
      simplify = FALSE
    )


  ## Add in condition variables for right-hand side
  var.list.rhs <-
    c(list(treatment = treatment), cov.list.na)

  ## Add in outcome variables
  var.list.all <-
    c(list(outcomes = outcomes), var.list.rhs)


  ## initialize formula vectors
  formulas.ols <- c()
  formulas.mlm <- c()

  if(has.ols.mods){

    ## Expand to a data frame where each case
    ## elements to go in a formula on right-hand side
    mod.grid.ols <-
      expand.grid(var.list.rhs,
                  stringsAsFactors = FALSE)

    ## elements to go in all of them:
    mod.grid.all.ols <-
      expand.grid(var.list.all,
                  stringsAsFactors = FALSE)

    ## Make right-hand side formula
    formulas.rhs.ols <-
      apply(mod.grid.ols, 1, function(x){
        ## dropping NA's (when category excluded):
        formula.rhs <- paste(na.omit(x), collapse= " + ")
        ## adding "1" to empty formula (to avoid errors):
        formula.rhs[formula.rhs == ""] <- 1
        formula.rhs
      })

    ## Probably don't want empty models, so throw warning:
    if( any(formulas.rhs.ols == "1") ) {
      warning(
        "NOTE: set includes one or more null models."
      )}


    ## Adds left-hand side (outcomes), make full formulas:
    formulas.ols <-
      lapply(
        as.vector(outer(outcomes, formulas.rhs.ols,
                        paste, sep=" ~ ")),
        as.formula)

  }


  if(has.mlm.mods){

    ## Add in condition variables for right-hand side
    var.list.rhs.mlm <-
      c(var.list.rhs, ranef.list)

    ## Add in outcome variables
    var.list.all.mlm <-
      c(list(outcomes = outcomes), var.list.rhs.mlm)

    ## Expand to a data frame where each case
    ## elements to go in a formula on right-hand side
    mod.grid.mlm <-
      expand.grid(var.list.rhs.mlm,
                  stringsAsFactors = FALSE)

    ## elements to go in all of them:
    mod.grid.all.mlm <-
      expand.grid(var.list.all.mlm,
                  stringsAsFactors = FALSE)



    ## .(2c) construct formulas ----

    ## Make right-hand side formula
    formulas.rhs.mlm <-
      apply(mod.grid.mlm, 1, function(x){
        ## dropping NA's (when category excluded):
        formula.rhs <- paste(na.omit(x), collapse= " + ")
        ## adding "1" to empty formula (to avoid errors):
        formula.rhs[formula.rhs == ""] <- 1
        formula.rhs
      })

    ## Probably don't want empty models, so throw warning:
    if( any(formulas.rhs.mlm == "1") ) {
      warning(
        "NOTE: set includes one or more null models."
      )}


    ## Adds left-hand side (outcomes), make full formulas:
    formulas.mlm <-
      lapply(
        as.vector(outer(outcomes, formulas.rhs.mlm,
                        paste, sep=" ~ ")),
        as.formula)

  }

  # Data frame of formulas:
  s.curve.mod$formulas.df <-
    data.frame(
        type = c(
          rep("single-level", length(formulas.ols)),
          rep("multi-level", length(formulas.mlm))
          ),
        # I() for taking formula list as-is:
        formulas = I(c(formulas.ols, formulas.mlm))
          )

  s.curve.mod$formulas <-
    s.curve.mod$formulas.df$formulas

  # ## Construct models specifications that incorporate single/multilevel and family:
  # spec.list <-
  #
  # if("lm" %in% model.family){
  #     expand.grid(
  #       mod.type = "lm",
  #       formula = formulas.ols,
  #       family = NA
  #     )
  # }


  ## .(2d) extra models ----

  ## Have a vector of treatment variables for use later:
  # treatment.vars <-
  #   mod.grid.all$treatment

  ## Add in any extra models, if applicable:
  if(!is.null(extra.models)){

    s.curve.mod$formulas <-
      c(s.curve.mod$formulas,
        lapply(extra.models, as.formula)
        )


  ## Converts formula names to character vector for easier review later:
  s.curve.mod$formula.names <-
    sapply(s.curve.mod$formulas, deparse, width.cutoff = 500)

    # if (is.null(extra.treatment)){
    #   if (length(treatment) == 1){
    #     ## fill it in with the presumed treatment variable
    #     ## (to be updates later for more robust checking):
    #     treatment.vars <-
    #       c(treatment.vars, rep(treatment, length(extra.models)))
    #
    #   } else {
    #     stop("need to specify treatment variables for extra models if more than 1")
    #   }
    # } else {
    #   treatment.vars <- c(treatment.vars, extra.treament)
    # }
  }



  ## Length of formula list (# initial specifications):
  s.curve.mod$n.formulas <-
    length(s.curve.mod$formulas)

  ## Initialize names of everything:
  full.names <-
    s.curve.mod$formula.names



  ## .(2e) subsets & weightings ----

  ## ..........subsets ----
  ## Make the data subsets:
  if (!is.null(subsets)){
    subsettings <-
      do.call(
        cbind.data.frame,
        lapply(subsets, function(singlesubset){
          with(dat, eval(parse(text = singlesubset)))
        }))

    ## Add prefix to subset for identification:
    subset.names <-
      paste0("SUBSET: ", subsets)

    names(subsettings) <-
      subset.names

    ## Add a no-subset "subset" if that's relevant
    if (subsets.exclude) {
      subsettings <-
        cbind.data.frame(
          setNames(data.frame(TRUE), "SUBSET: ALL DATA"),
          subsettings
        )

      subset.names <-
        names(subsettings)
    }


    ## Update the full names
    full.names <-
      as.vector(outer(full.names, subset.names,
                      paste, sep=", "))

    ## Write to model for curveRunner use:
    s.curve.mod$n.subsettings <-
      length(subset.names)
    s.curve.mod$subsettings <-
      subsettings
    s.curve.mod$subset.names <-
      subset.names

  } else {

    ## Just have a non-subsetting subset:
    subsettings <- setNames(NA, NA)

  }


  ## ..........weights ----


  ## Make the data weights:
  if (!is.null(weights)){

    ## Version if single weights vector is provided:
    if (is.numeric(weights)){

      ## Throw error if weights is incorrect length:
      if (length(weights) != nrow(dat)) stop("Weighting variable must be length of data")

      ## name weights specification component:
      weights.names <-
        "WEIGHT: WEIGHTED"

      ## Make a data frame of weights
      weightings <-
        data.frame(`WEIGHT: WEIGHTED` = weights)

    } ## End numeric approach

    ## Version if weights are provided in list vector form:
    if (is.list(weights)){

      ## Throw error if weights is incorrect length or not numeric:
      if (
        any(vapply(weights, length, 1) != nrow(dat)) ||
        any(vapply(weights, class, "class") != "numeric")
      ) stop("List of weights must all be numeric vectors equal to length of data.")


      ## Make a vector of names for the weights
      if(!is.null(names(weights))){
        ## use provided names if available:
        weights.names <- paste0("WEIGHT: ", names(weights))
      } else {
        ## otherwise make a simple sequence:
        weights.names <- paste0("WEIGHT: WEIGHT", seq_along(weights))
      }

      ## Make a data frame of weights
      weightings <-
        as.data.frame(weights)

      ## apply the names:
      names(weightings) <-
        weights.names

    } ## End list vector approach


    ## Version if names are provided (in character form):
    if (is.character(weights) && weights != "ipw.calc"){

      ## Make a vector of names for the weights
      if(!is.null(names(weights))){
        ## use provided names if available:
        weights.names <- paste0("WEIGHT: ", names(weights))
      } else {
        ## or just use the variable names from the dataframe:
        weights.names <- paste0("WEIGHT: ", weights)
      }

      ## Assign the weightings to a new dataframe:
      weightings <-
        do.call(
          cbind.data.frame,
          lapply(weights, function(wset){
            dat[wset]
          }))

      ## apply the names:
      names(weightings) <- weights.names

    } ## End character variable approach

    ## Version if "ipw.calc" is called:
    if (weights == "ipw.calc") {


      weights.names <- paste0("WEIGHT: IPW (calc)")

      weightings <-
        setNames(
          data.frame(rep(1, nrow(dat))),
          weights.names
        )
    } ## End dynamic calculation approach (actual reweights handled in
      ## curveRunner below)


    ## Add extra column if weights.exclude == TRUE
    if (weights.exclude) {
      weightings <-
        cbind.data.frame(
          setNames(data.frame(1), "WEIGHT: NONE"),
          weightings
        )
      weights.names <- names(weightings)
    }

    ## Construct new specification names based on the weights:
    full.names <-
      as.vector(
        outer(full.names, weights.names,
              paste, sep=", ")
      )

    ## Write to model for curveRunner use:
    s.curve.mod$weightings <- weightings
    s.curve.mod$weights.names <- weights.names
    s.curve.mod$n.weightings <- length(weights.names)

  } else {

    ## If no weights, just make a NA variable to indicate:
    weightings <- setNames(NA, NA)

  }



  ## CREATE SPECIFICATION INDEX LIST ----

  ## total number of specifications after variables and weightings:
  # s.curve.mod$n.specifications <- length(full.names)


  # Initialize specification list:
  s.curve.mod$spec.list <- data.frame()

  # Binding function for specification list:
  specbind <- function(model.type, multilevel, generalized, formulas, family){
    rbind(s.curve.mod$spec.list,
          expand.grid(
            model.type = model.type,
            family = family,
            weights = names(weightings),
            subset = names(subsettings),
            multilevel = multilevel,
            generalized = generalized,
            formula = formulas,
            stringsAsFactors = FALSE
          )
    )
    }


  if("lm" %in% model.type){

  ## Add linear models to the specification list:
  s.curve.mod$spec.list <-
    specbind(model.type = "lm",
             multilevel = FALSE,
             generalized = FALSE,
             formulas = formulas.ols,
             family = NA
    )

  }

  if("glm" %in% model.type){

    ## Add linear models to the specification list:
    s.curve.mod$spec.list <-
      specbind(model.type = "glm",
               multilevel = FALSE,
               generalized = FALSE,
               formulas = formulas.ols,
               family = mod.family
      )

  }

  if("lmer" %in% model.type){

    ## Add linear models to the specification list:
    s.curve.mod$spec.list <-
      specbind(model.type = "lmer",
               multilevel = TRUE,
               generalized = FALSE,
               formulas = formulas.mlm,
               family = NA
      )

  }

  if("glmer" %in% model.type){

    ## Add linear models to the specification list:
    s.curve.mod$spec.list <-
      specbind(model.type = "glmer",
               multilevel = TRUE,
               generalized = TRUE,
               formulas = formulas.mlm,
               family = mod.family
      )

  }

  if("glmmTMB" %in% model.type){

    ## Add linear models to the specification list:
    s.curve.mod$spec.list <-
      specbind(model.type = "glmmTMB",
               multilevel = TRUE,
               generalized = TRUE,
               formulas = formulas.mlm,
               family = mod.family
      )

  }


  ## Add in a treatment variable for the index (used for inverse probability weighting)
  s.curve.mod$spec.list$treatment <- treatment

  ## Add in count of specifications:
  s.curve.mod$n.specifiations <-
    nrow(s.curve.mod$spec.list)

  ## Give each specification a numeric index:
  s.curve.mod$spec.list$specification.index <-
    seq_len(s.curve.mod$n.specifiations)


  ## Just output the s-curve overall specification, if requested:
  if(model.only) return(s.curve.mod)




  ## (STEP 3) RUN MAIN SPECIFICATION CURVE --------------------------------
  results.list <-
    curveRunner(s.curve.mod)

  s.curve.mod$results <-
    do.call(results.list, what = bind_rows)


  ## (STEP 4) POST-RUN PROCESSING ------------------------------------------

  ## .(4a) duplicate handling -------
  ## Check for duplicated models, note positions:
  duplicates <-
    which(
      duplicated(
        lapply(s.curve.mod$results$final.specification,
               ## Moves around terms so that duplicates aren't missed due to
               ## R-style formatting:
               function(x){
                 str_split(x, ":") %>%
                   lapply(sort) %>%
                   vapply(function(x){paste0(x, collapse = ":")}, FUN.VALUE = "char") %>%
                   sort
               })
      ))

  ## Set a duplicate flag:
  s.curve.mod$has.duplicates <- FALSE

  ## Make note of duplicates, if existing, and update:
  if(length(duplicates) > 0){

    s.curve.mod$has.duplicates <- TRUE

    ## Moves results of any duplicated models to separate location:
    s.curve.mod$duplicated.models <-
      s.curve.mod$results[duplicates,]

    ## Moves duplicate specifications from specification list:
    s.curve.mod$duplicate.specs <-
      s.curve.mod$spec.list[duplicates,]

    ## Drops duplicates from results:
    s.curve.mod$results <-
      s.curve.mod$results[-duplicates,]

    ## Drops duplicates from specification list:
    s.curve.mod$spec.list <-
      s.curve.mod$spec.list[-duplicates,]

    ## Records the total number of speficiations:
    s.curve.mod$n.specifications.original <-
      s.curve.mod$n.specifications

    ## Updates the actual number of remaining specifications:
    s.curve.mod$n.specifications <-
      nrow(s.curve.mod$results)

    ## Print message: report on identified duplicates:
    message("Duplicate models dropped, ",
            s.curve.mod$n.specifications.original,
            " now ",
            s.curve.mod$n.specifications,
            " (less ",
            s.curve.mod$n.specifications.original -
            s.curve.mod$n.specifications,
            ")")
  }


  ## (pemutation test-not used in Bryan, Yeager, O'Brien) ----
  perm.dat <- NULL

  if(is.numeric(permutations)){

    # Set the RNG seed:
    if(!is.null(s.curve.mod$perm.seed)){
      # use specified seed if available:
      set.seed(s.curve.mod$perm.seed)
    } else {
      # Otherwise generate one and record it:
      s.curve.mod$perm.seed <- sample(1:1000, 1)
      set.seed(s.curve.mod$perm.seed)
    }

    ## Save a nested list of individual data frames
    ## like the one we just made above:
    perm.dat <- dat
      replicate(
        ## number of permutations
        permutations,
        {
          ## shuffle all treatment vars
          ## Same ordering for all variables in a particular
          ## iteration
          new.order <- sample(nrow(perm.dat))
          perm.dat[treatment] <-
            lapply(
              treatment,
              function(tvar){
                perm.dat[[tvar]][new.order]
              })

          ## run model for that permutation:
          curveRunner(s.curve.mod,
                      inc.full.models = FALSE,
                      inc.tidy.models = FALSE)
        },
        simplify=FALSE)
  } #

  ## .(4b & 4c) Model results and output ----

  ## Calls s-curve update to provide final output:
  s.curve.update(
    s.curve.mod = s.curve.mod,
    perm.dat = perm.dat
  )

}


