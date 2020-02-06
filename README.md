# SpecCurve
R package for constructing specification curves


This R package is intended to benefit research, data analysis, and open/reproducible science through helping analysts construct, run, and summarize "specificaion curves" (Simonsohn, Simmons, and Nelson, 2015) to examine the size, robustness, and sensitivity of analyzed treatment effects.   The tools in this package are designed to streamline the work to meet the three goals outlined in the original paper: 

i. identifying (and modeling) the set of theoretically justified, statistically valid, and non-redundant analytic specifications
ii. displaying alternative results graphically, allowing the identification of decisions producing different results
iii. conducting statistical tests to determine whether as a whole results are inconsistent with the null hypothesis

With the tools in this package, one can specify sets of choices from which to build models (currently, linear and generalized linear models using *lm()* and *glm()* respectively).  Potential choices currently include:

1. varying outcome variables (e.g., different transformations of outcome)
2. differing treatment variables
3. differing covariates (selecting 0 or 1 specifications of covariate(s) to include from each element of a list)
4. differing data subsets
5. weightings applied to individual cases in analysis

Given that, and an appropriate dataset, the main function s.curve() will construct a set of all the unique model specifications available from the choices above, running those, and summarizing the results.  Further tools allow for visualizations of the results to assisting goal ii, as well as testing the results of the specification set against a larger set of specification curves employing permutations of the treatment variable within the data to simulate an exact test.

This package is a work in progress, so I welcome any pull requests from interested parties or inquiries about collaboration.
