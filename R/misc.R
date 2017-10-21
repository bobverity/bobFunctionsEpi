#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.

#' @useDynLib bobFunctionsEpi
#' @importFrom Rcpp evalCpp
#' @importFrom odin odin
NULL

# -----------------------------------
#' functionTimer
#'
#' Calls an Rcpp function that includes a timer around an editable block. Useful for comparing speed of different function implementations.
#'
#' @param reps number of times to repeat block
#'
#' @export

functionTimer <- function(reps=1) {
    functionTimer_cpp(reps)
}

