#'@export
#'
#'@title Goodness of fit testing of CMR data for a CJS model
#'
#'@description Run goodness of fit tests from R2ucare on CJS processed data
#'
#'@param data.processed (required) a list with several variables
#'  (e.g., number of capture occasions, time intervals) returned by
#'  \code{\link[RMark]{process.data}()} based on live-only capture histories and
#'  model = "CJS".
#'
#'@param verbose (optional, default = \code{FALSE}) if TRUE, prints test details
#'
#'@param keep.intermediate (optional, default FALSE) this function creates and
#'  intermediate file (\code{CH.for.GOF.inp}) which is removed by default before
#'  returning. If this option is \code{TRUE}, then the intermediate file is not
#'  removed.
#'
#'@details This function produces a dataframe of summary stats from each R2ucare
#'  test with the option to print details of each test to the console.
#'  Note: future versions will allow the user to specify group covariates, but
#'  currently this feature is unavailable.
#'
#'
#'@return Returns a dataframe with goodness of fit testing results from overall
#'  CJS, and tests 2CT, 3SR, 3SM, 2CL
#'@author Sarah Gutowsky

CJS.GOF.testing <- function(data.processed, verbose = FALSE)
{
  # export capture history from the data.processed dataframe as .inp file and
  # import again to use with R2ucare
  # TODO: should really use a tempfile here so as not to clutter up the
  # working directory.
  RMark::export.chdata(data.processed, filename = "CH.for.GOF" , replace =
                         TRUE)
  data.inp = R2ucare::read_inp("CH.for.GOF.inp")
  data.hist = data.inp$encounter_histories
  data.freq = data.inp$sample_size
  GOF.tests.CJS <- data.frame(matrix(ncol = 5, nrow = 5))
  colnames(GOF.tests.CJS) <-
    c("test", "statistic", "df", "signed.test", "p.value")
  overall.CJS <- R2ucare::overall_CJS(data.hist, data.freq)
  GOF.tests.CJS[1, ] <-
    c(
      "overall CJS",
      overall.CJS$chi2,
      overall.CJS$degree_of_freedom,
      "not applicable",
      overall.CJS$p_value
    )
  test2ct <- R2ucare::test2ct(data.hist, data.freq, verbose = verbose)
  GOF.tests.CJS[2, ] <-
    c("Test2CT",
      test2ct$test2ct["stat"],
      test2ct$test2ct["df"],
      test2ct$test2ct["sign_test"],
      test2ct$test2ct["p_val"])
  test3sr <- R2ucare::test3sr(data.hist, data.freq, verbose = verbose)
  GOF.tests.CJS[3, ] <-
    c("Test3SR",
      test3sr$test3sr["stat"],
      test3sr$test3sr["df"],
      test3sr$test3sr["sign_test"],
      test3sr$test3sr["p_val"])
  test3sm <- R2ucare::test3sm(data.hist, data.freq, verbose = verbose)
  GOF.tests.CJS[4, ] <-
    c("Test3SM",
      test3sm$test3sm["stat"],
      test3sm$test3sm["df"],
      "not applicable",
      test3sm$test3sm["p_val"])
  test2cl <- R2ucare::test2cl(data.hist, data.freq, verbose = verbose)
  GOF.tests.CJS[5, ] <-
    c("Test2CL",
      test2cl$test2cl["stat"],
      test2cl$test2cl["df"],
      "not applicable",
      test2cl$test2cl["p_val"])

  if (isTRUE(verbose)) {
    message("Test2CT details:")
    print(test2ct$details)
    message("\nTest3SR details:")
    print(test3sr$details)
    message("\nTest3SM details:")
    print(test3sm$details)
    message("\nTest2CL details:")
    print(test2cl$details)
  }

  # Clean up
  if(unlink("CH.for.GOF.inp") == 1)
    warning("CJS.GOF.testing: Failed to remove temporary file 'CH.for.GOF.inp'")

  return(GOF.tests.CJS)

}
