#'
#'@export
#'
#'@title Force cleanup of temporary files
#'
#'@description Remove all mark .inp, .res, .out, .tmp, and .vcv files from a
#'    folder.
#'
#' @param folder (optional, character vector or length 1) folder to delete files from.
#'    If absent, defaults to the project root as determined by \link[here]{here}.
#'    If specified, folder is relative to project root. Note this function
#'    doesn't require the "here" package to be loaded, but rather just installed.
#'
#' @param prefix (optional, character vector of length 1) character string
#'     that all temporary files start with, default = "mark".
#'
#' @param recursive (optional, boolean, default = FALSE) If TRUE, searches
#'  recursively down through folder hierarchy for files to delete.
#'
#' @param test.only (optional, boolean, default = FALSE) If TRUE, returns the
#'  character vector of file names to be deleted without removing them.
#'
#'@details
#'   NOTE this is a little more liberal than RMark::cleanup() since it will match
#'  ANY file starting with "prefix" that ends with one of ".inp", ".out", etc.
#'
#'@author Dave Fifield
force.cleanup <- function(folder = "", prefix = "mark", recursive = FALSE,
                          test.only = FALSE) {
  folder <- here::here(folder)
  pat <- sprintf("^%s.*\\.inp$|^%s.*\\.res$|^%s.*\\.out$|^%s.*\\.vcv$|^%s.*\\.tmp$",
                 prefix, prefix, prefix, prefix, prefix)
  files <- list.files(path = folder, pattern = pat, full.names = TRUE, recursive = recursive)

  if (isFALSE(test.only))
    invisible(file.remove(files))
  else
    files
}
