## NOTE there is now a CMRfunctions.R file for functions specific
#to the CMR modeling process


# Remove all mark .inp, .res, .out, .tmp, and .vcv files from a folder.
#
# NOTE this is a little more liberal than RMark::cleanup() since it will match
# ANY file starting with "prefix" and ending with one of ".inp", ".out", etc.
#
# Arguments
#
# folder - folder to delete files from. If specified, folder is relative to
#    project root as determined using here(). Note this function doesn't require
#    the "here" package to be loaded, but rather just installed.
#
# prefix - character string that all temporary files start with, default =
#    "mark".
#
# recursive - if TRUE, searches recursively down through folder hierarchy for
#    files to delete.
#
# test.only - if TRUE, returns the list of files to be deleted without removing
#    them.
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
