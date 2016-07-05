#' Installs brainarray CDF and probe information packages
#'
#' This package downloads, installs, and/or loads the brainarray CDF and probe
#' enviroments from the Brainarray website.
#'
#' @param array A character of length 1 giving the giving the name of the
#'   array.
#' @param version A character of length 1 giving the version to download.
#' @param type A charcter giving the annotation type.
#' @param force.reinstall logical. Should the package be installed again?
#' @param force.download logical. Should the package be downloaded again?
#' @param use.temp.dir logical. Should a temporary dir be used?
#' @param path A charcter of length 1. The path describing where the downloaded
#'   files are saved.
#' @return Invisibly returns the location of the saved files.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau (at) gmail.com>
#' @seealso \code{\link{tempdir}}, \code{\link{install.packages}}
#' @references
#'   See \url{http://brainarray.mbni.med.umich.edu}.
#' @examples
#' tmp.path <- install.brainarray("hgu133a", version = "18.0.0", type = "ensg")
#' print(tmp.path)
#' @export

# Use a package (e.g. digest) to read the brainarray package and avoid most
# hardcoded functionality of this function!

install.brainarray <- function(array,
                               version = "19.0.0",
                               type = "entrezg",
                               force.reinstall = FALSE,
                               force.download = FALSE,
                               use.temp.dir = TRUE,
                               path) {
  array <- "hgu133plus2";
  #force.reinstall <- TRUE; force.download <- TRUE; use.temp.dir <- TRUE

  # Construct the file names used
  files <- paste0(array, "hs", type, c("cdf_", "probe_", ".db_"), version, ".tar.gz")

  # Load or install if not installed
  package.name <- paste0(array, "hs", type, c("cdf", "probe", ".db"))

  # Prepare to download
  # Base URL to the brainarray site
  base.url <-
    paste0("http://mbni.org/customcdf/",
           version, "/" , type, ".download")

    # Path to download files to
  if (missing(path)) {
    path <- ifelse(use.temp.dir,
                   tempdir(),
                   file.path(array, "Brainarray"))
  }

  # Loop through the packages
  for (i in seq_along(package.name)) {

    pkg <- package.name[i]

    # Is pkg installed... ?
    if (!(force.download | force.reinstall) &
          requireNamespace(pkg)) {

      # ... and is it the correct version?
      if (packageVersion(pkg) == version) {
        next  # if yes, then do the next pkg
      }

    } else {

      # Create the dir
      dir.create(path = path, recursive = TRUE, showWarnings = FALSE)

      # Download the files if they are not already
      if (!file.exists(file.path(path, files[i])) | force.download) {
        download.file(url = file.path(base.url, files[i]),
                      destfile = file.path(path, files[i]))
      }

      # Install package
      install.packages(pkgs = file.path(path, files[i]),
                       repos = NULL, type = "source")

      # Require package
      requireNamespace(pkg)
    }

  }

  return(invisible(file.path(path, files)))
}




