## Custom functions for setting up renv

# #### USAGE
# This script is intended to be sourced in the script as follows:
#
# source(file.path("utils", "setup-functions.R"))

setup_renv <- function(project_filepath = here::here(),
                       restore_packages = TRUE) {
  # Purpose: Load the project file and install any packages/dependecies needed
  #          to run the workflow
  #
  # Args:
  #   project_filepath: filepath to the project file to be loaded in, the
  #                     default is the `here::here()` function to identify the
  #                     path of the root directory, including whether or not
  #                     there is a .Rproj file present
  #   restore_packages: logical argument to determine whether or not packages
  #                     should be restored to what's in the renv.lock file;
  #                     default is TRUE

  # Source the .Rprofile file in case it was not read at launch
  source(file.path(project_filepath, ".Rprofile"))

  if (restore_packages == TRUE) {
    # install any necessary packages and dependencies from the renv.lock file
    renv::restore()
  }
}


check_r_bioc_versions <- function() {

  # Purpose: Check if R and Bioconductor versions are at appropriate minimums,
  #  and throw an error if they are not.
  # No input arguments are required and nothing is returned.

  # Check that R version is at least 4.2
  if (! (R.version$major == 4 && R.version$minor >= 2)){
    stop("R version must be at least 4.2.")
  }

  # Check that Bioconductor version is at least 3.15
  if (packageVersion("BiocVersion") < 3.15){
    stop("Bioconductor version is less than 3.15.")
  }

}

