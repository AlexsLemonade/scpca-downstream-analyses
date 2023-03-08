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
  
  # If we are in a snakemake/conda environment, do nothing
  is_conda <- grepl(".snakemake/conda", .libPaths(), fixed = TRUE)
  if(any(is_conda)){
    message("Using snakemake/conda: skipping renv")
    return()
  }
  
  # Check if there is an .Rprofile & source it
  rprofile <- file.path(project_filepath, ".Rprofile")
  if(file.exists(rprofile)){
    source(rprofile)
  }
  
  # check if renv is set up
  renv_unset <- is.null(options("renv.project.path"))
  if(renv_unset){
    warning("renv is not set up; package versions may differ from expectations")
  }else if(restore_packages) {
    # install any necessary packages and dependencies from the renv.lock file
    renv::restore()
  }
}


check_r_version <- function() {

  # Purpose: Check if R version is at appropriate minimum
  # No input arguments are required and nothing is returned.

  # Check that R version is at least 4.2
  if (! (R.version$major == 4 && R.version$minor >= 2)){
    stop("R version must be at least 4.2.")
  }

}

