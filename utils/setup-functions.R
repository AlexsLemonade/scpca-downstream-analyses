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
  #                     default here being the `here::here()` function
  #   restore_packages: binary argument to determine whether or not packages
  #                     should be restored to what's in the renv.lock file;
  #                     default here is TRUE
  
  # `here::here()` looks at a number of criteria to identify the root 
  # directory, including whether or not there is a .Rproj file present,
  # so we can pass this to `renv::load()` to load the project file
  renv::load(here::here())
  
  if (restore_packages == TRUE) {
    # install any necessary packages and dependecies from the renv.lock file
    renv::restore()
  }
}
