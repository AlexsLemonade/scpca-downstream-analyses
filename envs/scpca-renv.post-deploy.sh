#!/bin/sh

# install packages required for renv
Rscript --vanilla -e \
  "
   .libPaths(.Library)
   install.packages(c('jsonlite', 'purrr'), repos = 'https://cloud.r-project.org')
   library(jsonlite)
   library(purrr)
   source('.Rprofile')
   renv::restore()
  "

