#!/bin/sh

# install packages required for renv
Rscript --vanilla -e \
  "
   library(glue)
   library(jsonlite)
   library(purrr)
   source('.Rprofile')
   renv::restore()
  "

