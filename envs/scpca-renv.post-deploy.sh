#!/bin/sh

# install packages required for renv
Rscript --vanilla -e \
  "
   library(tidyverse)
   library(jsonlite)
   source('.Rprofile')
   renv::restore()
  " 
