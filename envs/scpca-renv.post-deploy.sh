#!/bin/sh

# install packages required for renv
Rscript --vanilla -e \
  "
   library(tidyverse)
   source('.Rprofile')
   renv::restore()
  "

