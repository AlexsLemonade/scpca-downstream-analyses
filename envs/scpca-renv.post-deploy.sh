#!/bin/sh

# set conda subdir if on Apple Silicon
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  conda config --env --set subdir osx-64
fi

SCPCATOOLS_VERS="v0.2.1"

# remove any settings for R_LIBS_USER that might have crept in (renv/RStudio?)
R_LIBS_USER="" 

# install github packages
Rscript --vanilla -e \
  "
  remotes::install_github('AlexsLemonade/scpcaTools', ref='${SCPCATOOLS_VERS}', upgrade='never')
  require(scpcaTools) # check installation
  "

