#!/bin/sh

# set conda subdir if on Apple Silicon
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  conda config --env --set subdir osx-64
fi


SCPCATOOLS_VERS='main'

# install github packages
Rscript --vanilla -e \
  "
  remotes::install_github('AlexsLemonade/scpcaTools', ref='${SCPCATOOLS_VERS}', upgrade='never')
  require(scpcaTools) # check installation
  "

