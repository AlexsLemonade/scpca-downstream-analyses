#!/bin/sh

# set conda subdir if on Apple Silicon
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  conda config --env --set subdir osx-64
fi

SCPCATOOLS_VERS="v0.2.1"

# install github packages
Rscript --no-init-file -e \
  "
  remotes::install_github('AlexsLemonade/scpcaTools', ref='${SCPCATOOLS_VERS}', upgrade='never')
  "

