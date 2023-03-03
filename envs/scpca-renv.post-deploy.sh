#!/bin/sh

# set conda subdir if on Apple Silicon
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  conda config --env --set subdir osx-64
fi

# install github packages
Rscript -e \
  "
  remotes::install_github('AlexsLemonade/scpcaTools', ref='v0.1.8', upgrade='never')
  "

