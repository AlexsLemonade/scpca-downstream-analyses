#!/bin/sh

# Deactivate system set user library
unset R_LIBS_USER 
# and make the change permanent
conda env config vars set R_LIBS_USER=''

# set conda subdir if on Apple Silicon
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  conda config --env --set subdir osx-64
fi


SCPCATOOLS_VERS='v0.2.1'

# install github packages
Rscript --vanilla -e \
  "
  remotes::install_github('AlexsLemonade/scpcaTools', ref='${SCPCATOOLS_VERS}', upgrade='never')
  require(scpcaTools) # check installation
  "

