#!/bin/sh

# set conda subdir if on Apple Silicon
if [[ "$(uname)" == 'Darwin' && "$(uname -m)" == 'arm64' ]]; then
  conda config --env --set subdir osx-64
fi


SCPCATOOLS_VERS='main'

# install packages not on conda-forge
Rscript --vanilla -e \
  "
  # install Harmony from CRAN
  remotes::install_version(
    'harmony', 
    version = '0.1.1', 
    repos = 'https://cloud.r-project.org', 
    upgrade = FALSE
  )

  # install ScPCA tools from Github
  remotes::install_github('AlexsLemonade/scpcaTools', ref='${SCPCATOOLS_VERS}', upgrade='never')
  require(scpcaTools) # check installation
  "

