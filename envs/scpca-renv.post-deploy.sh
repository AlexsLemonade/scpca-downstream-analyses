#!/bin/sh

# Record CONDA_SUBDIR if it is set to specify Intel on an Apple Silicon mac
# Used for compatibility with Bioconductor packages that may not work with Apple Silicon
if [[ $CONDA_SUBDIR == 'osx-64' ]]; then
  conda env config vars set CONDA_SUBDIR=osx-64
fi

# restore packages
Rscript -e 'renv::restore()'
