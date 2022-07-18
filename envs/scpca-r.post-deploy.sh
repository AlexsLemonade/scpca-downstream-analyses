#!/bin/bash
set -euo pipefail

Rscript -e 'renv::restore()'
