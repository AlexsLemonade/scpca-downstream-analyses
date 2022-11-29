#!/bin/bash
#
# Sets up conda environments for snakemake

set -euo pipefail

# Run from the script file location
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "Building conda environments"
snakemake --use-conda --conda-create-envs-only -c1
