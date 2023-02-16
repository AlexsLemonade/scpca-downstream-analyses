#! /bin/bash

# This script is used to establish symlinks for scpca-downstream-analyses to use 
# data stored in a shared directory.

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
# move back up to the scpca-downstream-analyses root
cd ..

# location for the shared data
share_base=/shared/data
repo_base=${share_base}/scpca-downstream-analyses
s3_base=s3://scpca-references/example-data/scpca-downstream-analyses

# create directory for example results
mkdir -p example-results/sample01
mkdir -p example-results/sample02

link_locs=(
  example-results/sample01/library01_processed_sce.rds
  example-results/sample01/library01_core_analysis_report.html
  example-results/sample02/library02_processed_sce.rds
  example-results/sample02/library02_core_analysis_report.html
)

for loc in ${link_locs[@]}
do
  # only make the links if replacing an old link or the file doesn't exist
  if [[ -L ${loc} || ! -e ${loc} ]]
  then
    ln -nsf ${modules_base}/${loc} ${loc}
  else
    echo "${loc} already exists and is not a link, delete or move it to create a link."
  fi
done

# zip and sync example results
zip -r core_example_results.zip $repo_base/example-results
aws s3 cp core_example_results.zip $s3_base/core_example_results.zip --acl public-read
rm ./core_example_results.zip
