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
mkdir -p clustering-example-results/sample01
mkdir -p clustering-example-results/sample02

# files only found in output from running core worfklow
core_link_locs=(
  sample01/library01_processed_sce.rds
  sample01/library01_core_analysis_report.html
  sample02/library02_processed_sce.rds
  sample02/library02_core_analysis_report.html
)

# clustering module output files
clustering_link_locs=(
  sample01/library01_clustered_sce.rds
  sample01/library01_clustering_report.html
  sample01/library01_clustering_stats
  sample02/library02_clustered_sce.rds
  sample02/library02_clustering_report.html
  sample02/library02_clustering_stats
)

for loc in ${core_link_locs[@]}
do
  # only make the links if replacing an old link or the file doesn't exist
  if [[ -L ${loc} || ! -e ${loc} ]]
  then
    ln -nsf ${repo_base}/example-results/${loc} example-results/${loc}
  else
    echo "${loc} already exists and is not a link, delete or move it to create a link."
  fi
done

for loc in ${clustering_link_locs[@]}
do
  # only make the links if replacing an old link or the file doesn't exist
  if [[ -L ${loc} || ! -e ${loc} ]]
  then
    ln -nsf ${repo_base}/example-results/${loc} clustering-example-results/${loc}
  else
    echo "${loc} already exists and is not a link, delete or move it to create a link."
  fi
done

# zip and sync core analysis example results
zip -r core_example_results.zip $repo_base/example-results
aws s3 cp core_example_results.zip $s3_base/core_example_results.zip --acl public-read
rm ./core_example_results.zip

# zip and sync clustering module example results
zip -r clustering_example_results.zip $repo_base/clustering-example-results
aws s3 cp clustering_example_results.zip $s3_base/clustering_example_results.zip --acl public-read
rm ./clustering_example_results.zip
