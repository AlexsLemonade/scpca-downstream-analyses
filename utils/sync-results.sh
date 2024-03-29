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
mkdir -p integration-example-results/

# files only found in output from running core worfklow
core_link_locs=(
  sample01/library01_processed.rds
  sample01/library01_core_analysis_report.html
  sample02/library02_processed.rds
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

# genes of interest module output files
goi_link_locs=(
  sample01/library01_goi_report.html
  sample01/library01_goi_stats
  sample02/library02_goi_report.html
  sample02/library02_goi_stats
)

# data integration module output files
integration_link_locs=(
  group01_integration_report.html
  group01_integrated_sce.rds
  group01_merged_sce.rds
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

for loc in ${goi_link_locs[@]}
do
  # only make the links if replacing an old link or the file doesn't exist
  if [[ -L ${loc} || ! -e ${loc} ]]
  then
    ln -nsf ${repo_base}/example-results/${loc} goi-example-results/${loc}
  else
    echo "${loc} already exists and is not a link, delete or move it to create a link."
  fi
done

for loc in ${integration_link_locs[@]}
do
  # only make the links if replacing an old link or the file doesn't exist
  if [[ -L ${loc} || ! -e ${loc} ]]
  then
    ln -nsf ${repo_base}/example-results/${loc} integration-example-results/${loc}
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

# zip and sync genes of interest module example results
zip -r goi_example_results.zip $repo_base/goi-example-results
aws s3 cp goi_example_results.zip $s3_base/goi_example_results.zip --acl public-read
rm ./goi_example_results.zip

# zip and sync integration module example results
zip -r integration_example_results.zip $repo_base/integration-example-results
aws s3 cp integration_example_results.zip $s3_base/integration_example_results.zip --acl public-read
rm ./integration_example_results.zip
