#!/bin/bash

NUMBER_OF_CORES=40

# add additional paths to run snakemake on here
PATHS_TO_RUN=( "." "seqc" "simulations" "pseudobam" "shredding" "personalized_simulation" "speed_k" "bias_and_error")

ARGUMENTS="-p -j ${NUMBER_OF_CORES} --dryrun"

########################################################################

BASE_PATH=$(pwd)

for path in "${PATHS_TO_RUN[@]}"
do
  cd ${path} || exit
  echo "In path: ${path}"
  snakemake ${ARGUMENTS}
  cd ${BASE_PATH} || exit
done
