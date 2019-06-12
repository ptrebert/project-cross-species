#!/usr/bin/env bash

REPOSITORY_PREFIX="/home/pebert/work/code/github/project-cross-species"

CONDA_PREFIX="/TL/epigenetics2/work/pebert/conda/envs/cse"

cd ${CONDA_PREFIX}
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
#touch ./etc/conda/activate.d/env_vars.sh
#touch ./etc/conda/deactivate.d/env_vars.sh