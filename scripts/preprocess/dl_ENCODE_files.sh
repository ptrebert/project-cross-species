#!/bin/bash

TARGET_DIR=/TL/deep/fhgfs/projects/pebert/thesis/projects/cross_species/rawdata/downloads

LISTING=/TL/deep/fhgfs/projects/pebert/thesis/projects/cross_species/rawdata

cd ${TARGET_DIR}

CMD="xargs -n 1 -a ${LISTING} curl -L -O -k -s -S"

${CMD}

exit $?

