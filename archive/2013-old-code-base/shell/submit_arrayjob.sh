#! /bin/bash

INSTRFILE=$1
LOGNAME=$2
JOBNAME=$3

commands2arrayjob.sh ${INSTRFILE} -cwd -j y -o /scratch/TL/pool0/pebert/splitfiles/${LOGNAME} -l h_rt=168:: -b y -N ${JOBNAME}

# options explained:
# -cwd: execute in cwd
# -j y: merge out/err streams into one
# -o: send stdout to...
# -l: set limits, in this case 7 day queue (168hrs)
# -b y: activate binary mode (expect binary, not shell script to be executed)
# -N: set job name
