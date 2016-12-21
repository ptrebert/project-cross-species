#!/bin/bash

INSTALLPATH="/TL/epigenetics2/work/pebert/conda/envs/pplcs"

echo "Conda environments can currently not be activated/deactivated in parallel (e.g., in a cluster environment)"
echo "See this github issue: https://github.com/conda/conda/issues/2502"
echo "Hence, the pipeline environment cannot be auto-activated from within Pied Piper"

exit 1

cd ${INSTALLPATH} &&
mkdir -p ./etc/conda/activate.d &&
mkdir -p ./etc/conda/deactivate.d &&
touch ./etc/conda/activate.d/env_vars.sh &&
touch ./etc/conda/deactivate.d/env_vars.sh

CRPBINPATH="/home/pebert/work/code/mpggit/creepiest"
CRPPYTHONPATH="/home/pebert/work/code/mpggit/creepiest"

echo "#!/bin/sh" > ./etc/conda/activate.d/env_vars.sh

echo "" >> ./etc/conda/activate.d/env_vars.sh

echo "export PATH=${CRPBINPATH}:\$PATH" >> ./etc/conda/activate.d/env_vars.sh

echo "export PYTHONPATH=${CRPPYTHONPATH}:\$PYTHONPATH" >> ./etc/conda/activate.d/env_vars.sh

echo "Environment adapted"
