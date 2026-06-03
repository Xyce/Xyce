#!/bin/bash -xe

cleanInstallDir() {
  local installDir="$1"
  if [ -d "${installDir}" ]; then
    cd "$(dirname ${installDir})"
    rm -rf "$(basename ${installDir})"
  else
    mkdir -p "${installDir}"
  fi
}

## The following are only valid when script is submitted via
## sbatch. For salloc these must be specified on the command line.
nodes=$SLURM_JOB_NUM_NODES          # Number of nodes

physCoreCount=$(lscpu -p | grep -v '^#' | sort -u -t, -k 2,4 | wc -l)
# divide by two assuming xyce jobs are using "mpirun -np 2 ..."
coreCount=$(( physCoreCount/2 ))

rawCmakeArgsList="$1"
xyceInstallDir="$2"
CI_PROJECT_DIR="$3"
ExcludeTestLabel="$4"

PATH="/projects/xyce/flexbison/bin:${PATH}"
export PATH

# activate python virtual for xyce build and test
source "/projects/xyce_user/pythonVirtEnv/Xyce/bin/activate"

# note use of the use of the pipeline's xyce-ctest.cmake file, NOT the
# build repos copy
ctest --timeout 1200 -DVERBOSITY=5 \
  -DNUM_PROCS=${coreCount} \
  -DUSE_GITLAB_CI_TESTING=ON \
  -DEXCLUDE_TESTS_WITH_LABEL="${ExcludeTestLabel}" \
  -DCMAKE_ARGS_LIST="${rawCmakeArgsList}" \
  -S "${CI_PROJECT_DIR}/cmake/ctest/xyce-ctest.cmake"
