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

# Number MPI processes to run on each node (a.k.a. PPN)
#
# the multiplier needs to be the smallest number of cores per node on
# the list of cts machines being used for testing. current list of
# cores per node on the relevant machines is:
#    amber|flight|hops: 112
#    manzano|jetcool: 48
#    attaway|ghost: 36
case ${SLURM_CLUSTER_NAME} in
  ghost|attaway|eclipse)
    cmult=36
    ;;
  amber|flight|hops)
    cmult=112
    ;;
  manzano|jetcool)
    cmult=48
    ;;
  *)
    cmult=12
    ;;
esac

# divide by two assuming xyce jobs are using "mpirun -np 2 ..."
cores=$(( ${nodes}*${cmult}/2 ))
rawCmakeArgsList="$1"
xyceInstallDir="$2"
CI_PROJECT_DIR="$3"

# note use of the use of the pipeline's xyce-ctest.cmake file, NOT the
# build repos copy
ctest --timeout 1200 -DVERBOSITY=5 \
      -DNUM_PROCS=28 \
      -DUSE_GITLAB_CI_TESTING=ON \
      -DCMAKE_ARGS_LIST="${rawCmakeArgsList}" \
      -S ${CI_PROJECT_DIR}/cmake/ctest/xyce-ctest.cmake
