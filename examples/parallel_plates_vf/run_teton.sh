#!/bin/bash
set -e

# Usage (run from this directory, examples/parallel_plates_vf/):
#   PROJ_ID=<your wckey> ./run_teton.sh parallel_plates 1 00:30
#
# Same generator pattern as nrsqsub_frontier / nrsqsub_teton-2.sh: writes an
# sbatch script to s.bin and submits it. parallel_plates_vf is a small
# radiation-only verification case (4400 hex elements, polynomialOrder=5,
# numSteps=5), so 1 node and a short walltime are enough -- see
# procedure.md in this folder.
#
# TASKS_PER_NODE default (32) is deliberately far below a full Turin node
# (192, nrsqsub_teton-2.sh's default): concentric_spheres_vf (1470
# elements) segfaulted at 192 ranks -- nekRS auto-sized the legacy
# Nek5000 per-rank element buffer too small for the resulting imbalanced
# partition. 32 ranks -> ~138 elements/rank here, comfortably safe.
# Override with TASKS_PER_NODE=<n> if you want more parallelism once this
# runs cleanly.

: ${PROJ_ID:="art"} # FIXME somehow not used??
: ${QUEUE:="short"} # short / general
: ${TASKS_PER_NODE:=32}

: ${NEKRS_HOME:=$HOME/.local/nekrs}
export NEKRS_HOME
export PATH="$NEKRS_HOME/bin:$PATH"

source $NEKRS_HOME/bin/nrsqsub_utils
setup $# 1

#cpu_per_node=384
cpu_per_node=$TASKS_PER_NODE
let nn=$nodes*$cpu_per_node
let ntasks=nn
backend=serial

chk_case $ntasks

# sbatch
SFILE=s.bin
echo "#!/bin/bash" > $SFILE
#echo "#SBATCH -A $PROJ_ID" >>$SFILE
echo "#SBATCH -J $jobname" >>$SFILE
echo "#SBATCH -o %x-%j.out" >>$SFILE
echo "#SBATCH -t ${time}:00" >>$SFILE
echo "#SBATCH -N $qnodes" >>$SFILE
echo "#SBATCH -p $QUEUE" >>$SFILE
echo "#SBATCH --wckey=${PROJ_ID}" >> $SFILE
echo "#SBATCH --exclusive" >>$SFILE
echo "#SBATCH --mem=0" >>$SFILE
echo "#SBATCH --ntasks-per-node=$cpu_per_node" >>$SFILE
echo "#SBATCH --cpus-per-task=1" >>$SFILE
echo "#SBATCH --hint=nomultithread" >>$SFILE


echo "module purge" >> $SFILE
echo "module load PrgEnv-gnu" >> $SFILE
echo "module load craype-x86-turin" >> $SFILE
echo "module load craype-network-ofi" >> $SFILE
echo "module load cray-mpich" >> $SFILE
echo "module load cmake" >> $SFILE
echo "module list" >> $SFILE

echo "squeue -u \$USER" >>$SFILE

echo "export MPICH_GPU_SUPPORT_ENABLED=0" >>$SFILE

echo "ulimit -s unlimited " >>$SFILE
echo "export NEKRS_HOME=$NEKRS_HOME" >>$SFILE
echo "export NEKRS_GPU_MPI=0" >>$SFILE

echo "export MPICH_MPIIO_STATS=1" >>$SFILE
echo "export MPICH_OFI_NIC_POLICY=NUMA" >>$SFILE

echo "export FI_CXI_RX_MATCH_MODE=hybrid" >> $SFILE

echo "export PMI_MMAP_SYNC_WAIT_TIME=600" >> $SFILE

echo "" >> $SFILE
echo "date" >>$SFILE
echo "" >> $SFILE

echo "ldd $bin" >> $SFILE
#echo "lscpu" >> $SFILE
#echo "lstopo-no-graphics" >> $SFILE
#echo "OMP_NUM_THREADS=1 srun --mpi=cray_shasta ./hello_jobstep" >>$SFILE

if [ $RUN_ONLY -eq 0 ]; then
  echo -e "\n# precompilation" >>$SFILE
  CMD_build="srun --nodes 1 --ntasks=$cpu_per_node --mpi=cray_shasta $bin --backend $backend --device-id 0 $extra_args --setup \$case_tmp --build-only \$ntasks_tmp"
  add_build_CMD "$SFILE" "$CMD_build" "$ntasks"
fi


if [ $BUILD_ONLY -eq 0 ]; then
  link_neknek_logfile "$SFILE"
  echo -e "\n# actual run" >>$SFILE
  echo "srun --mpi=cray_shasta $bin --backend $backend --device-id 0 $extra_args --setup $case" >>$SFILE
fi
sbatch $SFILE

# clean-up
