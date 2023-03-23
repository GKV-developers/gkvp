#!/bin/sh

###  NOTE  ###
###  Flow supercomputer Type I sub-system, PRIMEHPC FX1000 (Nagoya Univ, 2020)
###
###  - Computation nodes(total 2304 nodes)
###      CPU: A64FX (2.0GHz, 12coresx4CMG=48cores, 512bit SIMD) x1 per node
###      Peak performance: DP 3.379 TFLOPS per node (Boost: 3.3792 TFLOPS)
###      Cache L1: 64 KiB, 4 way
###      Cache L1 Bandwidth: 230+ GB/s(load), 115+ GB/s (store)
###      Cache L2: 8 MiB, 16 way per CMG(NUMA), 4CMG per node
###      Cache L2 Bandwidth: 3.6+ TB/s per node
###                          115+ GB/s(load), 57+ GB/s(store) per core
###      Memory: HBM2 32 GiB
###      Memory Bandwidth: 1024 GB/s per node
###
###      Therefore, a recommended GKV parallelization may be 
###          (MPI Processes)x(12 OpenMP Threads)
###          =(12 cores per CMG)x(4 CMG)x(Node numbers)
###      1 MPI process should be assigined to 1 CMG.
###
###  - Interconnect
###      Tofu Interconnect D (28 Gbps x 2 lane x 10 port)
###      [Performance] 8B Put latency: 0.49-0.54 usec
###                    1MiB Put throughput: 6.35 GB/s
###
###  - Job class (May 2020)
###      fx-debug  :  1 - 36  nodes,   1 hour,  50 run/300 submit
###      fx-small  :  1 - 24  nodes, 168 hour, 100 run/300 submit
###      fx-middle : 12 - 96  nodes,  72 hour,  50 run/300 submit
###      fx-large  : 96 - 192 nodes,  72 hour,  25 run/300 submit
###      fx-xlarge : 96 - 768 nodes,  24 hour,   5 run/300 submit
###
###  - Commands
###      (Submit a batch job : "pjsub sub.q") Use shoot script for GKV.
###      Check job status    : "pjstat" or "pjstat -E" for step jobs
###      Delete job          : "pjdel JOBID"
###      Show budget info    : "charge"
###      Show disk usage     : "lfs quota -u (YOUR ACCOUNT ID) /home"
###                          : "lfs quota -u (YOUR ACCOUNT ID) /data"
##############

#PJM --rsc-list "rscgrp=fx-debug"
#PJM --rsc-list "node=8"       
#### --rsc-list "node=5x8x8"       
#PJM --rsc-list "elapse=00:10:00"
#PJM --mpi "proc=32"           
#### --mpi "rank-map-bynode"
#### --mpi "rank-map-hostfile=rankmapfile.dat"
#PJM -j                          
#PJM -s                           

NUM_NODES=${PJM_NODE}             # Nodes
NUM_CORES=12                      # Cores per node
NUM_PROCS=$(( ${NUM_NODES} * 4 )) # MPI processes
export OMP_NUM_THREADS=12         # OpenMP threads per MPI


echo "                  Nodes: ${NUM_NODES}"
echo "         Cores per node: ${NUM_CORES}"
echo "          MPI Processes: ${NUM_PROCS}"
echo " OpenMP threads per MPI: ${OMP_NUM_THREADS}"



### Working directory 
DIR=%%DIR%%
LDM=gkvp.exe
NL=gkvp_namelist.%%%

export XOS_MMM_L_PAGING_POLICY=demand:demand:demand # For Largepage

export PLE_MPI_STD_EMPTYFILE="off" # Suppress stdout of filesize-0

module load fftw-tune phdf5 netcdf-c netcdf-fortran
###module unload tcs
###module load fftw/3.3.8
###export PATH=/opt/FJSVxtclanga/tcsds-1.2.25/bin:$PATH
###export LD_LIBRARY_PATH=/opt/FJSVxtclanga/tcsds-1.2.25/lib64:$LD_LIBRARY_PATH
###export OPAL_PREFIX=/opt/FJSVxtclanga/tcsds-1.2.25


#### Run
date
cd ${DIR}
export fu05=${DIR}/${NL}
mpiexec -n ${NUM_PROCS} ${DIR}/${LDM}
   # -n        "Total number of MPI processes"
date


##### Run with Fujitsu profiler fipp (re-compile with -Nfjprof option)
#date
#cd ${DIR}
#export fu05=${DIR}/${NL}
#fipp -C -d ${DIR}/fjprof_dir/pa0 -Icpupa -Impi -Sregion  mpiexec -n ${NUM_PROCS} ${DIR}/${LDM}
#date
#echo "#!/bin/sh" > ${DIR}/fjprof_dir/fugaku_fipppx.sh
#echo "set -Ceu" >> ${DIR}/fjprof_dir/fugaku_fipppx.sh
#echo "set -x" >> ${DIR}/fjprof_dir/fugaku_fipppx.sh
#echo "fipppx -A -d pa0 -Icpupa     -p0,limit=4 -o prof_cpupa.txt" >> ${DIR}/fjprof_dir/fugaku_fipppx.sh
#echo "fipppx -A -d pa0 -Ibalance   -p0,limit=4 -o prof_balance.txt" >> ${DIR}/fjprof_dir/fugaku_fipppx.sh
#echo "#fipppx -A -d pa0 -Icall      -p0,limit=4 -o prof_call.txt" >> ${DIR}/fjprof_dir/fugaku_fipppx.sh
#echo "fipppx -A -d pa0 -Isrc:./src -p0,limit=4 -o prof_src.txt" >> ${DIR}/fjprof_dir/fugaku_fipppx.sh


##### Run with Fujitsu profiler fapp (re-compile with -Nfjprof option)
#date
#cd ${DIR}
#export fu05=${DIR}/${NL}
#Npa=1  # Elementary report
##Npa=5  # Simple report
##Npa=11 # Standard report
##Npa=17 # Detailed report
#for i in `seq 1 ${Npa}`; do
#  echo "pa"${i} `date`
#  fapp -C -d ${DIR}/fjprof_dir/pa${i} -Hevent=pa${i} -Sregion  mpiexec -n ${NUM_PROCS} ${DIR}/${LDM}
#done
#date
#
#echo "#!/bin/sh" > ${DIR}/fjprof_dir/fugaku_fapppx.sh
#for i in `seq 1 ${Npa}`; do
#  echo "fapppx -A -d ./pa${i} -Icpupa,mpi -tcsv -o pa${i}.csv" >> ${DIR}/fjprof_dir/fugaku_fapppx.sh
#done
#echo "cp /opt/FJSVxtclanga/tcsds-1.2.25/misc/cpupa/cpu_pa_report.xlsm ./" >> ${DIR}/fjprof_dir/fugaku_fapppx.sh
#
#
