#!/bin/sh

#### Available HPC systems
# HELIOS  : Helios supercomputer, Bullx B510 (IFERC)
# PS      : Plasma simulator, Fujitsu FX100 (NIFS)
# NU      : NU FX100, Fujitsu FX100 (Nagoya Univ.)
# OAKLEAF : Oakleaf-fx, Fujitsu FX10 (Univ. Tokyo)
# K       : K computer, Fujitsu K computer (RIKEN)
HPC=


#### For IGS, set IGSDIR including METRIC*.OUT
#IGSDIR=/home/maeyama/gkvp/f0.52/AUGEM/aug_equilibrium/metriclist/igs_nss065ntheta049/field/
IGSDIR=


#### Set
if [[ -z ${HPC} ]]; then
  echo "Choose HPC system: HELIOS, PS, NU, OAKLEAF, K"
  exit

elif [[ ${HPC} = "HELIOS" ]]; then
  DIR=/csc/workdir3/smaeyama/gkvp/f0.52/rev13/
  LDM=gkvp_mpifft.exe
  NL=gkvp_f0.52_namelist.001
  SC=sbatch
  JS=sub.q

elif [[ ${HPC} = "PS" ]]; then
  DIR=/data/lng/maeyama/gkvp/f0.52/rev13/
  LDM=gkvp_mpifft.exe
  NL=gkvp_f0.52_namelist.001
  SC=pjsub
  JS=sub.q

elif [[ ${HPC} = "NU" ]]; then
  DIR=/large2/z41049a/gkvp/f0.52/rev13/
  LDM=gkvp_mpifft.exe
  NL=gkvp_f0.52_namelist.001
  SC=pjsub
  JS=sub.q

elif [[ ${HPC} = "OAKLEAF" ]]; then
  DIR=/group2/gc15/c15064/gkvp/f0.52/rev13/
  LDM=gkvp_mpifft.exe
  NL=gkvp_f0.52_namelist.001
  SC=pjsub
  JS=sub.q

elif [[ ${HPC} = "K" ]]; then
  DIR=/data/hp120011/maeyama/gkvp/f0.52/rev13/
  LDM=gkvp_mpifft.exe
  NL=gkvp_f0.52_namelist.001
  SC=pjsub
  JS=sub.q
  #JS=sub.q_cnt

else
  echo "The HPC system is not available. HPC = "${HPC}
fi

mkdir -p ${DIR}/
mkdir -p ${DIR}/log
mkdir -p ${DIR}/hst
mkdir -p ${DIR}/phi
mkdir -p ${DIR}/fxv
mkdir -p ${DIR}/cnt

cp ./${LDM} ${DIR}
cp ./${NL} ${DIR}
cp ./${JS} ${DIR}
cp ./Makefile ${DIR}
cp -r ../src ${DIR}
cp -r ../lib ${DIR}

if [[ -n ${IGSDIR} ]]; then
  mkdir -p ${DIR}/eqdsk
  cp ${IGSDIR}/METRIC*.OUT ${DIR}/eqdsk
fi

${SC} ${JS}
