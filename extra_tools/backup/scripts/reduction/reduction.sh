#!/bin/sh

#DIAGparameters
snum=1
enum=1
loopskip=10
flag_cnt=F
flag_ffinzv=F
flag_mom=T
flag_mominkxky=T
flag_mominxy=F
flag_mominz=F
flag_trninkxky=F

#GKVparameters
nxw=128
nyw=32
nx=80
global_ny=19
global_nz=32
global_nv=48
global_nm=15
nzb=2
nprocw=4
nprocz=8
nprocv=12
nprocm=2
nprocs=2
vmax=4.00000000000000 


date

### Reduction. INPUT:mominkxky_t*, OUTPUT:mominkx_t* ###
loopsta=0
loopend=1800

  echo "Reduction. INPUT:mominkxky_t\*, OUTPUT:mominkx_t\*"
  echo "loopsta = ${loopsta}, loopend = ${loopend}"
  loop=${loopsta}
  while [ $loop -le ${loopend} ]
  do
    cloop=`printf "%08d" $loop`
    awk -f ./reduction_xy2x.awk \
        -v ncom=2 nx=`expr 2 \* ${nx} + 1` ny=`expr ${global_ny} + 1` nd=`expr 2 + 8 \* ${nprocs}` \
        ./plt/mominkxky_t${cloop}.dat > ./plt/mominkx_t${cloop}.dat
    loop=`expr ${loop} + ${loopskip}`
  done


### Reduction. INPUT:mominkxky_t*, OUTPUT:mominky_t* ###
loopsta=0
loopend=1800

  echo "Reduction. INPUT:mominkxky_t\*, OUTPUT:mominky_t\*"
  echo "loopsta = ${loopsta}, loopend = ${loopend}"
  loop=${loopsta}
  while [ $loop -le ${loopend} ]
  do
    cloop=`printf "%08d" $loop`
    awk -f ./reduction_xy2y.awk \
        -v ncom=2 nx=`expr 2 \* ${nx} + 1` ny=`expr ${global_ny} + 1` nd=`expr 2 + 8 \* ${nprocs}` \
        ./plt/mominkxky_t${cloop}.dat > ./plt/mominky_t${cloop}.dat
    loop=`expr ${loop} + ${loopskip}`
  done


### Time average. INPUT:mominkxky_t*, OUTPUT:ave_mominkxky.dat ###
loopsta=1000
loopend=1800

  echo "Time average. INPUT:mominkxky_t*, OUTPUT:ave_mominkxky.dat"
  echo "loopsta = ${loopsta}, loopend = ${loopend}"
  loop=${loopsta}
  filename=
  filenum=
  while [ $loop -le ${loopend} ]
  do
    cloop=`printf "%08d" $loop`
    filename=${filename}" ./plt/mominkxky_t${cloop}.dat"
    filenum=`expr ${filenum} + 1`
    loop=`expr ${loop} + ${loopskip}`
  done
 # echo "Number of files = ${filenum}"
 # echo ${filename}
  awk -f ./timeaverage.awk \
      -v ncom=2 nrow=`expr \( 2 \* ${nx} + 1 + 1 \) \* \( ${global_ny} + 1 \)` ncol=`expr 2 + 2 + 8 \* ${nprocs}` \
      nfile=${filenum} ${filename} > ./plt/ave_mominkxky.dat
     # nfile=`ls -1 ./plt/mominkxky_t* | wc -l` ./plt/mominkxky_t* > ./plt/ave_mominkxky.dat


### Time average. INPUT:mominkx_t*, OUTPUT:ave_mominkx.dat ###
loopsta=1000
loopend=1800

  echo "Time average. INPUT:mominkx_t*, OUTPUT:ave_mominkx.dat"
  echo "loopsta = ${loopsta}, loopend = ${loopend}"
  loop=${loopsta}
  filename=
  filenum=
  while [ $loop -le ${loopend} ]
  do
    cloop=`printf "%08d" $loop`
    filename=${filename}" ./plt/mominkx_t${cloop}.dat"
    filenum=`expr ${filenum} + 1`
    loop=`expr ${loop} + ${loopskip}`
  done
 # echo "Number of files = ${filenum}"
 # echo ${filename}
  awk -f ./timeaverage.awk \
      -v ncom=2 nrow=`expr 2 \* ${nx} + 1` ncol=`expr 1 + 2 + 8 \* ${nprocs}` \
      nfile=${filenum} ${filename} > ./plt/ave_mominkx.dat
 # awk -f ./reduction_xy2x.awk \
 #     -v ncom=0 nx=`expr 2 \* ${nx} + 1` ny=`expr ${global_ny} + 1` nd=`expr 2 + 8 \* ${nprocs}` \
 #     ./plt/ave_mominkxky.dat > ./plt/ave_mominkx.dat


### Time average. INPUT:mominky_t*, OUTPUT:ave_mominky.dat ###
loopsta=1000
loopend=1800

  echo "Time average. INPUT:mominky_t*, OUTPUT:ave_mominky.dat"
  echo "loopsta = ${loopsta}, loopend = ${loopend}"
  loop=${loopsta}
  filename=
  filenum=
  while [ $loop -le ${loopend} ]
  do
    cloop=`printf "%08d" $loop`
    filename=${filename}" ./plt/mominky_t${cloop}.dat"
    filenum=`expr ${filenum} + 1`
    loop=`expr ${loop} + ${loopskip}`
  done
 # echo "Number of files = ${filenum}"
 # echo ${filename}
  awk -f ./timeaverage.awk \
      -v ncom=2 nrow=`expr ${global_ny} + 1` ncol=`expr 1 + 2 + 8 \* ${nprocs}` \
      nfile=${filenum} ${filename} > ./plt/ave_mominky.dat
 # awk -f ./reduction_xy2y.awk \
 #     -v ncom=0 nx=`expr 2 \* ${nx} + 1` ny=`expr ${global_ny} + 1` nd=`expr 2 + 8 \* ${nprocs}` \
 #     ./plt/ave_mominkxky.dat > ./plt/ave_mominky.dat

date

