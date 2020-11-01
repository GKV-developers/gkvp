#!/bin/sh

CURRENTDIR=`pwd`

param="1.0"

OLD_RUNDIR=`echo ${param} | awk '{printf("%s%05.2f","ky",$1)}'`
OLD_WKDIR="\/data\/group1\/z43460z\/gkvp\/f0.57\/bpp\/lin02_full\/${OLD_RUNDIR}"
kymin=`echo ${param} | awk '{printf("%5.2fd0",$1)}'`
OLD_PARAM="kymin = ${kymin},"

grep "$OLD_WKDIR" shoot &&\
grep "$OLD_PARAM" gkvp_namelist &&\
for param in `seq 0.4 0.1 1.0`; do # ky scan  (ky = param)

 
  ### change working dir ###                                                                           ### for debug ###
  NEW_RUNDIR=`echo ${param} | awk '{printf("%s%05.2f","ky",$1)}'`;                                     echo ${NEW_RUNDIR}
  NEW_WKDIR="\/data\/group1\/z43460z\/gkvp\/f0.57\/bpp\/lin02_full\/${NEW_RUNDIR}"
  sed -i -e "s/${OLD_WKDIR}/${NEW_WKDIR}/g" shoot
 
  ### change parameter ###
  kymin=`echo ${param} | awk '{printf("%5.2fd0",$1)}'`
  NEW_PARAM="kymin = ${kymin},";                                                                       echo ${NEW_PARAM}
  sed -i -e "s/${OLD_PARAM}/${NEW_PARAM}/g" gkvp_namelist
 
  ### execution ###                                                                                    #################
  ./shoot 13 14

  OLD_WKDIR=${NEW_WKDIR}
  OLD_PARAM=${NEW_PARAM}
  

done
