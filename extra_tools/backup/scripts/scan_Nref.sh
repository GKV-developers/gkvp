#!/bin/sh

CURRENTDIR=`pwd`

param="9.0"

OLD_RUNDIR=`echo ${param} | awk '{printf("%s%05.2fd20","Nref",$1)}'`
OLD_WKDIR="\/data\/group1\/z43460z\/gkvp\/f0.57\/bpp\/lin04_full_collscan\/${OLD_RUNDIR}"
Nref=`echo ${param} | awk '{printf("%5.2fd20",$1)}'`
OLD_PARAM="Nref = ${Nref},"

grep "$OLD_WKDIR" shoot &&\
grep "$OLD_PARAM" gkvp_namelist &&\
for param in `seq 3.0 1.0 9.0`; do # Nref scan  (Nref = param)

 
  ### change working dir ###                                                                           ### for debug ###
  NEW_RUNDIR=`echo ${param} | awk '{printf("%s%05.2fd20","Nref",$1)}'`;                                echo ${NEW_RUNDIR}
  NEW_WKDIR="\/data\/group1\/z43460z\/gkvp\/f0.57\/bpp\/lin04_full_collscan\/${NEW_RUNDIR}"
  sed -i -e "s/${OLD_WKDIR}/${NEW_WKDIR}/g" shoot
 
  ### change parameter ###
  Nref=`echo ${param} | awk '{printf("%5.2fd20",$1)}'`
  NEW_PARAM="Nref = ${Nref},";                                                                         echo ${NEW_PARAM}
  sed -i -e "s/${OLD_PARAM}/${NEW_PARAM}/g" gkvp_namelist
 
  ### execution ###                                                                                    #################
  ./shoot 8 8

  OLD_WKDIR=${NEW_WKDIR}
  OLD_PARAM=${NEW_PARAM}
  

done
