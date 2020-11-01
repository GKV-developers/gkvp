#!/bin/sh

CURRENTDIR=`pwd`


  ### copy run dir ###                                                                                                ### for debug ###
  OLD_RUNDIR="shat0286beta0250"
  NEW_RUNDIR="shat0786beta0100";                                                                                      echo ${NEW_RUNDIR}
  mkdir -p ${NEW_RUNDIR}/run
  cp ${OLD_RUNDIR}/run/gkvp_mpifft.exe         ${NEW_RUNDIR}/run
  cp ${OLD_RUNDIR}/run/go.sh                   ${NEW_RUNDIR}/run
  cp ${OLD_RUNDIR}/run/sub.q                   ${NEW_RUNDIR}/run
  cp ${OLD_RUNDIR}/run/gkvp_f0.52_namelist.001 ${NEW_RUNDIR}/run
  
  ### change work dir name ###
  OLD_WKDIR="\/large2\/z41049a\/gkvp\/f0.52\/applegate2007ppcf_salpha_lb\/linky0070_shatbetadep\/shat0286beta0250"
  NEW_WKDIR="\/large2\/z41049a\/gkvp\/f0.52\/applegate2007ppcf_salpha_lb\/linky0070_shatbetadep\/shat0786beta0100";   #echo ${NEW_WKDIR}
  sed -i -e "s/${OLD_WKDIR}/${NEW_WKDIR}/g" ${NEW_RUNDIR}/run/go.sh
  sed -i -e "s/${OLD_WKDIR}/${NEW_WKDIR}/g" ${NEW_RUNDIR}/run/sub.q
  sed -i -e "s/${OLD_WKDIR}/${NEW_WKDIR}/g" ${NEW_RUNDIR}/run/gkvp_f0.52_namelist.001
  
  ### change parameter ###
  OLD_PARAM1="s_hat    = 0.286d0"
  NEW_PARAM1="s_hat    = 0.786d0";                                                                                    #echo ${NEW_PARAM1}
  sed -i -e "s/${OLD_PARAM1}/${NEW_PARAM1}/g" ${NEW_RUNDIR}/run/gkvp_f0.52_namelist.001
  
  ### change parameter ###
  OLD_PARAM2="beta = 2.5d-2"
  NEW_PARAM2="beta = 0.01d0";                                                                                         #echo ${NEW_PARAM2}
  sed -i -e "s/${OLD_PARAM2}/${NEW_PARAM2}/g" ${NEW_RUNDIR}/run/gkvp_f0.52_namelist.001

  ### execution ###                                                                                                   #################
  cd ${NEW_RUNDIR}/run
  #./go.sh


cd ${CURRENTDIR}
