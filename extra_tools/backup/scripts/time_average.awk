#!/bin/awk
# <Note>
#   Time average of 1D data
#   The 1D data, data[1:NR,1:NF], consists of time[1:NR], f[1:NR,2:NF] and comment or blank lines.
#
# <How to use>
#   awk -f time_average.awk -v tsta=*** tend=*** datafile > out.dat
#
# <Example of datafile>
#   # It may contain comment or blank lines.
#         t[ 1] f[ 1, 2] ... f[ 1,NF]
#         t[ 2] f[ 2, 2] ... f[ 2,NF]
#           |  
#         t[NR] f[NR, 2] ... f[NR,NF]
#

BEGIN{
#  tsta = 100  # Time averaged duration
#  tend = 200  #        tsta < t < tend
  iflag = 0
}

(tsta<$1)&&($1<tend){

  ### Read data during tsta < t < tend ###
  for ( ic=1; ic<=NF; ic++ ) {
    data[ic] = $ic
  }

  if (iflag==0){ ### Initialization ###
    tsta_a = data[1]
    for ( ic=1; ic<=NF; ic++ ) {
      sum[ic] = 0.0
      data0[ic] = data[ic]
    }
    iflag = 1

  } else { ### Time integration by trapezoid rule ###
    tend_a = data[1]
    dt = data[1] - data0[1]
    sum[1] = sum[1] + dt
    for ( ic=2; ic<=NF; ic++ ) {
      sum[ic] = sum[ic] + 0.5 * (data[ic] + data0[ic]) * dt
    }
    for ( ic=1; ic<=NF; ic++ ) {
      data0[ic] = data[ic]
    }

  }
}

END{ ### Output ###
  print "# Column, Averaged(" tsta_a "<t<" tend_a ")"
  for ( ic=2; ic<=NF; ic++ ) {
    printf( "%8i%17.7e\n", ic, sum[ic]/sum[1] )
  }
}

