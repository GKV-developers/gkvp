#!/bin/awk
# <Note>
#   Sum of 1D data.
#   The 1D data, wk[1:NR,1:NF], consists of x[1:nx], data[0:nx-1,0:nd-1] and comment or blank lines.
#
# <How to use>
#   awk -f reduction.awk datafile > out.dat
#
# <Example of datafile>
#       # Comment                                      <- No blank line between comment and data.
#         x[   0] data[   0,   0] ... data[   0,nd-1]
#         x[   1] data[   1,   0] ... data[   1,nd-1]
#           |  
#         x[nx-1] data[nx-1,   0] ... data[nx-1,nd-1]
#

BEGIN{
#  ncom = -1  # Number of comment lines
#  nx = -1    # Number of the first column
#  nd = -1    # Number of data
  for ( id=0; id<=nd-1; id++ ) {
    sum[id] = 0.0
  }
}

(NR<ncom){ print $0 }
(NR==ncom){ printf("%s%16s",$1,$3); for(ic=4;ic<NF;ic++){printf("%17s",$ic)}; printf("%17s\n",$NF) }
(ncom<NR)&&((NR-ncom-1)<nx){ 
  x[NR-ncom-1] = $1
  for ( id=0; id<=nd-1; id++ ) {
  #  data[NR-ncom-1,id] = $(id+1+1)
    sum[id] = sum[id] + $(id+1+1)
  }
}

END{
  for ( id=0; id<=nd-1; id++ ) {
    printf( "%17.7e", sum[id] )
  #  printf( "%17.7e", sum[id]/nx ) # Average
  }
  printf( "\n" )
}
