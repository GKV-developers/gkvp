#!/bin/awk
# <Note>
#   Sum of 2D data over the second column.
#   The 2D data, wk[1:NR,1:NF], consists of x[0:nx-1], y[0:ny-1], data[0:nx-1,0:ny-1,0:nd-1], comment and blank lines.
#
# <How to use>
#   awk -f reduction.awk datafile > out.dat
#
# <Example of datafile>
#       # Comment                                                        <- No blank line between comment and data.
#         x[   0] y[   0] data[   0,   0,   0] ... data[   0,   0,nd-1]
#         x[   0] y[   0] data[   1,   0,   0] ... data[   1,   0,nd-1]
#           |     |
#         x[nx-1] y[   0] data[nx-1,   0,   0] ... data[nx-1,   0,nd-1]
#                                                                        <- A blank line between data blocks.
#         x[   0] y[   1] data[   0,   1,   0] ... data[   0,   1,nd-1]
#         x[   1] y[   1] data[   1,   1,   0] ... data[   1,   1,nd-1]
#           |     |
#         x[nx-1] y[   1] data[nx-1,   1,   0] ... data[nx-1,   1,nd-1]
#
#                          |
#
#         x[   0] y[ny-1] data[   0,ny-1,   0] ... data[   0,ny-1,nd-1]
#         x[   1] y[ny-1] data[   1,ny-1,   0] ... data[   1,ny-1,nd-1]
#           |     |
#         x[nx-1] y[ny-1] data[nx-1,ny-1,   0] ... data[nx-1,ny-1,nd-1]
#                                                                        <- A blank line at the end.
#

BEGIN{
#  ncom = -1         # Number of comment lines
#  nx = -1           # Number of the first column
#  ny = -1           # Number of the second column
#  nd = -1           # Number of data
  for ( ix=0; ix<=nx-1; ix++ ) {
    for ( id=0; id<=nd-1; id++ ) {
      sum[ix,id] = 0.0
    }
  }
}

(NR<ncom){ print $0 }
(NR==ncom){ printf("%s%16s",$1,$2); for(ic=4;ic<NF;ic++){printf("%17s",$ic)}; printf("%17s\n",$NF) }
(ncom<NR)&&((NR-ncom-1)<nx){ x[NR-ncom-1] = $1 }
(ncom<NR)&&((NR-ncom-1)%(nx+1)==0){ y[int((NR-ncom-1)/(nx+1))] = $2 }
(ncom<NR)&&((NR-ncom-1)%(nx+1)!=nx){
  for ( id=0; id<=nd-1; id++ ) {
  #  data[(NR-ncom-1)%(nx+1),int((NR-ncom-1)/(nx+1)),id] = $(id+1+2)
    sum[(NR-ncom-1)%(nx+1),id] = sum[(NR-ncom-1)%(nx+1),id] + $(id+1+2)
  }
}

END{
  for ( ix=0; ix<=nx-1; ix++ ) {
    printf( "%17.7e", x[ix] )
    for ( id=0; id<=nd-1; id++ ) {
      printf( "%17.7e", sum[ix,id] )
    #  printf( "%17.7e", sum[ix,id]/ny ) # Average
    }
    printf( "\n" )
  }
}
