#!/bin/awk
# <Note>
#   Average of data over input files.
#   Al file must have the same form. The data consists of data[1:nrow,1:ncol], comment and blank lines.
#
# <How to use>
#   awk -f reduction.awk datafile > out.dat
#
# <Example of datafile>
#       # Comment
#         data[   1,   1] data[   1,   2] ... data[   1,ncol]
#         data[   2,   1] data[   2,   2] ... data[   2,ncol]
#                          |
#         data[  ir,   1] data[  ir,   2] ... data[  ir,ncol]
#                                                                   <- Blank lines may appear.
#                          |
#         data[nrow,   1] data[nrow,   2] ... data[nrow,ncol]
#

BEGIN{
#  ncom = -1         # Number of comment lines
#  nfile = -1        # Number of files
#  nrow = -1         # Number of rows ( including blank lines, but not comment lines )
#  ncol = -1         # Number of columns
  for ( ir=1; ir<=nrow; ir++ ) {
    for ( ic=1; ic<=ncol; ic++ ) {
      sum[ir,ic] = 0.0
    }
  }
}

#(FNR<=ncom){ print $0 }
(ncom<FNR){
  if ( NF==0 ) {
    sum[FNR-ncom,1] = "blankline"
  } else {
    for ( ic=1; ic<=ncol; ic++ ) {
    #  data[FNR-ncom,ic] = $ic
      sum[FNR-ncom,ic] = sum[FNR-ncom,ic] + $ic
    }
  }
}

END{
  for ( ir=1; ir<=nrow; ir++ ) {
    if ( sum[ir,1]=="blankline" ) {
      print
    } else {
      for ( ic=1; ic<=ncol; ic++ ) {
      #  printf( "%17.7e", sum[ir,ic] )
        printf( "%17.7e", sum[ir,ic]/nfile ) # Average
      }
      printf( "\n" )
    }
  }
}
