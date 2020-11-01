#!/usr/bin/env python3
nprocw = 4
nprocz = 5
nprocv = 4
nprocm = 4
nprocs = 4

ll  = 4   # process number per node
ls1 = 1   # segment length in z
ls2 = 1   # segment length in v
ls3 = 1   # segment length in m
lg1 = 5   # number of segments in z
lg2 = 4   # number of segments in v
lg3 = 4   # number of segments in m
lp2 = 2   # number of species in v
lp3 = 2   # number of species in m

#------------------------------------------------------------------------------
#
#   Create rank-map file for 3D torus network
#
#
#                   -----      
#                 /     /      
#                 -----  | ls3   "a segment"  -  Parallel FFT is performed    
#                |    |  |                       in a segment.
#                |    | / ls2                    ( ll*ls1*ls2*ls3 == nprocw )
#                 -----          
#                  ls1         
#
#
#
#           -----    -----    ------       
#         /     /  /     /  /      /       
#         -----  | -----  | ------  |      
#        |     | ||     | ||     |  |      
#      -----    -----    -----   | /       "a segment group for each species"
#    /     /  /     /  /     /  --    lg3  
#    -----  | -----  | -----  |    /          -  Segments are set up so that
#   |    |  ||    |  ||    |  | --  |            communications in zz, vl, mu
#   |    | / |    | / |    | /   |  |            are performed between
#    -----    -----    -----     | /             adjacent segments.
#    /     /  /     /  /     /  --               ( lg1 == nprocz,
#    -----  | -----  | -----  |                    lg2 == nprocv,
#   |    |  ||    |  ||    |  |   lg2              lg3 == nprocm  )
#   |    | / |    | / |    | /             
#    -----    -----    -----               
#              lg1                         
#
#
#
#          ------------------------- 
#        /                         / 
#       /                         / |    
#      /                         /  |     "species set up"                      
#     /                         /   |                                           
#    /                         /   /         -  Segment groups are arranged
#    -------------------------    /             in vl (or mu) direction
#   |                         |  / /            so that spc_comm_world is
#   |                         | / / |           performed in a cross section. 
#   |                         |/ /  |           ( lp2*lp3 == nprocs )
#    -------------------------  /   |
#    /                         /   /
#    -------------------------    /
#   |                         |  /
#   |                         | /
#   |                         |/
#    ------------------------- 
#
#
#
#------------------------------------------------------------------------------
import sys

if (ll*ls1*ls2*ls3*lg1*lg2*lg3*lp2*lp3 != nprocw*nprocz*nprocv*nprocm*nprocs):
    print("Total process number is illegal.")
    sys.exit()
elif (ll*ls1*ls2*ls3 != nprocw):
    print("Segment size is illegal.")
    sys.exit()
elif (lg1 != nprocz):
    print("Number of segments is illegal.")
    sys.exit()
elif (lg2 != nprocv):
    print("Number of segments is illegal.")
    sys.exit()
elif (lg3 != nprocm):
    print("Number of segments is illegal.")
    sys.exit()
elif (lp2*lp3 != nprocs):
    print("Number of species is illegal.")
    sys.exit()


mapping=[]
with open("for_gnuplot.dat", mode="w") as f:
    rank = 0
    for p in range(lp3):    # mapping for each species
        for o in range(lp2):
  
           for n in range(lg3):    # mapping of segments
               for m in range(lg2):
                   for l in range(lg1):
  
                       for k in range(ls3):    # mapping in a segment
                           for j in range(ls2):
                               for i in range(ls1):
                                   for h in range(ll):    # mapping in a node
                                       p1 = i+ls1*l
                                       p2 = j+ls2*m+ls2*lg2*o
                                       p3 = k+ls3*n+ls3*lg3*p
                                       mapping.append([p1,p2,p3])
                                       f.write("{} {} {} # rank = {}\n".format(p1, p2, p3, rank))
                                       rank = rank + 1
                               #f.write("\n")
                           f.write("\n")
                       f.write("\n")
  
                   f.write("\n")
               f.write("\n")
           f.write("\n")
  
#import numpy as np
#print(np.array(mapping).shape)

with open("rankmapfile.dat", mode="w") as f:
    for rank in range(len(mapping)):
       f.write("({},{},{})\n".format(mapping[rank][0],mapping[rank][1],mapping[rank][2]))
