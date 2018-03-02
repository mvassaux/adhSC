import os,sys

import re
import numpy as np
import math
import subprocess
import fileinput

# Setting up the class used to store a focal adhesion data
class body:
     def __init__(self, num, ux, uy, uz, fx, fy, fz, fr, tr):
         self.num = num
         self.ux = ux
         self.uy = uy
         self.uz = uz
         self.fx = fx
         self.fy = fy
         self.fz = fz
         self.fr = fr
         self.tr = tr

# Setting the timestep to analyze
sRTS = sys.argv[1]
iRTS = int(sRTS)

# Setting output file name
outfile = sys.argv[2]

# Opening the file containing focal adhesions location
BDAT = open("./DATBOX/BODIES.DAT", "r")
BDATl = BDAT.readlines()
BDAT.close
BDATnl = len(BDATl)

# Counting and listing the focal adhesions
Npf = 0
Lpf = []

for i in range(0,BDATnl,1):
   if '$bdyty' in BDATl[i]:
      if 'INTEs' in BDATl[i+9]:
         Bline1 = BDATl[i+1].split()
         bpf = int(Bline1[1])
         Lpf.append(bpf)
         Npf = Npf + 1

## Initializing the structure to store focal adhesion data
pf = np.empty(Npf+1, dtype=object)
for i in range(0,Npf,1):
   pf[i] = body(0,0.,0.,0.,0.,0.,0.,0.,0.)

# Retrieving the force vector on each focal adhesion
for i in range(0,Npf,1):
    pf[i].num = Lpf[i]
    name='./POSTPRO/REAC_%07d.DAT'%(Lpf[i])

    RDAT = open(name, "r")
    RDATl = RDAT.readlines()
    RDAT.close
    RDATnl = len(RDATl)

    Rline = RDATl[iRTS-1].split()
    pf[i].fx = -1*float(Rline[1].replace('D','E'))
    pf[i].fy = -1*float(Rline[2].replace('D','E'))
    pf[i].fz = -1*float(Rline[3].replace('D','E'))
    pf[i].fr = ((pf[i].fx)**2+(pf[i].fy)**2+(pf[i].fz)**2)**0.5

# Computing the total resulting forces
Rtot = 0.
Rhtot = 0.
Ttot = 0.
for i in range(0,Npf,1):
   Rtot = Rtot + (((pf[i].fx)**2+(pf[i].fy)**2+(pf[i].fz)**2)**0.5)
   Rhtot = Rhtot + (((pf[i].fx)**2+(pf[i].fy)**2)**0.5)
   Ttot = Ttot + pf[i].tr

# Printing the results to be picked up by the calling bash script
# print(iRTS,",",Rtot,",",Rhtot)
with open(outfile, "a") as myfile:
    myfile.write(str(iRTS)+","+str(Rtot)+","+str(Rhtot)+"\n")
