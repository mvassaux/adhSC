import os,sys

import re
import numpy as np
import math
import subprocess
import fileinput

# Setting up the class used to store the bodies data
class body:
     def __init__(self, ox, oy, oz, ux, uy, uz):
         self.ux = ox
         self.uy = oy
         self.uz = oz
         self.fx = ux
         self.fy = uy
         self.fz = uz

# Handling the BODIES file of the previous phase
BDAT = open("./BODIES.OLD", "r")
BDATl = BDAT.readlines()
BDAT.close
BDATnl = len(BDATl)

# Handling the DOF file at the end of the last phase
DDAT = open("./DOF.LAST", "r")
DDATl = DDAT.readlines()
DDAT.close
DDATnl = len(DDATl)

# Counting number of bodies/particles
nbpt = 0
for i in range(0,DDATnl,1):
   if '$bdyty' in DDATl[i]:
      nbpt = nbpt + 1

# Initilizating the storage of the bodies data
bod = np.empty(nbpt+1, dtype=object)
for i in range(0,nbpt+1,1):
   bod[i] = body(0.,0.,0.,0.,0.,0.)

# Sorting the data of the DOF file
j = 0
for i in range(0,DDATnl,1):
   if '$bdyty' in DDATl[i]:
      j = j + 1
      mx = re.search('X.1.=(.{14})', DDATl[i+3])
      if mx:
         nux = mx.group(1)
         bod[j].ux = float(nux.replace('D','E'))
      my = re.search('X.2.=(.{14})', DDATl[i+3])
      if my:
         nuy = my.group(1)
         bod[j].uy = float(nuy.replace('D','E'))
      mz = re.search('X.3.=(.{14})', DDATl[i+3])
      if mz:
         nuz = mz.group(1)
         bod[j].uz = float(nuz.replace('D','E'))

# Updating the initial position of the bodies with the previous phase
# displacements
j = 0
for i in range(0,BDATnl,1):
   if '$bdyty' in BDATl[i]:
      j = j + 1
      mx = re.search('coo1=(.{14})', BDATl[i+6])
      if mx:
         nox = mx.group(1)
         bod[j].ox = float(nox.replace('D','E')) + bod[j].ux
      my = re.search('coo2=(.{14})', BDATl[i+6])
      if my:
         noy = my.group(1)
         bod[j].oy = float(noy.replace('D','E')) + bod[j].uy
      mz = re.search('coo3=(.{14})', BDATl[i+6])
      if mz:
         noz = mz.group(1)
         bod[j].oz = float(noz.replace('D','E')) + bod[j].uz

      for k in range(1,j,1):
         while ((abs(bod[k].ox-bod[j].ox)<1.0e-6) and (abs(bod[k].oy-bod[j].oy)<1.0e-6) and (abs(bod[k].oz-bod[j].oz)<1.0e-6)):
             print('SUPERPOSED PARTICLES > Seprating with offset...', k, j)
             bod[j].ox = bod[j].ox + 1.0e-6
             bod[j].oy = bod[j].oy + 1.0e-6
             bod[j].oz = bod[j].oz + 1.0e-6

      BDATl[i+6] = ' NO6xx      1                coo1=%14.7E  coo2=%14.7E  coo3=%14.7E  \n' % (bod[j].ox, bod[j].oy, bod[j].oz)

# Writing the BODIES file for the next phase
B = open("./BODIES.DAT", "w")
for i in range(0,BDATnl,1):
   B.write(BDATl[i])
B.close()
