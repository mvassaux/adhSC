import os,sys

import numpy as np
import math
import subprocess
import fileinput
import re

# Setting up the class used to store the bodies data
class body:
     def __init__(self, upd, x1, x2, x3):
         self.upd = upd
         self.x1 = x1
         self.x2 = x2
         self.x3 = x3

# Handling the VLOC_RLOC file
VDAT = open("./Vloc_Rloc.LAST", "r")
VDATl = VDAT.readlines()
VDAT.close
VDATnl = len(VDATl)
## Reinitialization of the time
Vline = VDATl[3].split()
VDATl[3] = '$steps      %s              time=%14.7E  \n' % (0, 0.)
## Reinitialization of the interactions
for i in range(0,VDATnl,1):
   if '$icdan' in VDATl[i]:
      Vline = VDATl[i+2].split()
      if ('PTPT3' in VDATl[i]):
         mgaptt = re.search('gapTT(.{14})', VDATl[i+5])
         if mgaptt:
            gapTT = mgaptt.group(1)
         VDATl[i+3] = '                             rlt/H 0.0000000D+00  rln/H 0.0000000D+00  rls/H 0.0000000D+00 \n'
         VDATl[i+4] = '                             vlt = 0.0000000D+00  vln = 0.0000000D+00  vls = 0.0000000D+00 \n'
         VDATl[i+5] = '                             gapTT 0.0000000D+00 \n'
         VDATl[i+6] = '                             n(1)= 0.0000000D+00  n(2)= 0.0000000D+00  n(3)= 0.0000000D+00 \n'
         VDATl[i+9] = ' %s \n'% (gapTT)
      else:
         for j in range(0,11,1):
            VDATl[i+j] = ''
## Writing to the VLOC_RLOC initial file for the next phase
V = open("./Vloc_Rloc.INI", "w")
for i in range(0,VDATnl,1):
   V.write(VDATl[i])
V.close()

## Reinitialization of the time
DDAT = open("./DOF.LAST", "r")
DDATl = DDAT.readlines()
DDAT.close
DDATnl = len(DDATl)
## Reinitialization of the interactions
Dline = DDATl[3].split()
DDATl[3] = '$steps     %s              time=%14.7E  \n' % (0, 0.)
## Writing to the DOF initial file for the next phase
D = open("./DOF.INI", "w")
for i in range(0,6,1):
   D.write(DDATl[i])
D.close()
