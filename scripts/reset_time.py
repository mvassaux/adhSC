import os,sys

import numpy as np
import math
import subprocess
import fileinput
import re

# Handling the VLOC_RLOC file
VDAT = open("./Vloc_Rloc.LAST", "r")
VDATl = VDAT.readlines()
VDAT.close
VDATnl = len(VDATl)
## Reinitialization of the time
Vline = VDATl[3].split()
VDATl[3] = '$steps      %s              time=%14.7E  \n' % ('0', 0.)
## Saving file
V = open("./Vloc_Rloc.INI", "w")
for i in range(0,VDATnl,1):
   V.write(VDATl[i])
V.close()

# Handling the DOF file
DDAT = open("./DOF.LAST", "r")
DDATl = DDAT.readlines()
DDAT.close
DDATnl = len(DDATl)
## Reinitialization of the time
Dline = DDATl[3].split()
DDATl[3] = '$steps     %s              time=%14.7E  \n' % ('0', 0.)
## Saving file
D = open("./DOF.INI", "w")
for i in range(0,DDATnl,1):
   D.write(DDATl[i])
D.close()
