# Running the simulation of one of the adhesion phases based on the configuration
# provided in the DATBOX directory

# importing chipy module
from pylmgc90.chipy import *
import numpy

# Initializing
#Initialize()

def compute_poteng(field1, field2):
    W=0.
    for ibdyty in range(RBDY3_GetNbRBDY3()):
        ff1 = RBDY3_GetBodyVector(field1,ibdyty+1)
        ff2 = RBDY3_GetBodyVector(field2,ibdyty+1)
        for ii in range(len(ff1)):
            W += -0.5*ff1[ii]*ff2[ii]

    print(field1,"x",field2,W)
    return W

# checking/creating mandatory subfolders
checkDirectories()

# logMes
# utilities_DisableLogMes()

#
# defining some variables
#

# space dimension
dim = 3

# modeling hypothesis ( 1 = plain strain, 2 = plain stress, 3 = axi-symmetry)
mhyp = 1

# time evolution parameters
nb_steps = int(sys.argv[1])
dt = float(sys.argv[2])

# theta integrator parameter
theta = 0.5

# deformable  yes=1, no=0
deformable = 0

# interaction parameters
freq_detect = 1
Rloc_tol = 5.e-2

# nlgs parameters
tol = 1e-4
relax = 1.0
norm = 'Quad '
gs_it1 = 50
gs_it2 = 10
type='Stored_Delassus_Loops         '

# write parameter
noutp = int(sys.argv[3])
freq_write   =  max(1,int(nb_steps/noutp))

# display parameters
ndisp = int(sys.argv[4])
freq_display =  max(1,int(nb_steps/ndisp))
diam_cell = 30.
lref = diam_cell/60

#
# read and load
#

# Set space dimension
SetDimension(dim,mhyp)
#
utilities_logMes('INIT TIME STEPPING')
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)
#
utilities_logMes('READ BEHAVIOURS')
ReadBehaviours()
if deformable: ReadModels()
#
utilities_logMes('READ BODIES')
ReadBodies()
#
utilities_logMes('LOAD BEHAVIOURS')
LoadBehaviours()
if deformable: LoadModels()
#
utilities_logMes('READ INI DOF')
ReadIniDof()
#
if deformable:
  utilities_logMes('READ INI GPV')
  ReadIniGPV()
#
utilities_logMes('READ DRIVEN DOF')
ReadDrivenDof()
#
utilities_logMes('LOAD TACTORS')
LoadTactors()
#
utilities_logMes('READ INI Vloc Rloc')
ReadIniVlocRloc()

#
# paranoid writes
#
utilities_logMes('WRITE BODIES')
WriteBodies()
utilities_logMes('WRITE BEHAVIOURS')
WriteBehaviours()
utilities_logMes('WRITE DRIVEN DOF')
WriteDrivenDof()

#
# open display & postpro
#

utilities_logMes('DISPLAY & WRITE')
OpenDisplayFiles()
OpenPostproFiles()
PT3Dx_SetReferenceRadius(lref/10.)


#
# simulation part ...
#

# ... calls a simulation time loop
# since constant compute elementary mass once
ComputeMass()
# set gravity acceleration
grav=numpy.zeros(3)
grav[2]=-9.81*1e-6
bulk_behav_SetGravity(grav)

for k in range(0,nb_steps):
  #
  IncrementStep()

  ComputeFext()
  ComputeBulk()
  ComputeFreeVelocity()

  SelectProxTactors(freq_detect)

#  RecupRloc(Rloc_tol)

  ExSolver(type, norm, tol, relax, gs_it1, gs_it2)
  UpdateTactBehav()

  StockRloc()

  ComputeDof()
  #RBDY3_FatalDamping(1)
  UpdateStep()

  WriteOutDof(freq_write)
  WriteOutVlocRloc(freq_write)

  WriteDisplayFiles(freq_display,lref)
  WritePostproFiles()

  # Check energies in the system
  KE = postpro_3D_GetKineticEnergy()
  print(KE)

  # Computes PE = -0.5*field1*field2
  PE = compute_poteng("X____","Reac_")

  # if PE>PE_tsh:
  #     break

WriteLastDof()
WriteLastVlocRloc()
#
# close display & postpro
#
CloseDisplayFiles()
ClosePostproFiles()

# this is the end
Finalize()
