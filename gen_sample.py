# Cell adhesion model based on the software LMGC90, coded in Python 3
# Developed by Maxime Vassaux (m.vassaux@ucl.ac.uk)

# Adhesion is computed in two stages:
#   - a spreading phase from a elliptical cell configuration
#   - a stabilisation phase where the real cell forces are computed

#Units: nanonewton = micrometre * picog * (microsecondss)^-2
#Units: nanonewton = micrometre * g * (secondss)^-2

# Importing required python general libraries
import os,sys
import numpy
import math
import random
import subprocess
import fileinput

# Importing the LMGC90 system generation python library
from pylmgc90.pre import *

# Setting the seeds of the random processes used for the generation
# of intracellular components
seed_fo = 12345
seed_mf = 12345
seed_mt = 12345
seed_if = 12345

# Creation of directories where built systems are stored
if not os.path.isdir('./DATBOX_SPRD'):
  os.mkdir('./DATBOX_SPRD')
if not os.path.isdir('./DATBOX_STBL'):
  os.mkdir('./DATBOX_STBL')

# Setting the simulation to 3D
dim = 3

# Containers creation
## for particles bodies
bodies1 = avatars()
bodies2 = avatars()
ns_ls = []
ce_ls = []
fa_ls = []
## for material types
mat = materials()
## for model types
mod = models()
## for vibility tables
svs = see_tables()
## for contact laws
tacts_stif = tact_behavs()
tacts_nrml = tact_behavs()
tactsfinal = tact_behavs()

# Processing the gen_sample script inputs
nsteps_sprd = int(sys.argv[1])
dt_sprd = float(sys.argv[2])
nsteps_stbl = int(sys.argv[3])
dt_stbl = float(sys.argv[4])
noutp = int(sys.argv[5])

# Computing the loading increment speed: 1/(nstep*dt) > seconds-1
charinc_sprd =1/(nsteps_sprd*dt_sprd)
charinc_stbl =1/(nsteps_stbl*dt_stbl)

# Defining material types
tdur = material(name='TDURx',materialType='RIGID',density=1.)
cell = material(name='CELLx',materialType='RIGID',density=1.0e-8) # densite eau (1e-3 pg/microm3) ou (1e-15 g/microm3)
mat.addMaterial(tdur,cell)

# Defining model types
modr = model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)
mod.addModel(modr)

# Defining geometrical parameters of the cell and constituents (micrometers)
## Cell dimensions
diam_cell = 30.
diam_core = diam_cell/1.5
lref = diam_cell/60.
xctr_cell = 0.
yctr_cell = 0.
zctr_cell = 0.
center_cell=numpy.array([xctr_cell, yctr_cell, zctr_cell])
ctsl_part_radii = 3.0*lref
celb_part_radii = 3.0*lref
celb_part_gap = 1.0*celb_part_radii
radii_celb = diam_cell/2. + celb_part_radii

## Nucleus dimensions
core_part_radii = 1.5*lref
core_part_gap = 2.0*core_part_radii
corb_part_radii = lref
corb_part_gap = 2.0*corb_part_radii
radii_corb = diam_core/2. + corb_part_radii

## Substrate dimensions
lxsubst = 200.
lysubst = 200.
lzsubst = 5.0*lref
subst_part_radii = lzsubst
subst_part_gap = 2.0*subst_part_radii
xctr_subst = 0.
yctr_subst = 0.
zctr_subst = -0.75*diam_cell
center_subst=numpy.array([xctr_subst, yctr_subst, zctr_subst])

## Filaments networks parameters
nb_nodes_mf=400
nb_nodes_mt=400
nb_nodes_if=200
radii_part_networks = 1.0*lref
alert_intra_network = 4.0*radii_part_networks
alert_extra_network = 2.0*radii_part_networks

## Stress fibres parameters
fibers_part_radii = 0.5*lref
sf_inc_min = 5.*fibers_part_radii
sf_inc_max = 10.*fibers_part_radii
shoot_cone_angle = math.pi/2.

## Focal adhesions parameters
shift_amp_prot = 3.0
nb_focals = 8
nb_focals_ptour = 10
nb_fibres = 0
coor_focals = numpy.zeros((3*nb_focals+1), dtype='f')
coor_protei = numpy.zeros((3*nb_focals+1), dtype='f')
seen_focals = numpy.zeros((nb_focals+1), dtype='i')
numfibre_focals = numpy.zeros((nb_focals+1), dtype='i')

# Generating nucleus membrane particles
nphiinc = int(math.pi*radii_corb/corb_part_gap)
phi=-1*math.pi/2.
for i in range(0,nphiinc,1):
   nthetainc = 2*int(math.pi*radii_corb*math.cos(phi)/corb_part_gap)
   theta=0.
   for j in range(0,nthetainc,1):
      body=rigidSphere(r=corb_part_radii,
            center=numpy.array([(xctr_cell+(radii_corb)*(math.cos(phi)*math.cos(theta))),
                                (yctr_cell+(radii_corb)*(math.cos(phi)*math.sin(theta))),
                                (zctr_cell+(radii_corb)*(math.sin(phi)))]),
            model=modr, material=cell, color='BNOYs')
      body.addContactors(shape='PT3Dx', color='BNOYp')
      bodies1 += body
      bodies2 += body
      ns_ls.append(body)
      ce_ls.append(body)
      theta += 2*math.pi/nthetainc
   phi += math.pi/nphiinc

# Generating inner nucleus particles
core_part_nb = 0
while (core_part_nb < 100):
   core_part_rlst=granulo_Random(10000, 0.75*core_part_radii, 1.25*core_part_radii)
   core_part_nb, core_part_coor=depositInBox3D(core_part_rlst, diam_core, diam_core, diam_core)
for i in range(0,core_part_nb,1):
   px = xctr_cell + core_part_coor[3*i + 0]
   py = yctr_cell + core_part_coor[3*i + 1]
   pz = zctr_cell + core_part_coor[3*i + 2] - diam_core/2.
   dist_cell_cntr = ((px-xctr_cell)**2.0+(py-yctr_cell)**2.0+(pz-zctr_cell)**2.0)**0.5
   if((dist_cell_cntr + core_part_rlst[i] < diam_core/2.)):
      body = rigidSphere(r=core_part_rlst[i], center=numpy.array([px, py, pz]),
                      model=modr, material=cell, color='INOYs')
      body.addContactors(shape='PT3Dx', color='INOYp')
      bodies1 += body
      bodies2 += body
      ns_ls.append(body)
      ce_ls.append(body)

# Generating cell membrane particles
nphiinc = int(math.pi*radii_celb/celb_part_gap)
phi=-1*math.pi/2.
for i in range(0,nphiinc,1):
   nthetainc = 2*int(math.pi*radii_celb*math.cos(phi)/celb_part_gap)
   theta=0.
   for j in range(0,nthetainc,1):
      if((zctr_cell+(radii_celb)*(math.sin(phi))) < zctr_cell):
         body=rigidSphere(r=celb_part_radii,
            center=numpy.array([(xctr_cell+(radii_celb)*(math.cos(phi)*math.cos(theta))),
                                (yctr_cell+(radii_celb)*(math.cos(phi)*math.sin(theta))),
                                (zctr_cell+(radii_celb)*(math.sin(phi)))]),
            model=modr, material=cell, color='BCELs')
         body.addContactors(shape='PT3Dx', color='BCELp')
         body.addContactors(shape='SPHER', color='BCELt', byrd=celb_part_radii/10.)
         bodies1 += body
         bodies2 += body
         ce_ls.append(body)
      theta += 2*math.pi/nthetainc
   phi += math.pi/nphiinc
nphiinc = int(math.pi*radii_celb/celb_part_gap)
phi=-1*math.pi/2.
for i in range(0,nphiinc,1):
   nthetainc = 2*int(math.pi*radii_celb*math.cos(phi)/celb_part_gap)
   theta=0.
   for j in range(0,nthetainc,1):
      if((zctr_cell+(radii_celb)*(math.sin(phi))) > zctr_cell):
         body=rigidSphere(r=celb_part_radii,
            center=numpy.array([(xctr_cell+(radii_celb)*(math.cos(phi)*math.cos(theta))),
                                (yctr_cell+(radii_celb)*(math.cos(phi)*math.sin(theta))),
                                (zctr_cell+(radii_celb)*(math.sin(phi)))]),
            model=modr, material=cell, color='BCELs')
         body.addContactors(shape='PT3Dx', color='BCELp')
         body.addContactors(shape='SPHER', color='BCELt', byrd=celb_part_radii/10.)
         bodies1 += body
         bodies2 += body
         ce_ls.append(body)
      theta += 2*math.pi/nthetainc
   phi += math.pi/nphiinc

# Generating cytosol particles
ctsl_part_nb = 0
while (ctsl_part_nb < 100):
   ctsl_part_rlst=granulo_Random(100000, 1.0*ctsl_part_radii, 2.0*ctsl_part_radii)
   ctsl_part_nb, ctsl_part_coor=depositInBox3D(ctsl_part_rlst, 3.*diam_cell, 3.*diam_cell, 1.*diam_cell)
for i in range(0,ctsl_part_nb,1):
   px = xctr_cell + ctsl_part_coor[3*i + 0]
   py = yctr_cell + ctsl_part_coor[3*i + 1]
   pz = zctr_cell + ctsl_part_coor[3*i + 2]
   dist_cell_cntr = ((px-xctr_cell)**2.0+(py-yctr_cell)**2.0+(pz-zctr_cell)**2.0)**0.5
   if((dist_cell_cntr - ctsl_part_rlst[i] > diam_core/2. + core_part_radii) and
      (dist_cell_cntr + ctsl_part_rlst[i] < diam_cell/2.)):
      body = rigidSphere(r=ctsl_part_rlst[i], center=numpy.array([px, py, pz]),
                      model=modr, material=cell, color='CTSLs')
      body.addContactors(shape='PT3Dx', color='CTSLp')
      bodies1 += body
      bodies2 += body
      ce_ls.append(body)

#Proposition de modif du gen sampl pour deformer la cellule en cube : 

phi=math.pi/4.
theta=math.pi/2.
rad=radii_celb + celb_part_gap

i=-1
j=1
k=0
coor_focals[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j)
coor_protei[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j) + ((rad*(2./3.)**0.5*math.cos(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.cos(theta*j))*2.*celb_part_gap 
coor_focals[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j)
coor_protei[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j) + ((rad*(2./3.)**0.5*math.sin(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.sin(theta*j))*2.*celb_part_gap
coor_focals[3*k + 2] = rad*-1.*(1./3.)**0.5
coor_protei[3*k + 2] = rad*-1.*(1./3.)**0.5 + ((rad*-1.*(1./3.)**0.5)**2.0)**0.5/(rad*-1.*(1./3.)**0.5)*2.*celb_part_gap

 
i=-1
j=2
k=1
coor_focals[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j)
coor_protei[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j) + ((rad*(2./3.)**0.5*math.cos(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.cos(theta*j))*2.*celb_part_gap 
coor_focals[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j)
coor_protei[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j) + ((rad*(2./3.)**0.5*math.sin(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.sin(theta*j))*2.*celb_part_gap
coor_focals[3*k + 2] = rad*-1.*(1./3.)**0.5
coor_protei[3*k + 2] = rad*-1.*(1./3.)**0.5 + ((rad*-1.*(1./3.)**0.5)**2.0)**0.5/(rad*-1.*(1./3.)**0.5)*2.*celb_part_gap
 
i=-1
j=3
k=2
coor_focals[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j)
coor_protei[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j) + ((rad*(2./3.)**0.5*math.cos(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.cos(theta*j))*2.*celb_part_gap 
coor_focals[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j)
coor_protei[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j) + ((rad*(2./3.)**0.5*math.sin(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.sin(theta*j))*2.*celb_part_gap
coor_focals[3*k + 2] = rad*-1.*(1./3.)**0.5
coor_protei[3*k + 2] = rad*-1.*(1./3.)**0.5 + ((rad*-1.*(1./3.)**0.5)**2.0)**0.5/(rad*-1.*(1./3.)**0.5)*2.*celb_part_gap
 
i=-1
j=4
k=3
coor_focals[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j)
coor_protei[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j) + ((rad*(2./3.)**0.5*math.cos(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.cos(theta*j))*2.*celb_part_gap 
coor_focals[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j)
coor_protei[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j) + ((rad*(2./3.)**0.5*math.sin(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.sin(theta*j))*2.*celb_part_gap
coor_focals[3*k + 2] = rad*-1.*(1./3.)**0.5
coor_protei[3*k + 2] = rad*-1.*(1./3.)**0.5 + ((rad*-1.*(1./3.)**0.5)**2.0)**0.5/(rad*-1.*(1./3.)**0.5)*2.*celb_part_gap
 
i=1
j=1
k=4
coor_focals[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j)
coor_protei[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j) + ((rad*(2./3.)**0.5*math.cos(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.cos(theta*j))*2.*celb_part_gap 
coor_focals[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j)
coor_protei[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j) + ((rad*(2./3.)**0.5*math.sin(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.sin(theta*j))*2.*celb_part_gap
coor_focals[3*k + 2] = rad*(1./3.)**0.5
coor_protei[3*k + 2] = rad*(1./3.)**0.5 + ((rad*(1./3.)**0.5)**2.0)**0.5/(rad*(1./3.)**0.5)*2.*celb_part_gap

i=1
j=2
k=5
coor_focals[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j)
coor_protei[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j) + ((rad*(2./3.)**0.5*math.cos(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.cos(theta*j))*2.*celb_part_gap 
coor_focals[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j)
coor_protei[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j) + ((rad*(2./3.)**0.5*math.sin(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.sin(theta*j))*2.*celb_part_gap
coor_focals[3*k + 2] = rad*(1./3.)**0.5
coor_protei[3*k + 2] = rad*(1./3.)**0.5 + ((rad*(1./3.)**0.5)**2.0)**0.5/(rad*(1./3.)**0.5)*2.*celb_part_gap

i=1
j=3
k=6
coor_focals[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j)
coor_protei[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j) + ((rad*(2./3.)**0.5*math.cos(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.cos(theta*j))*2.*celb_part_gap 
coor_focals[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j)
coor_protei[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j) + ((rad*(2./3.)**0.5*math.sin(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.sin(theta*j))*2.*celb_part_gap
coor_focals[3*k + 2] = rad*(1./3.)**0.5
coor_protei[3*k + 2] = rad*(1./3.)**0.5 + ((rad*(1./3.)**0.5)**2.0)**0.5/(rad*(1./3.)**0.5)*2.*celb_part_gap

i=1
j=4
k=7
coor_focals[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j)
coor_protei[3*k + 0] = rad*(2./3.)**0.5*math.cos(theta*j) + ((rad*(2./3.)**0.5*math.cos(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.cos(theta*j))*2.*celb_part_gap 
coor_focals[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j)
coor_protei[3*k + 1] = rad*(2./3.)**0.5*math.sin(theta*j) + ((rad*(2./3.)**0.5*math.sin(theta*j))**2.0)**0.5/(rad*(2./3.)**0.5*math.sin(theta*j))*2.*celb_part_gap
coor_focals[3*k + 2] = rad*(1./3.)**0.5
coor_protei[3*k + 2] = rad*(1./3.)**0.5 + ((rad*(1./3.)**0.5)**2.0)**0.5/(rad*(1./3.)**0.5)*2.*celb_part_gap


 
for i in range(0,nb_focals,1):
   for phase in range(0,2,1):
      body=rigidSphere(r=1.*celb_part_radii,
                    center=numpy.array([coor_focals[3*i+0], coor_focals[3*i+1], coor_focals[3*i+2]]),
                    model=modr, material=cell, color='INTEs')
      nami='I{1:{0}}'.format("04d",i)
      body.addContactors(shape='PT3Dx', color=nami)
      body.addContactors(shape='PT3Dx', color='FOCAL')
      if (phase == 0):
         ce_ls.append(body)
         fa_ls.append(body)
      if (phase == 0):
         charinc = charinc_sprd
#         vspread = 1.
#         hspread = 1.
         bodies1 += body
      elif (phase == 1):
         charinc = charinc_stbl
#         vspread = 0.
#         hspread = 0.
         bodies2 += body
 
for i in range(0,nb_focals,1):
   for phase in range(0,2,1):
      body=rigidSphere(r=1.*celb_part_radii,
                    center=numpy.array([coor_protei[3*i+0], coor_protei[3*i+1], coor_protei [3*i+2]]),
                    model=modr, material=cell, color='PROTs')
      nami='I{1:{0}}'.format("04d",i)
      body.addContactors(shape='PT3Dx', color='PROTp')
      if (phase == 0):
         ce_ls.append(body)
         fa_ls.append(body)
      if (phase == 0):
         charinc = charinc_sprd
         vspread = 1.
         hspread = 1.
         bodies1 += body
      elif (phase == 1):
         charinc = charinc_stbl
         vspread = 0.
         hspread = 0.
         bodies2 += body
      body.imposeDrivenDof(component=1, description='predefined', ct=charinc*coor_protei[3*i + 0], amp=0., omega=0., phi=0., rampi=hspread, ramp=0., dofty='vlocy')
      body.imposeDrivenDof(component=2, description='predefined', ct=charinc*coor_protei[3*i + 1], amp=0., omega=0., phi=0., rampi=hspread, ramp=0., dofty='vlocy')
      body.imposeDrivenDof(component=3, description='predefined', ct=charinc*coor_protei[3*i + 2], amp=0., omega=0., phi=0., rampi=vspread, ramp=0., dofty='vlocy')
 

# Generating the particles of the microfilaments network (randomly oriented and
# homogeneously distributed in the radial direction)
random.seed(seed_mf)
for i in range(0,nb_nodes_mf,1):
  while 1:
     rad = random.uniform(diam_core/2.+ 2.*corb_part_radii+radii_part_networks, diam_cell/2-radii_part_networks)
     theta = random.uniform(0., 2.*math.pi)
     phi = random.uniform(-math.pi/2., math.pi/2.)
     px = xctr_cell+(rad)*(math.cos(phi)*math.cos(theta))
     py = yctr_cell+(rad)*(math.cos(phi)*math.sin(theta))
     pz = zctr_cell+(rad)*(math.sin(phi))
     dist_cell_cntr = ((px-xctr_cell)**2.0+(py-yctr_cell)**2.0+(pz-zctr_cell)**2.0)**0.5
     if((dist_cell_cntr - radii_part_networks > diam_core/2. + core_part_radii) and
        (dist_cell_cntr + radii_part_networks < diam_cell/2.)):
        break
  center_mf_node=numpy.array([px, py, pz])

  body=rigidSphere(r=radii_part_networks, center=center_mf_node,
        model=modr, material=cell, color='MFsxx')
  body.addContactors(shape='PT3Dx', color='MFpxx')
  bodies1 += body
  bodies2 += body
  ce_ls.append(body)

# Generating the particles of the microtubules network (oriented in prefered
# direction and homogeneously distributed in the radial direction)
random.seed(seed_mt)
dir_phi = random.uniform(0., 2.*math.pi)
dir_theta = random.uniform(-math.pi/2., math.pi/2.)

for i in range(0,nb_nodes_mt,1):
  while 1:
     rad = random.uniform(diam_core/2.+ 2.*corb_part_radii+radii_part_networks, diam_cell/2.-radii_part_networks)
     while 1:
         theta = random.gauss(dir_theta, math.pi/2.)
         if (theta>dir_theta-math.pi and theta<dir_theta+math.pi):
            break
     while 1:
         phi = random.gauss(dir_phi, math.pi/2.)
         if (phi>dir_phi-math.pi/2. and phi<dir_phi+math.pi/2.):
            break
     px = xctr_cell+(rad)*(math.cos(phi)*math.cos(theta))
     py = yctr_cell+(rad)*(math.cos(phi)*math.sin(theta))
     pz = zctr_cell+(rad)*(math.sin(phi))
     dist_cell_cntr = ((px-xctr_cell)**2.0+(py-yctr_cell)**2.0+(pz-zctr_cell)**2.0)**0.5
     if((dist_cell_cntr - radii_part_networks > diam_core/2. + core_part_radii) and
        (dist_cell_cntr + radii_part_networks < diam_cell/2.)):
        break
  center_mt_node=numpy.array([px, py, pz])

  body=rigidSphere(r=radii_part_networks, center=center_mt_node,
        model=modr, material=cell, color='MTsxx')
  body.addContactors(shape='PT3Dx', color='MTpxx')

  bodies1 += body
  bodies2 += body
  ce_ls.append(body)

# Generating the particles of the intermediate filaments network (randomly oriented and
# normally distributed in the radial direction centered on the nucleus diameter, only keeping particle outside the nucleus)
random.seed(seed_if)

for i in range(0,nb_nodes_if,1):
  while 1:
     while 1:
        rad = random.gauss(diam_core/2.+2*corb_part_radii, diam_cell/2.)
        if (rad>diam_core/2 + 2.*corb_part_radii + radii_part_networks and rad<max(diam_cell,diam_cm2)/2-radii_part_networks):
            break
     theta = random.uniform(0., 2.*math.pi)
     phi = random.uniform(-math.pi/2., math.pi/2.)
     px = xctr_cell+(rad)*(math.cos(phi)*math.cos(theta))
     py = yctr_cell+(rad)*(math.cos(phi)*math.sin(theta))
     pz = zctr_cell+(rad)*(math.sin(phi))
     dist_cell_cntr = ((px-xctr_cell)**2.0+(py-yctr_cell)**2.0+(pz-zctr_cell)**2.0)**0.5
     dist_cell_cntr = ((px-xctr_cell)**2.0+(py-yctr_cell)**2.0+(pz-zctr_cm1)**2.0)**0.5
     dist_cm2_cntr = ((px-xctr_cell)**2.0+(py-yctr_cell)**2.0+(pz-zctr_cm2)**2.0)**0.5
     if((dist_cell_cntr - radii_part_networks > diam_core/2. + core_part_radii) and
        (dist_cell_cntr + radii_part_networks < diam_cell/2.) and (dist_cm2_cntr + radii_part_networks < diam_cm2/2.)):
        break
  center_if_node=numpy.array([px, py, pz])

  body=rigidSphere(r=radii_part_networks, center=center_if_node,
        model=modr, material=cell, color='IFsxx')
  body.addContactors(shape='PT3Dx', color='IFpxx')

  bodies1 += body
  bodies2 += body
  ce_ls.append(body)

# Generating the stress fibres as random chains of particles
seen_focals = numpy.zeros((nb_focals+1), dtype='i')
numfibre_focals = numpy.zeros((nb_focals+1), dtype='i')
nfibre = 0
for i in range(0,nb_focals,1):
  npart_fibre = 0
  if(seen_focals[i] == 0):
     numfibre_focals[i] = nfibre
     seen_focals[i] = 1

     ax = coor_focals[3*i+0]
     ay = coor_focals[3*i+1]
     az = coor_focals[3*i+2]

     dist_cell_cntr = ((ax-xctr_cell)**2.0+(ay-yctr_cell)**2.0+(az-zctr_cell)**2.0)**0.5
     phi_ref = math.asin((zctr_cell-az)/(dist_cell_cntr))

     cos_theta_ref = (xctr_cell-ax)/(dist_cell_cntr*math.cos(phi_ref))
     sin_theta_ref = (yctr_cell-ay)/(dist_cell_cntr*math.cos(phi_ref))
     if (sin_theta_ref > 0.):
        theta_ref = math.acos(cos_theta_ref)
     else:
        theta_ref = 2.*math.pi - math.acos(cos_theta_ref)

     while 1:
        ox = ax
        oy = ay
        oz = az

        ntry = 0
        while (ntry < 10000):
           rad = random.uniform(sf_inc_min, sf_inc_max)
           theta = random.uniform(-1*shoot_cone_angle, shoot_cone_angle)
           phi = random.uniform(-1*shoot_cone_angle, shoot_cone_angle)

           ax = ox + (rad)*(math.cos(phi+phi_ref)*(math.cos(theta+theta_ref)))
           ay = oy + (rad)*(math.cos(phi+phi_ref)*(math.sin(theta+theta_ref)))
           az = oz + (rad)*(math.sin(phi+phi_ref))

           dist_cell_cntr = ((ax-xctr_cell)**2.0+(ay-yctr_cell)**2.0+(az-zctr_cell)**2.0)**0.5
           dist_cell_cntr = ((ax-xctr_cell)**2.0+(ay-yctr_cell)**2.0+(az-zctr_cm1)**2.0)**0.5
           dist_cm2_cntr = ((ax-xctr_cell)**2.0+(ay-yctr_cell)**2.0+(az-zctr_cm2)**2.0)**0.5
           if((dist_cell_cntr - fibers_part_radii > diam_core/2. + core_part_radii) and
              (dist_cell_cntr + fibers_part_radii < diam_cell/2.) and (dist_cm2_cntr + fibers_part_radii < diam_cm2/2.)):
              break
           if (ntry == 5001):
               theta_ref = -1.*theta_ref
               phi_ref = -1.*phi_ref

        center_sf_node = numpy.array([ax, ay, az])
        name='Ss{1:{0}}'.format("03d",nfibre)
        body=rigidSphere(r=1.*fibers_part_radii, center=center_sf_node,
           model=modr, material=cell, color=name)

        dist_min_unseen_focal = diam_cell/2.
        closest_focal = i
        for j in range(0,nb_focals,1):
           if(seen_focals[j] == 0):
              dist_unseen_focal = ((ax-coor_focals[2*j])**2.0+(ay-coor_focals[2*j+1])**2.0)**0.5
              if(dist_unseen_focal < dist_min_unseen_focal):
                 dist_min_unseen_focal = dist_unseen_focal
                 closest_focal = j
        dist_cell_cntr = ((ax-xctr_cell)**2.0+(ay-yctr_cell)**2.0+(az-zctr_cm1)**2.0)**0.5
        dist_cm2_cntr = ((ax-xctr_cell)**2.0+(ay-yctr_cell)**2.0+(az-zctr_cm2)**2.0)**0.5
        dist_env1 = diam_cell/2. - dist_cell_cntr - fibers_part_radii
        dist_env2 = diam_cm2/2. - dist_cm2_cntr - fibers_part_radii

        if(ntry < 10000):
           if((dist_env1 < celb_part_radii or dist_env2 < celb_part_radii) and npart_fibre > 20): ## pour eviter la fin trop rapide de la generation
              name='Se{1:{0}}'.format("03d",nfibre) ## point enveloppe dorsale
              body.addContactors(shape='PT3Dx', color=name)
              break
           else:
              if (npart_fibre == 0):
                 name='Sf{1:{0}}'.format("03d",nfibre) ## point sur point focal
                 body.addContactors(shape='PT3Dx', color=name)
              else:
                 name='Sp{1:{0}}'.format("03d",nfibre) ## point interne
                 body.addContactors(shape='PT3Dx', color=name)

              bodies1 += body
              bodies2 += body
              ce_ls.append(body)
           npart_fibre += 1
        else:
           sys.exit("Error message ntry exceed 10000")

     bodies1 += body
     bodies2 += body
     ce_ls.append(body)
     nfibre +=1
nb_fibres = nfibre


# Generating particles interactions in both phases
for i in range(0,2,1):
   tacts_tmp = tact_behavs()
## The spreading phase is performed with extremely stiff interactions to avoid
## vibrations and divergence of the solver
   if (i == 0):
      visco_pen = 0.00000001
      stiff_pen = 5000.00
      pstr_pen = 0.00
## The stabilization phase is performed with calibrated and accurate parameters
## to produce realistic forces in the cell
   elif (i == 1):
      visco_pen = 1.0
      stiff_pen = 0.20
      pstr_pen = 2.00

   ctnrm = tact_behav(name='ictn0',law='IQS_CLB',fric=0.07*visco_pen)
   tacts_tmp+=ctnrm

   sprmt=tact_behav(name='ismt0',law='ELASTIC_ROD' , stiffness= 22.8e+00*stiff_pen, prestrain=+0.00)
   tacts_tmp+=sprmt

   cabif=tact_behav(name='icif0',law='ELASTIC_WIRE', stiffness= 15.7e+00*stiff_pen, prestrain=+0.10)
   tacts_tmp+=cabif

   cabmf=tact_behav(name='icmf0',law='ELASTIC_WIRE', stiffness= 1.00e+00*stiff_pen, prestrain=-0.02*pstr_pen)
   tacts_tmp+=cabmf

   cabsf=tact_behav(name='icsf0',law='ELASTIC_WIRE', stiffness= 9.14e+00*stiff_pen, prestrain=-0.20*pstr_pen)
   tacts_tmp+=cabsf

   cabnm=tact_behav(name='icnm0',law='ELASTIC_ROD', stiffness= 1.00e+00*stiff_pen, prestrain=+0.00)
   tacts_tmp+=cabnm

   cabcm=tact_behav(name='iccm0',law='ELASTIC_ROD', stiffness= 5.00e+00*stiff_pen, prestrain=-0.02*pstr_pen)
   tacts_tmp+=cabcm

   cabfa=tact_behav(name='icfa0',law='ELASTIC_WIRE', stiffness= 10.0e+00*stiff_pen, prestrain=+0.00)
   tacts_tmp+=cabfa

   cabci=tact_behav(name='icci0',law='ELASTIC_WIRE', stiffness= 100.e+00*stiff_pen, prestrain=+0.00)
   tacts_tmp+=cabci

   caban=tact_behav(name='adhes',law='ELASTIC_WIRE', stiffness= 1.0e+08*stiff_pen, prestrain=-0.00*pstr_pen)
   tacts_tmp+=caban

   if (i == 0):
      tacts_stif = tacts_tmp
   elif (i == 1):
      tacts_nrml = tacts_tmp

# Creating the visibility tables that determine pairs of interacting particle interact
## for contact interaction
vinin = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='INOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='INOYs',alert=core_part_radii)
svs+=vinin

vinbn = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='INOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='BNOYs',alert=max(core_part_radii,corb_part_radii))
svs+=vinbn

vincs = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='INOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='CTSLs',alert=max(core_part_radii,ctsl_part_radii))
svs+=vincs

vbnsu = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BNOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='SUBST',alert=0.1*corb_part_radii)
svs+=vbnsu

vbtbt = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BCELt',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='BCELt',alert=celb_part_radii/10.)
svs+=vbtbt

vbnbc = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BNOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='BCELs',alert=max(corb_part_radii,celb_part_radii))
svs+=vbnbc

vbnmf = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BNOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='MFsxx',alert=corb_part_radii)
svs+=vbnmf

vbnmt = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BNOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='MTsxx',alert=corb_part_radii)
svs+=vbnmt

vbnif = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BNOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='IFsxx',alert=corb_part_radii)
svs+=vbnif

vbcsu = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BCELs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='SUBST',alert=0.1*celb_part_radii)
svs+=vbcsu

vbcmf = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BCELs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='MFsxx',alert=celb_part_radii)
svs+=vbcmf
##
vbcif = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BCELs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='IFsxx',alert=celb_part_radii)
svs+=vbcif

vbcmt = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BCELs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='MTsxx',alert=celb_part_radii)
svs+=vbcmt

vbncs = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BNOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='CTSLs',alert=ctsl_part_radii)
svs+=vbncs

vbccs = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='BCELs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='CTSLs',alert=ctsl_part_radii)
svs+=vbccs

vcscs = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='CTSLs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='CTSLs',alert=ctsl_part_radii)
svs+=vcscs

vcssu = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
   colorCandidat='CTSLs',behav=ctnrm, CorpsAntagoniste='RBDY3',
   antagoniste='SPHER',colorAntagoniste='SUBST',alert=ctsl_part_radii)
svs+=vcssu

for i in range(0,nb_fibres,1):
   name='Ss{1:{0}}'.format("03d",i)
   vbnsf = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
      colorCandidat='BNOYs',behav=ctnrm, CorpsAntagoniste='RBDY3',
      antagoniste='SPHER',colorAntagoniste=name,alert=max(corb_part_radii,fibers_part_radii))
   svs+=vbnsf

for i in range(0,nb_fibres,1):
   name='Ss{1:{0}}'.format("03d",i)
   vsfbc = see_table(CorpsCandidat='RBDY3',candidat='SPHER',
      colorCandidat=name,behav=ctnrm, CorpsAntagoniste='RBDY3',
      antagoniste='SPHER',colorAntagoniste='BCELs',alert=max(corb_part_radii,fibers_part_radii))
   svs+=vsfbc

## for elastic_rod interaction (microtubules)
vmtmt = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='MTpxx',behav=sprmt, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='MTpxx',alert=alert_intra_network)
svs+=vmtmt

## for elastic_wire interaction (intermediate filaments)
vifif = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='IFpxx',behav=cabif, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='IFpxx',alert=alert_intra_network)
svs+=vifif

## for elastic_wire interaction (microfilaments)
vmfmf = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='MFpxx',behav=cabmf, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='MFpxx',alert=alert_intra_network)
svs+=vmfmf

## for elastic_wire interaction (stress fibres)
for i in range(0,nb_fibres,1):
   nam1='Sp{1:{0}}'.format("03d",i)
   vspsp = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
      colorCandidat=nam1,behav=cabsf, CorpsAntagoniste='RBDY3',
      antagoniste='PT3Dx',colorAntagoniste=nam1,alert=2.*sf_inc_max)
   svs+=vspsp

   nam2='Sf{1:{0}}'.format("03d",i)
   vspsf = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
      colorCandidat=nam1,behav=cabsf, CorpsAntagoniste='RBDY3',
      antagoniste='PT3Dx',colorAntagoniste=nam2,alert=2.*sf_inc_max)
   svs+=vspsf

   nam2='Se{1:{0}}'.format("03d",i)
   vspse = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
      colorCandidat=nam1,behav=cabsf, CorpsAntagoniste='RBDY3',
      antagoniste='PT3Dx',colorAntagoniste=nam2,alert=2.*sf_inc_max)
   svs+=vspse

## for elastic_wire interaction (nucleus membrane)
vbnbn = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='BNOYp',behav=cabnm, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='BNOYp',alert=2.*corb_part_gap)
svs+=vbnbn

## for elastic_wire interaction (cell membrane)
vbcbc = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='BCELp',behav=cabcm, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='BCELp',alert=2.*celb_part_gap)
svs+=vbcbc

## for elastic_wire interaction (focal adhesions)
for i in range(0,nb_focals,1):
   name='I{1:{0}}'.format("04d",i)
   vbcfo = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
      colorCandidat='BCELp',behav=cabfa, CorpsAntagoniste='RBDY3',
      antagoniste='PT3Dx',colorAntagoniste=name,alert=5.*celb_part_gap)
   svs+=vbcfo

   nam2='Sf{1:{0}}'.format("03d",numfibre_focals[i])
   vfosf = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
      colorCandidat=name,behav=cabsf, CorpsAntagoniste='RBDY3',
      antagoniste='PT3Dx',colorAntagoniste=nam2,alert=200.*sf_inc_max)
   svs+=vfosf

vancr = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='FOCAL',behav=caban, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='PROTp',alert=8.*celb_part_gap)
svs+=vancr


## for elastic_wire interaction (networks interactions)
vmtmf = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='MTpxx',behav=cabci, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='MFpxx',alert=alert_extra_network)
svs+=vmtmf

for i in range(0,nb_fibres,1):
   name='Sp{1:{0}}'.format("03d",i)
   vifsf = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
      colorCandidat='IFpxx',behav=cabci, CorpsAntagoniste='RBDY3',
      antagoniste='PT3Dx',colorAntagoniste=name,alert=alert_extra_network)
   svs+=vifsf

   vmfsf = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
      colorCandidat='MFpxx',behav=cabci, CorpsAntagoniste='RBDY3',
      antagoniste='PT3Dx',colorAntagoniste=name,alert=alert_extra_network)
   svs+=vmfsf

vmfbn = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='MFpxx',behav=cabci, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='BNOYp',alert=alert_extra_network)
svs+=vmfbn

vmfbc = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='MFpxx',behav=cabci, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='BCELp',alert=alert_extra_network)
svs+=vmfbc

vifbn = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
   colorCandidat='IFpxx',behav=cabci, CorpsAntagoniste='RBDY3',
   antagoniste='PT3Dx',colorAntagoniste='BNOYp',alert=alert_extra_network)
svs+=vifbn

for i in range(0,nb_fibres,1):
   name='Se{1:{0}}'.format("03d",i)
   vbcsf = see_table(CorpsCandidat='RBDY3',candidat='PT3Dx',
      colorCandidat='BCELp',behav=cabci, CorpsAntagoniste='RBDY3',
      antagoniste='PT3Dx',colorAntagoniste=name,alert=alert_extra_network)
   svs+=vbcsf

# Writing the configuration files for the spreading phase
tactsfinal = tacts_stif
svsfinal = svs
writeBodies(bodies1,chemin='DATBOX_SPRD/')
writeBulkBehav(mat,chemin='DATBOX_SPRD/')
writeTactBehav(tactsfinal,svsfinal,chemin='DATBOX_SPRD/')
writeDrvDof(bodies1,chemin='DATBOX_SPRD/')
writeDofIni(bodies1,chemin='DATBOX_SPRD/')
writeVlocRlocIni(chemin='DATBOX_SPRD/')

post = postpro_commands()
nlgs = postpro_command(name='SOLVER INFORMATIONS', step=1)
post.addCommand(nlgs)
ener = postpro_command(name='KINETIC ENERGY', step=max(1,int(nsteps_sprd/noutp)))
post.addCommand(ener)
rset = postpro_command(name='NEW RIGID SETS', step=max(1,int(nsteps_sprd/noutp)), rigid_sets=[ns_ls, ce_ls])
post.addCommand(rset)
torq = postpro_command(name='TORQUE EVOLUTION', step=max(1,int(nsteps_sprd/noutp)), rigid_set=fa_ls)
post.addCommand(torq)
trac = postpro_command(name='BODY TRACKING', step=max(1,int(nsteps_sprd/noutp)), rigid_set=ns_ls)
post.addCommand(trac)
writePostpro(commands=post, parts=bodies1, path='DATBOX_SPRD/')

# Writing the configuration files for the stabilization phase
tactsfinal = tacts_nrml
svsfinal = svs
writeBodies(bodies2,chemin='DATBOX_STBL/')
writeBulkBehav(mat,chemin='DATBOX_STBL/')
writeTactBehav(tactsfinal,svsfinal,chemin='DATBOX_STBL/')
writeDrvDof(bodies2,chemin='DATBOX_STBL/')
writeDofIni(bodies2,chemin='DATBOX_STBL/')
writeVlocRlocIni(chemin='DATBOX_STBL/')

post = postpro_commands()
nlgs = postpro_command(name='SOLVER INFORMATIONS', step=1)
post.addCommand(nlgs)
ener = postpro_command(name='KINETIC ENERGY', step=max(1,int(nsteps_stbl/noutp)))
post.addCommand(ener)
rset = postpro_command(name='NEW RIGID SETS', step=max(1,int(nsteps_stbl/noutp)), rigid_sets=[ns_ls, ce_ls])
post.addCommand(rset)
torq = postpro_command(name='TORQUE EVOLUTION', step=max(1,int(nsteps_stbl/noutp)), rigid_set=fa_ls)
post.addCommand(torq)
trac = postpro_command(name='BODY TRACKING', step=max(1,int(nsteps_stbl/noutp)), rigid_set=ce_ls)
post.addCommand(trac)
writePostpro(commands=post, parts=bodies2, path='DATBOX_STBL/')
