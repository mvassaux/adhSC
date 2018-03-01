# adhSC

Mechanical model of a single adherent Stem Cell based on the tensegrity theory. The cell, its constituents (nucleus, cytoskeleton, cytosol, and membrane) are described as particles interacting by means of contact, cable or spring. Calibration and checking against experiments has been featured in the following publication:

> Maxime Vassaux and Jean-Louis Milan. [*Stem cell mechanical behaviour modelling: substrate’s curvature influence during adhesion.*](https://doi.org/10.1007/s10237-017-0888-4) Biomechanics and modeling in mechanobiology, 16 (2017): 1295-1308. 

The model is built and simulated using the [LMGC90](https://git-xen.lmgc.univ-montp2.fr/lmgc90/lmgc90_user/wikis/home) software written in Python. Additional scripts are provided to manage the workflow

<img src="https://github.com/mvassaux/mvassaux.github.io/raw/master/static/intra-cell_forces_blutens-redcomp.png" width="300" />

## Execution of the model

Beforehand, LMGC90 must be installed on the machine following these [instructions](https://git-xen.lmgc.univ-montp2.fr/lmgc90/lmgc90_user/wikis/download_and_install).

Generating the input files (`DATBOX/`) is done by the `gen_sample.py` script, and the simulation of the model by the `command.py` script.

All of this is made easier, by gathering common parameters to both scripts and executing them using a single bash script: 

    bash launcher.sh
    
## Post-treatment of the results

Two types of results are produced:
- data stored in plain text format as `.dat` files, additional scripts in `post/` can be found to help extract relevant information
- visualisation in `.vtk` are intended to be opened using software such as [ParaView](https://www.paraview.org/)
