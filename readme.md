Below scripts are written in Fortran 95. This code is used to compute the trajectory of a photon
in a curved spacetime around Black Holes.

Files main and main_more_flags contain instructions to compile whole code.

File main.f95 is the core of the code.
In this file we may set a few grids of parameters (the code will compute many trajectories 
depending on the grids size) or set parameters for a single trajectory.
In near future I will create a python file to navigate all options and the whole code will be
used only to calculate a single trajectory.

Files:
calki.f95, 
const.f95, 
funkcjecarlsona.f95, 
funkcjejacobiego.f95, 
mukoncowe.f95, 
phitekoncowe.f95, 
pierwiastki.f95, 
pierwiastkimu.f95, 
procgeo.f95, 
przesuniecia.f95, 
redshift.f95, 
ukoncowe.f95, 
wspl.f95 
 contain definitions of many functions, subroutines and variables.

A full explanation of all significant variables is included in file "zmienne_i_parametry" (in Polish for now).

File "moduly" contains the full explanation of functions and/or subroutines in particular files (in Polish for now). 
The description for last 3 files will appear soon.
