Below scripts are written in Fortran 95. While code is used to compute the trajectory of a photon
in a curved spacetime around Black Holes.

Files main or main_more_flags contain instruction to compile whole code.

File main.f95 is the core of the whole code.
In this file we may set few grids of parameters (the computation will calculate many trajectories)
or set parameters for a single trajectory.
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
wspl.f95, 

contain plenty of functions, subroutines, declarations of variables.

File "zmienne_i_parametry" contains the whole description of variables (in polish for now).

File "moduly" contains the descriptions of particular files (in polish for now). The description
for last 3 files will appear soon.
