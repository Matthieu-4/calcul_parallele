#Projet de calcul parallèle

## Before anything

make install


## Compile -> bin/main.exe

make


## Cleanning binaries and object files

make clean


## Compile and run

make run

make run H=<overlap> V=[0,1,2]

0 : CG

1 : Schwarz (Dirichlet)

2 : Schwarz (Dirichlet-Neumann)

## Compile and create graphs

make graph


## Create a 3D image of solution

make plot
