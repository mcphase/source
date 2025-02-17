#<!--mcphase.mcphas.j-->
#***************************************************************
# Lattice and Exchange Parameter file for
# mcphas version 3.0
# - program to calculate static magnetic properties
# reference: M. Rotter JMMM 272-276 (2004) 481
# mcdisp version 3.0
# - program to calculate the dispersion of magnetic excitations
# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751
#***************************************************************
#
# W1 (2) 
#
# Lattice Constants (A)
#
#! a= 3.164800 b= 3.164800 c= 3.164800  alpha= 90.000000 beta= 90.000000 gamma= 90.000000
#
#! r1a=   1 r2a= 0 r3a=  0
#! r1b=   0 r2b= 1 r3b=  0   primitive lattice vectors [a][b][c]
#! r1c=   0 r2c= 0 r3c=  1
#
#! nofatoms= 2  nofcomponents=3  number of atoms in primitive unit cell/number of components of each spin
#*************************************************************************
#ATOM TYPE W1 ; number of the atom in the UNIT CELL = 1 ; number of the atom within this type = 2
#! da= 0.0000 [a] db= 0.0000 [b] dc= 0.0000 [c] nofneighbours=0 diagonalexchange=1  sipffilename= W2p1.sipf
#-------------------------------------------------------------------------------------
#da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]
#*************************************************************************
#ATOM TYPE W1 ; number of the atom in the UNIT CELL = 2 ; number of the atom within this type = 2
#! da= 0.5000 [a] db= 0.5000 [b] dc= 0.5000 [c] nofneighbours=0 diagonalexchange=1 sipffilename= W2p2.sipf
#-------------------------------------------------------------------------------------
#da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]
