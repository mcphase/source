cd ../examples/tmnic2

cif2mcphas -rl edvardson -pc 10 -sp -nm Ni -ch Tm=2.15939417656068,Ni=-0.338086596354515,C=3.22959620989701,E1=-2.6945,E2=-1.5506,E3=-0.3955,E4=-0.2008,E5=-0.2084,E6=-0.2369 TmNiC2_400K_mod.cif

# <h2> Calculate Crystal Field Phonon Coupling and Physical Properties </h2> 

cp  Tm1_1.sipf Tm_new.sipf

setvariable MAGNETIC 0 Ni*.sipf

# from cif2mcphas the extended crystal unit cell is created
#  ... the primitive
# unit cell should be used to minimize calculation time
#! r1a= 1   r2a= 0 r3a=   0
#! r1b= 0 r2b=   1   r3b=   0.5   primitive lattice vectors [a][b][c]
#! r1c= 0 r2c=   0 r3c=   0.5 

setvariable r3b 0.5 mcphas_all_atoms.j
setvariable r3c 0.5 mcphas_all_atoms.j

reduce_unitcell mcphas_all_atoms.j > mcphas.j
#remove all sipf files which are not needed any more
perl -l -n -e "unlink" reduce_unitcell_sipf.del


# <h3> cfph interaction with rmax large enough to reach convergence </h3>

# second the cfph interaction is calculated from a pointcharge model 
# with rmax=20 Angstroem
# crystal field parameters are calculated too, in this step and 
# stored in file results/makenn.a18.sipf (which models the Tm 4f electron shell)
# use option -r to put phonon module in Tm*.sipf
 makenn 10 -cfph -r


# 0) distribute C-C bonding electron equally on both C, put Ni-C bond electron on Ni, Tm-C on Tm, Tm-Ni on Tm
# 
 reduce_unitcell -delatoms 5:3+-0.5:4+-0.5,6:2+-1,7:2+-1,8:2+-1,9:2+-1,10:1+-1,11:1+-1,12:1+-1,13:1+-1,14:1+-1,15:1+-1,16:1+-1,17:1+-1 results/makenn.j > mnn.j

getvalue.pl -c -2.376716e-01 0 7 1 0 mnn.j


cp mcphas.j mcphas_all_atoms.j
# put to mcphas.j also a file without the E charges - so phonons can be modelled now ...
reduce_unitcell -delatoms 5,6,7,8,9,10,11,12,13,14,15,16,17 mcphas_all_atoms.j > mcphas.j 

# <h3> Set-up of the Phonon Interactions using a Born v Karman Model</h3>

# create table of bvk springs and modify it to reproduce high energy optical phonons in RNiC2
 makenn 5.1 -bvk > bvk_springs.dat
 
 fillcol 4 '1000xexp(-0.5xc3xc3)' bvk_springs.dat

# compute makenn.j file
 makenn 5.1 -bvk bvk_springs.dat

cp results/makenn.j mcphasph.j


getvalue.pl -c +2.698873023e+02 0 4 1 0 mcphasph.j

# <h3> Sum Phonon and CF Phonon interaction parameters and calculate physical properties </h3>
#  sum of Born v Karman model parameters (as listed in mcphasph.j)
# and crystal field phonon interaction (gener ated by makenn) by program addj


# set maximum displacement for all phonons to umax=0.5 a0 with a0=0.52 Angsgtroem
 setvariable MODPAR8 0.5 *.sipf


# second with cf-Phonon interaction parameters 
 addj  mcphasph.j mnn.j > mcphas.j
 mcphasit -v -doeps  -prefix test_ 
getvalue.pl -c 52.85719 0 9 1 0 results/test_mcphas.fum

# <h2> Magnetic Scattering </h2>

 setup_mcdisp_mf -prefix test_ 12 0 0 1
 mcdispit -prefix test_ -pinit 0.1

range 7 0.2 0.4 results/test_mcdisp.qei
getvalue.pl -c 0.005526 9 10 13.9967 0 results/test_mcdisp.qei
getvalue.pl -c 0.005499 9 11 13.9967 0 results/test_mcdisp.qei
getvalue.pl -c 0.000816 9 12 13.9967 0 results/test_mcdisp.qei
  

cd ../../demo
