cd ../examples/tmcu2_cf_phonon

call cif2mcphas -nm Cu  -pc 4.4 -ch Tm=0.8,Cu=-0.4 TmCu2.cif


call setvariable r1a 0.5 mcphas_all_atoms.j
call setvariable r1b 0.5 mcphas_all_atoms.j
call setvariable r1c 0.5 mcphas_all_atoms.j

call reduce_unitcell mcphas_all_atoms.j > mcphas.j

REM remove all sipf files which are not needed any more
call perl -l -n -e "unlink" reduce_unitcell_sipf.del


REM <h2> Calculate Crystal Field Phonon Coupling </h2> 

REM second the cfph interaction is calculated from a pointcharge model 
REM with rmax=15 Angstroem
REM use option -r to put phonon module into Tm*.sipf

call makenn 15 -cfph -r
call copy results/makenn.j results/mnn.j

REM <h2> PHONON CALCULATION </h2>

echo "create mcphas.j with interactions and sipf file K matrices from "
echo "Born van Karman springs in bvk_springs.dat"

call makenn 5.1 -bvk > bvk_springs.dat
REM leave transversal springs zero and setup phonon mcphas.j
call makenn 5.1 -bvk bvk_springs.dat

copy results/makenn.j mcphas.j

copy mcphas.j mcphasph.j

REM get a zero mcdisp.mf file to calculate phonons
copy mcdispph.mf mcdisp.mf

call mcdispit -prefix 001

range 7 0.49 0.51 results/001mcdisp.qei
getvalue.pl -c 0.000234  9 12 11.5574       0    results/001mcdisp.qei

REM <h2> DMD Method for CF PHONON INTERACTION PROBLEM </h2>


REM generate the sum of Born v Karman model parameters (as listed in mcphasph.j)
REM and crystal field phonon interaction (gener ated by makenn) by program addj
REM and start mcdispit to calculate dispersion - hkl specified in mcdisp.par
REM Temperature T=2K specified in mcdisp.mf

REM make maximum displacement for all phonons to umax=1.0 a0 with a0=0.52 Angsgtroem
call setvariable MODPAR8 0.1 *.sipf

call setvariable xmax 20 mcphas.ini


REM second with cf-Phonon interaction parameters
call addj mcphasph.j results/mnn.j > mcphas.j
call mcphasit -v -doeps
call setup_mcdisp_mf 20  0 0 0    
call mcdispit -prefix test_ -pinit 0.1

range 7 0.49 0.51 results/test_mcdisp.qei
getvalue.pl -c 0.000529  9 12 12.9345     0    results/test_mcdisp.qei
getvalue.pl -c 0.01138  9 11 12.9345      0    results/test_mcdisp.qei
getvalue.pl -c 0.01141   9 10 12.9345      0    results/test_mcdisp.qei
 
REM  <h2> diffraction pattern </h2>

REM remove autocreated mcdiff.in with nonmagnetic atoms (otherwise these will be counted twice)
call del mcdiff.in

call setup_mcdiff_in 20 0 0 0
call setvariable lambda 2.45 mcdiff.in
call setvariable lorentz 2 mcdiff.in
call setvariable thetamax 30 mcdiff.in
call mcdiff
average -av -dmin=0.001 6 results/mcdiff.out
getvalue.pl -c 0.01028  6 7 39.715     0    results/mcdiff.out

cd ../../demo
