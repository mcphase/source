cd ../examples/gdni2b2c

call cif2mcphas -ph -nm Ni1 GdNi2B2C.cif

REM from cif2mcphas the extended crystal unit cell is created     
REM the primitive unit cell should be used to minimize calculation time     
REM r1a= 0.5 r2a=  0.5 r3a=   0
REM r1b=-0.5 r2b=  0.5 r3b=   0   primitive lattice vectors [a][b][c]
REM r1c= 0.5 r2c=  0.5 r3c=   1

call setvariable r1a  0.5 mcphas*.j          
call setvariable r1b -0.5 mcphas*.j          
call setvariable r1c  0.5 mcphas*.j          
call setvariable r2a 0.5 mcphas*.j          
call setvariable r2b 0.5 mcphas*.j          
call setvariable r2c 0.5 mcphas*.j     
call setvariable r3a 0 mcphas*.j          
call setvariable r3b 0 mcphas*.j          
call setvariable r3c 1 mcphas*.j          

   
 reduce_unitcell mcphas_magnetic_atoms.j > mcphas.j  
copy mcphas.j mcphas_magnetic_atoms.j
 reduce_unitcell mcphas_all_atoms.j > mcphas.j  
copy mcphas.j mcphas_all_atoms.j 
perl -l -n -e "unlink" reduce_unitcell_sipf.del

REM <h2> get elastic constants for this material ... </h2>

REM get elastic constants into file makenn.Cel by using makenn
REM do a BvK phonon model to calculate the elastic constants:
call copy mcphas_all_atoms.j mcphas.j
call makenn 7 -bvk > bvk_parameter.tbl
call factcol 4 1.1 bvk_parameter.tbl
call makenn 7 -bvk bvk_parameter.tbl

REM now elastic constants are in results/makenn.Cel


REM <h2> Set up magnetic interactions </h2>

copy mcphas.sipf Gd0.sipf
copy mcphas_magnetic_atoms.j mcphas.j

REM ... several parametrisations for yielding usable J(Q) from J(R)
REM from fit using calcsta.forfit
REM call makenn 10.4 -rkky3d 2100 2.838 2.838 2.192
call makenn 12 -rkky3d -116 1.162 1.162 2.912

copy results/makenn.j mcphas_rkky.j

REM ... remove first neighbour magnetoelastic interactions in order
REM to simulate in accordance with experiment the  dc/c Tdep and Hdep ...
REM and its position derivatives for the first neighbour only
call makenn 3.6 -rkky3d -116 1.162 1.162 2.912 -djdx
call makenn 3.6 -rkky3d -116 1.162 1.162 2.912 -djdy
call makenn 3.6 -rkky3d -116 1.162 1.162 2.912 -djdz

cd results
copy makenn.djdx makenn.1.djdx
copy makenn.djdy makenn.1.djdy
copy makenn.djdz makenn.1.djdz
cd ..

REM and its position derivatives 
call makenn 12 -rkky3d -116 1.162 1.162 2.912 -djdx
call makenn 12 -rkky3d -116 1.162 1.162 2.912 -djdy
call makenn 12 -rkky3d -116 1.162 1.162 2.912 -djdz

call addj -s -1.0 results/makenn.1.djdx results/makenn.djdx > mcphas.djdx
call addj -s -1.0 results/makenn.1.djdy results/makenn.djdy > mcphas.djdy
call addj -s -1.0 results/makenn.1.djdz results/makenn.djdz > mcphas.djdz

REM compute classical dipolar interaction
call makenn 40 
REM add exchang and classical dipolar interaction
addj results/makenn.j mcphas_rkky.j > mcphas.j

call insertfile 21 mcphas.j results/makenn.Cel
REM ... now the calculated elastic constants are finally in mcphas.j

REM check that maximum in J(Q) is really at Q=(0.55 0 0)
call mcdispit -jq -minE -1000

getvariable.pl -c 0.5625 hmax results/mcdisp.jq
getvariable.pl -c 0.3257 jqmax results/mcdisp.jq

REM <h2> compute magnetic properties </h2>
call mcphasit -doeps -prefix test_
getvalue.pl -c -1.93896 0 8 1 0 results/test_mcphas.fum
call setup_mcdiff_in -prefix test_ 2 0 0 0
call setvariable thetamax 20 test_mcdiff.in
call mcdiff -prefix test_
rem select (0.571 0 0)  and (0 0.571 0) 
range -d 6 2.1737E+01 2.1739E+01  results/test_mcdiff.out
sortf 8  results/test_mcdiff.out
getvalue.pl -c 7.4767E-03 0 8 1 0  results/test_mcdiff.out
getvalue.pl -c 7.4767E-03 0 8 2 0  results/test_mcdiff.out
getvalue.pl -c 3.7677E-02 0 8 3 0  results/test_mcdiff.out
getvalue.pl -c 3.7677E-02 0 8 4 0  results/test_mcdiff.out

cd ../../demo

