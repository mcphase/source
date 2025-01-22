cd ../examples/Cu2OSeO3/4interacting_tetrahedrons_noDM


rem the positions in mcphas_magnetc_atoms.j have been created by powdercell2j
rem ... mind: the output of powdercell positional coordinates had to be shifted
rem by 0.25 0.25 0.25 in order to bring all tetrahedrons into the unit cell.
copy mcphas_magnetic_atoms.j mcphas.j
call makenn 6.4 -rkky 1 1 
call fform 4 6 10.8f results/makenn.j
rem JFM_W
call substitute 0.00436137 4.31 results/makenn.j
rem JAFM_S
call substitute 0.00431320 -14.654 results/makenn.j
rem JFM_S
call substitute 0.00369989 11.0336 results/makenn.j
call substitute 0.00369988 11.0336 results/makenn.j
rem JAF_W
call substitute 0.00330466 -2.3274  results/makenn.j
rem JAF_O..O
call substitute  0.00048307 -3.879   results/makenn.j

rem zero other couplings
call substitute -0.00018414 0   results/makenn.j
call substitute -0.00005257 0   results/makenn.j
call substitute  0.00051291 0   results/makenn.j
call substitute  0.00050946 0   results/makenn.j
call substitute  -0.00002940 0   results/makenn.j
call substitute  -0.00022272 0   results/makenn.j
call substitute  -0.00010512 0   results/makenn.j
call substitute  0.00000811 0   results/makenn.j
call substitute  0.00002439 0   results/makenn.j
call substitute   0.00006034 0   results/makenn.j
call substitute   0.00051850 0   results/makenn.j
call substitute  0.00050987 0   results/makenn.j
call substitute  0.00048384 0   results/makenn.j
call substitute  -0.00022272 0   results/makenn.j


call clusterize results/makenn.j 1 5 9 10 0 2 6 13 16 0 3 7 12 15 0 4 8 11 14

call singleion -r cluster1.sipf 2 0 0 0  0 0 0  0 0 0  0 0 0  0 0 0
call singleion -r cluster2.sipf 2 0 0 0  0 0 0  0 0 0  0 0 0  0 0 0
call singleion -r cluster3.sipf 2 0 0 0  0 0 0  0 0 0  0 0 0  0 0 0
call singleion -r cluster4.sipf 2 0 0 0  0 0 0  0 0 0  0 0 0  0 0 0
getvariable.pl -c -26.592  Eigenvalues results/cluster1.sipf.levels.cef

call mcphasit -v
range 2 4.9 5.1 results/mcphas.fum
getvalue.pl -c -32.90033  1 8 52 0 results/mcphas.fum
getvalue.pl -c 1.699 1 10 52 0 results/mcphas.fum
call setup_mcdisp_mf 2 0 0 0
rem create mcdisp.trs for quasielastic excitations
call mcdispit -c  -maxE 40   -pinit 0.1
getvalue.pl -c 9.35207 0 6 2 0 results/mcdisp.trs
getvalue.pl -c 0.204099  0 8 2 0 results/mcdisp.trs

call mcdispit -t -prefix 011
getvalue.pl -c -39.16   5 9 0.6 0 results/011mcdisp.qom
getvalue.pl -c -35.77   5 13 0.6 0 results/011mcdisp.qom
range 5 0.4 0.55 results/011mcdisp.qei
 average -dmin=0.1 -av  9  results/011mcdisp.qei
getvalue.pl -c 0.2972018 9 10 9.07597 0  results/011mcdisp.qei
getvalue.pl -c 0.23144 9 11 9.07597 0  results/011mcdisp.qei
cd ../../../demo
