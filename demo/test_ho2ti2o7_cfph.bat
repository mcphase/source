cd ../examples/ho2ti2o7

cif2mcphas -pc 4 -sp -rl edvardson -ch Ho=1.5568106174469,Ti=0.957270234823227,O1=-0.819195047020912,O2=-0.267656281590462,E1=-0.0899769198149443,E2=-0.186926786787808 Ho2Ti2O7_mod.cif

call setvariable r1a 0 mcphas*.j          
call setvariable r1b 0.5 mcphas*.j          
call setvariable r1c 0.5 mcphas*.j          
call setvariable r2a 0.5 mcphas*.j          
call setvariable r2b 0 mcphas*.j          
call setvariable r2c 0.5 mcphas*.j          
call setvariable r3a 0.5 mcphas*.j          
call setvariable r3b 0.5 mcphas*.j          
call setvariable r3c 0 mcphas*.j     

   
 reduce_unitcell mcphas_magnetic_atoms.j > mcphas.j  
copy mcphas.j mcphas_magnetic_atoms.j
 reduce_unitcell mcphas_all_atoms.j > mcphas.j  
copy mcphas.j mcphas_all_atoms.j 
perl -l -n -e "unlink" reduce_unitcell_sipf.del

REM ***************************************************************
 copy mcphas_all_atoms.j mcphas.j

REM the cfph interaction is calculated from a pointcharge model    
REM  also the crystal field will      
REM be calculated using this model         
call makenn 6 -cfph -r

reduce_unitcell -delatoms 23:9+-0.5:1+-0.5,24:9+-0.5:2+-0.5,25:9+-0.5:3+-0.5,26:9+-0.5:4+-0.5,27:10+-0.5:1+-0.5,28:10+-0.5:2+-0.5,29:10+-0.5:3+-0.5,30:10+-0.5:4+-0.5,31:11+-0.5:1+-0.5,32:11+-0.5:2+-0.5,33:12+-0.5:3+-0.5,34:12+-0.5:4+-0.5,35:13+-0.5:3+-0.5,36:13+-0.5:1+-0.5,37:14+-0.5:2+-0.5,38:14+-0.5:4+-0.5,39:15+-0.5:2+-0.5,40:15+-0.5:3+-0.5,41:16+-0.5:1+-0.5,42:16+-0.5:4+-0.5,43:17+-0.5:1+-0.5,44:17+-0.5:2+-0.5,45:18+-0.5:3+-0.5,46:18+-0.5:4+-0.5,47:19+-0.5:3+-0.5,48:19+-0.5:1+-0.5,49:20+-0.5:2+-0.5,50:20+-0.5:4+-0.5,51:21+-0.5:2+-0.5,52:21+-0.5:3+-0.5,53:22+-0.5:1+-0.5,54:22+-0.5:4+-0.5  results/makenn.j > cfph.j
REM now in cfph.j there is the CF-phonon interaction
perl -l -n -e "unlink" reduce_unitcell_sipf.del

getvariable.pl -c 26 nofatoms cfph.j
getvariable.pl -c 48  nofcomponents cfph.j


REM add all interactions ...
call addj  mcphasph.j cfph.j > mcphas_cfph.j
call addj mcphas_cfph.j mcphasex_cd.j > mcphas.j


call mcphasit -v -doeps -prefix test_ 
getvalue.pl -c -7.1003 1 8 1 0 results/test_mcphas.fum

cd ../../demo
