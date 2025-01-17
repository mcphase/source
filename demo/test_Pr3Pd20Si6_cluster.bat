cd ../examples/Pr3Pd20Si6
call singleion -nt 20 -r cluster.sipf 2 0 0 0 0 0 0
getvariable.pl -c  -6.51653 Eigenvalues results/cluster.sipf.levels.cef

rem calculate the corresponding specific heat
cd results
call cpsingleion 0.1 10 0.01 cluster.sipf.levels.cef > cp.clc

getvalue.pl -c  6.4965507 1 2 0.14 0 cp.clc
cd ../../../demo


