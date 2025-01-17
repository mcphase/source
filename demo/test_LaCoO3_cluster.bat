cd ../examples/LaCoO3_podlesnyak_polaron/2ions

REM cluster with perl in cluster.sipf

REM calculate the clusters energy levels
call singleion -nt 20 -r cluster.sipf 2 0 0 0 0 0 0
getvariable.pl -c -19.6053 Eigenvalues results/cluster.sipf.levels.cef
REM create mcdisp.trs for quasielastic excitations
call mcdispit -c  -maxE 0.9 -minE 0.6  -pinit 0.1
getvalue.pl -c 53.9858 5 7 8 0 results/mcdisp.trs
REM possibly do with powdermagnon to get powder average !!
call powdermagnon 0.1 3 0.1 10 0.7 0.9

REM calculate the quasielastic scattering powder average
call mcdispit -t
call powdermagnon -r results/mcdisp.qei 0.7 0.9 0.3 > results/powdermagnon.clc
getvalue.pl -c 0.262316 5 8 0.2 0 results/powdermagnon.clc

cd ../7ions/

REM with  #!noperl option in cluster.sipf

REM try to do - might not work due to too small memory ... !!!
REM calculate the clusters energy levels
call singleion -nt 20 -r cluster.sipf 2 0 0 0 0 0 0
getvalue.pl -c 0.8633 0 6 4 0 results/cluster.sipf.trs
REM create mcdisp.trs for quasielastic excitations
call mcdispit -c  -maxE 0.9 -minE 0.6  -pinit 0.1
getvalue.pl -c 11.358 5 7 6 0 results/mcdisp.trs
REM possibly do with powdermagnon to get powder average !!
call powdermagnon 0.1 0.2 0.1 3 0.7 0.9
REM calculate the quasielastic scattering powder average
call mcdispit -t

call powdermagnon -r results/mcdisp.qei 0.7 0.9 0.3 

cd ../../../demo

