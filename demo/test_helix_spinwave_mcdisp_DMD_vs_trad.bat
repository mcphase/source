cd ../examples/helix_spinwave

REM calculate the magnetic structure ***************
call mcphasit

getvalue.pl -c 0.64445  1 8 6 0 results/mcphas.fum

REM calculate the excitations at T=2 H=(0 0 0) *****
call setup_mcdisp_mf 2 0 0 0

REM -max 6 ... at least 6 transitions have to be included to get the 
REM soft mode correctly 
call mcdispit -max 6

range 7 0.43 0.45 results/mcdisp.qei
getvalue.pl -c 0.04562 9 10 2.47221  0 results/mcdisp.qei

call mcdispit -r 0.5 -max 6
range -d 7 0.444 0.445 results/mcdisp.dsigma

getvalue.pl -c 0.535 8 9 1.751 0 results/mcdisp.dsigma
getvalue.pl -c 5.98  8 10 1.751 0 results/mcdisp.dsigma
getvalue.pl -c 11.9   8 11 1.751 0 results/mcdisp.dsigma
getvalue.pl -c -7.95 8 12 1.751 0 results/mcdisp.dsigma
getvalue.pl -c -3.17 8 13 1.751 0 results/mcdisp.dsigma


cd ../../demo
