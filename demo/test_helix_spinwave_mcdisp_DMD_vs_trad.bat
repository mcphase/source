cd ../examples/helix_spinwave

REM calculate the magnetic structure ***************
call mcphasit

getvalue.pl -c 0.64445  1 8 6 0 results/mcphas.fum

REM calculate the excitations at T=2 H=(0 0 0) *****
call setup_mcdisp_mf 2 0 0 0

REM -max 6 ... at least 6 transitions have to be included to get the 
REM soft mode correctly 
call mcdispit -max 6

range 7 1.22 1.23 results/mcdisp.qei
getvalue.pl -c 0.63 9 10 2.472  0 results/mcdisp.qei

call mcdispit -r 0.5 -max 6
range -d 7 1.22 1.23 results/mcdisp.dsigma

getvalue.pl -c 0.7402 9 10 1.751 0 results/mcdisp.dsigma
getvalue.pl -c 6.07  9 11 1.751 0 results/mcdisp.dsigma
getvalue.pl -c 12.1  9 12 1.751 0 results/mcdisp.dsigma
# chixy sign is determined by turning sense of helix, therefore we omit this check of column 13 and 14
# getvalue.pl -c -7.95 9 13 1.751 0 results/mcdisp.dsigma
# getvalue.pl -c -3.24 9 14 1.751 0 results/mcdisp.dsigma


cd ../../demo
