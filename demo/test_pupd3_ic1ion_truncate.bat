cd ../examples/pupd3

REM in mcphas 5.5 there was definitely something buggy in
REM output - see emails with Duc autumn 2024 and some bugs were removed
REM (ic1ion and icf1ion were not consistent) 
REM 16.1.2025 tested with respect to XTLS on pu3p.sipf 
REM with alpha=beta=gamma=0  see 

ic1ion pu3p.sipf
getvariable.pl -c -10206 E0 results/ic1ion.out
getvalue.pl -c 0  0 1 1 0  results/ic1ion.out
getvalue.pl -c 0  0 1 2 0  results/ic1ion.out
getvalue.pl -c 18.93  0 1 3 0  results/ic1ion.out
getvalue.pl -c 18.93  0 1 4 0  results/ic1ion.out
getvalue.pl -c 18.93  0 1 5 0  results/ic1ion.out
getvalue.pl -c 18.93  0 1 6 0  results/ic1ion.out
getvalue.pl -c 279.4  0 1 7 0  results/ic1ion.out
getvalue.pl -c 279.4  0 1 8 0  results/ic1ion.out
getvalue.pl -c 462  0 1 9 0  results/ic1ion.out
getvalue.pl -c 462  0 1 10 0  results/ic1ion.out
getvalue.pl -c 583  0 1 11 0  results/ic1ion.out
getvalue.pl -c 583  0 1 12 0  results/ic1ion.out
getvalue.pl -c 583  0 1 13 0  results/ic1ion.out
getvalue.pl -c 583  0 1 14 0  results/ic1ion.out



icf1ion pu3p.sipf
getvariable.pl -c -1195 E0 results/icf1ion.out
getvalue.pl -c 0  0 1 1 0  results/icf1ion.out
getvalue.pl -c 0  0 1 2 0  results/icf1ion.out
getvalue.pl -c 7.36  0 1 3 0  results/icf1ion.out
getvalue.pl -c 7.36  0 1 4 0  results/icf1ion.out
getvalue.pl -c 7.36  0 1 5 0  results/icf1ion.out
getvalue.pl -c 7.36  0 1 6 0  results/icf1ion.out
getvalue.pl -c 93.3  0 1 7 0  results/icf1ion.out
getvalue.pl -c 93.3  0 1 8 0  results/icf1ion.out
getvalue.pl -c 381  0 1 9 0  results/icf1ion.out
getvalue.pl -c 381  0 1 10 0  results/icf1ion.out
getvalue.pl -c 409  0 1 11 0  results/icf1ion.out
getvalue.pl -c 409  0 1 12 0  results/icf1ion.out
getvalue.pl -c 491  0 1 13 0  results/icf1ion.out
getvalue.pl -c 491  0 1 14 0  results/icf1ion.out


mcphasit -prefix test_

getvalue.pl -c -10206.9 2 8 1 0 results/test_mcphas.fum
getvalue.pl -c 0.695 2 10 21 0 results/test_mcphas.fum

spins -c -prefix test_ 1 0 0 81 > dd
getvariable.pl -c 1.410495   "a(0,0)" dd 
getvariable.pl -c 0.000000   "a(2,-2)" dd 
getvariable.pl -c 0.000000   "a(2,-1)" dd 
getvariable.pl -c 0.004368   "a(2,0)" dd 
getvariable.pl -c 0.000000   "a(2,1)" dd 
getvariable.pl -c 0.000000   "a(2,2)" dd 
getvariable.pl -c 0.000000   "a(4,-4)" dd 
getvariable.pl -c 0.000000   "a(4,-3)" dd 
getvariable.pl -c 0.000000   "a(4,-2)" dd 
getvariable.pl -c 0.000000   "a(4,-1)" dd 
getvariable.pl -c 0.035164   "a(4,0)" dd 
getvariable.pl -c 0.000000   "a(4,1)" dd 
getvariable.pl -c 0.000000   "a(4,2)" dd 
getvariable.pl -c 0.000000   "a(4,3)" dd 
getvariable.pl -c 0.029023   "a(4,4)" dd 
getvariable.pl -c 0.000000   "a(6,-6)" dd 
getvariable.pl -c 0.000000   "a(6,-5)" dd 
getvariable.pl -c 0.000000   "a(6,-4)" dd 
getvariable.pl -c 0.000000   "a(6,-3)" dd 
getvariable.pl -c 0.000000   "a(6,-2)" dd 
getvariable.pl -c 0.000000   "a(6,-1)" dd 
getvariable.pl -c -0.074555   "a(6,0)" dd 
getvariable.pl -c 0.000000   "a(6,1)" dd 
getvariable.pl -c 0.000000   "a(6,2)" dd 
getvariable.pl -c 0.000000   "a(6,3)" dd 
getvariable.pl -c 0.202595   "a(6,4)" dd 
getvariable.pl -c 0.000000   "a(6,5)" dd 
getvariable.pl -c 0.000000   "a(6,6)" dd 

spins -s -M -prefix test_ 1 0 0 81 > dd

getvariable.pl -c 0.000000  "aS1(0,0)" dd 
getvariable.pl -c 0.000000  "aS1(1,-1)" dd 
getvariable.pl -c 0.000000  "aS1(1,0)" dd 
getvariable.pl -c 0.000000  "aS1(1,1)" dd 
getvariable.pl -c 0.000000  "aS1(2,-2)" dd 
getvariable.pl -c 0.000000  "aS1(2,-1)" dd 
getvariable.pl -c 0.000000  "aS1(2,0)" dd 
getvariable.pl -c 0.003722  "aS1(2,1)" dd 
getvariable.pl -c 0.000000  "aS1(2,2)" dd 
getvariable.pl -c 0.000000  "aS1(3,-3)" dd 
getvariable.pl -c 0.000000  "aS1(3,-2)" dd 
getvariable.pl -c 0.000000  "aS1(3,-1)" dd 
getvariable.pl -c 0.000000  "aS1(3,0)" dd 
getvariable.pl -c 0.000000  "aS1(3,1)" dd 
getvariable.pl -c 0.000000  "aS1(3,2)" dd 
getvariable.pl -c 0.000000  "aS1(3,3)" dd 
getvariable.pl -c 0.000000  "aS1(4,-4)" dd 
getvariable.pl -c 0.000000  "aS1(4,-3)" dd 
getvariable.pl -c 0.000000  "aS1(4,-2)" dd 
getvariable.pl -c 0.000000  "aS1(4,-1)" dd 
getvariable.pl -c 0.000000  "aS1(4,0)" dd 
getvariable.pl -c 0.070528  "aS1(4,1)" dd 
getvariable.pl -c 0.000000  "aS1(4,2)" dd 
getvariable.pl -c 0.006511  "aS1(4,3)" dd 
getvariable.pl -c 0.000000  "aS1(4,4)" dd 
getvariable.pl -c 0.000000  "aS1(5,-5)" dd 
getvariable.pl -c 0.000000  "aS1(5,-4)" dd 
getvariable.pl -c 0.000000  "aS1(5,-3)" dd 
getvariable.pl -c 0.000000  "aS1(5,-2)" dd 
getvariable.pl -c 0.000000  "aS1(5,-1)" dd 
getvariable.pl -c 0.000000  "aS1(5,0)" dd 
getvariable.pl -c 0.000000  "aS1(5,1)" dd 
getvariable.pl -c 0.000000  "aS1(5,2)" dd 
getvariable.pl -c 0.000000  "aS1(5,3)" dd 
getvariable.pl -c 0.000000  "aS1(5,4)" dd 
getvariable.pl -c 0.000000  "aS1(5,5)" dd 
getvariable.pl -c 0.000000  "aS1(6,-6)" dd 
getvariable.pl -c 0.000000  "aS1(6,-5)" dd 
getvariable.pl -c 0.000000  "aS1(6,-4)" dd 
getvariable.pl -c 0.000000  "aS1(6,-3)" dd 
getvariable.pl -c 0.000000  "aS1(6,-2)" dd 
getvariable.pl -c 0.000000  "aS1(6,-1)" dd 
getvariable.pl -c 0.000000  "aS1(6,0)" dd 
getvariable.pl -c 0.007933  "aS1(6,1)" dd 
getvariable.pl -c 0.000000  "aS1(6,2)" dd 
getvariable.pl -c 0.001331  "aS1(6,3)" dd 
getvariable.pl -c 0.000000  "aS1(6,4)" dd 
getvariable.pl -c 0.071479  "aS1(6,5)" dd 
getvariable.pl -c 0.000000  "aS1(6,6)" dd 
getvariable.pl -c 0.000000  "aS2(0,0)" dd 
getvariable.pl -c 0.000000  "aS2(1,-1)" dd 
getvariable.pl -c 0.000000  "aS2(1,0)" dd 
getvariable.pl -c 0.000000  "aS2(1,1)" dd 
getvariable.pl -c 0.000000  "aS2(2,-2)" dd 
getvariable.pl -c 0.003722  "aS2(2,-1)" dd 
getvariable.pl -c 0.000000  "aS2(2,0)" dd 
getvariable.pl -c 0.000000  "aS2(2,1)" dd 
getvariable.pl -c 0.000000  "aS2(2,2)" dd 
getvariable.pl -c 0.000000  "aS2(3,-3)" dd 
getvariable.pl -c 0.000000  "aS2(3,-2)" dd 
getvariable.pl -c 0.000000  "aS2(3,-1)" dd 
getvariable.pl -c 0.000000  "aS2(3,0)" dd 
getvariable.pl -c 0.000000  "aS2(3,1)" dd 
getvariable.pl -c 0.000000  "aS2(3,2)" dd 
getvariable.pl -c 0.000000  "aS2(3,3)" dd 
getvariable.pl -c 0.000000  "aS2(4,-4)" dd 
getvariable.pl -c -0.006511  "aS2(4,-3)" dd 
getvariable.pl -c 0.000000  "aS2(4,-2)" dd 
getvariable.pl -c 0.070528  "aS2(4,-1)" dd 
getvariable.pl -c 0.000000  "aS2(4,0)" dd 
getvariable.pl -c 0.000000  "aS2(4,1)" dd 
getvariable.pl -c 0.000000  "aS2(4,2)" dd 
getvariable.pl -c 0.000000  "aS2(4,3)" dd 
getvariable.pl -c 0.000000  "aS2(4,4)" dd 
getvariable.pl -c 0.000000  "aS2(5,-5)" dd 
getvariable.pl -c 0.000000  "aS2(5,-4)" dd 
getvariable.pl -c 0.000000  "aS2(5,-3)" dd 
getvariable.pl -c 0.000000  "aS2(5,-2)" dd 
getvariable.pl -c 0.000000  "aS2(5,-1)" dd 
getvariable.pl -c 0.000000  "aS2(5,0)" dd 
getvariable.pl -c 0.000000  "aS2(5,1)" dd 
getvariable.pl -c 0.000000  "aS2(5,2)" dd 
getvariable.pl -c 0.000000  "aS2(5,3)" dd 
getvariable.pl -c 0.000000  "aS2(5,4)" dd 
getvariable.pl -c 0.000000  "aS2(5,5)" dd 
getvariable.pl -c 0.000000  "aS2(6,-6)" dd 
getvariable.pl -c 0.071479  "aS2(6,-5)" dd 
getvariable.pl -c 0.000000  "aS2(6,-4)" dd 
getvariable.pl -c -0.001331  "aS2(6,-3)" dd 
getvariable.pl -c 0.000000  "aS2(6,-2)" dd 
getvariable.pl -c 0.007933  "aS2(6,-1)" dd 
getvariable.pl -c 0.000000  "aS2(6,0)" dd 
getvariable.pl -c 0.000000  "aS2(6,1)" dd 
getvariable.pl -c 0.000000  "aS2(6,2)" dd 
getvariable.pl -c 0.000000  "aS2(6,3)" dd 
getvariable.pl -c 0.000000  "aS2(6,4)" dd 
getvariable.pl -c 0.000000  "aS2(6,5)" dd 
getvariable.pl -c 0.000000  "aS2(6,6)" dd 
getvariable.pl -c 0.272735  "aS3(0,0)" dd 
getvariable.pl -c 0.000000  "aS3(1,-1)" dd 
getvariable.pl -c 0.000000  "aS3(1,0)" dd 
getvariable.pl -c 0.000000  "aS3(1,1)" dd 
getvariable.pl -c 0.000000  "aS3(2,-2)" dd 
getvariable.pl -c 0.000000  "aS3(2,-1)" dd 
getvariable.pl -c 0.029124  "aS3(2,0)" dd 
getvariable.pl -c 0.000000  "aS3(2,1)" dd 
getvariable.pl -c 0.000000  "aS3(2,2)" dd 
getvariable.pl -c 0.000000  "aS3(3,-3)" dd 
getvariable.pl -c 0.000000  "aS3(3,-2)" dd 
getvariable.pl -c 0.000000  "aS3(3,-1)" dd 
getvariable.pl -c 0.000000  "aS3(3,0)" dd 
getvariable.pl -c 0.000000  "aS3(3,1)" dd 
getvariable.pl -c 0.000000  "aS3(3,2)" dd 
getvariable.pl -c 0.000000  "aS3(3,3)" dd 
getvariable.pl -c 0.000000  "aS3(4,-4)" dd 
getvariable.pl -c 0.000000  "aS3(4,-3)" dd 
getvariable.pl -c 0.000000  "aS3(4,-2)" dd 
getvariable.pl -c 0.000000  "aS3(4,-1)" dd 
getvariable.pl -c 0.090362  "aS3(4,0)" dd 
getvariable.pl -c 0.000000  "aS3(4,1)" dd 
getvariable.pl -c 0.000000  "aS3(4,2)" dd 
getvariable.pl -c 0.000000  "aS3(4,3)" dd 
getvariable.pl -c -0.022049  "aS3(4,4)" dd 
getvariable.pl -c 0.000000  "aS3(5,-5)" dd 
getvariable.pl -c 0.000000  "aS3(5,-4)" dd 
getvariable.pl -c 0.000000  "aS3(5,-3)" dd 
getvariable.pl -c 0.000000  "aS3(5,-2)" dd 
getvariable.pl -c 0.000000  "aS3(5,-1)" dd 
getvariable.pl -c 0.000000  "aS3(5,0)" dd 
getvariable.pl -c 0.000000  "aS3(5,1)" dd 
getvariable.pl -c 0.000000  "aS3(5,2)" dd 
getvariable.pl -c 0.000000  "aS3(5,3)" dd 
getvariable.pl -c 0.000000  "aS3(5,4)" dd 
getvariable.pl -c 0.000000  "aS3(5,5)" dd 
getvariable.pl -c 0.000000  "aS3(6,-6)" dd 
getvariable.pl -c 0.000000  "aS3(6,-5)" dd 
getvariable.pl -c 0.000000  "aS3(6,-4)" dd 
getvariable.pl -c 0.000000  "aS3(6,-3)" dd 
getvariable.pl -c 0.000000  "aS3(6,-2)" dd 
getvariable.pl -c 0.000000  "aS3(6,-1)" dd 
getvariable.pl -c -0.001015  "aS3(6,0)" dd 
getvariable.pl -c 0.000000  "aS3(6,1)" dd 
getvariable.pl -c 0.000000  "aS3(6,2)" dd 
getvariable.pl -c 0.000000  "aS3(6,3)" dd 
getvariable.pl -c 0.076607  "aS3(6,4)" dd 
getvariable.pl -c 0.000000  "aS3(6,5)" dd 
getvariable.pl -c 0.000000  "aS3(6,6)" dd 

setup_mcdisp_mf -prefix test_ 1 0 0 81
mcdispit -prefix test_ -maxE 100

getvalue.pl -c   0.0483285 9 10 6.38398 0 results/test_mcdisp.qei
getvalue.pl -c    0.0501943 9 11 6.38398 0  results/test_mcdisp.qei

spins  -c -M -prefix test_ 1 0 0 81 0 0 1 6.38398

spins  -s -M -prefix test_ 1 0 0 81 0 0 1 6.38398

rm dd
cd ../../demo
