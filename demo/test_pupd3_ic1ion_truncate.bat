cd ../examples/pupd3

REM in mcphas 5.5 there was definitely something buggy in
REM output - see emails with Duc autumn 2024 and some bugs were removed
REM (ic1ion and icf1ion were not consistent) 
REM 16.1.2025 tested with respect to XTLS on pu3p.sipf 
REM with alpha=beta=gamma=0  see 

# 200 sec
ic1ion pu3p.sipf
getvariable.pl -c -7496 E0 results/ic1ion.out
getvalue.pl -c 0  0 1 1 0  results/ic1ion.out
getvalue.pl -c 0  0 1 2 0  results/ic1ion.out
getvalue.pl -c 12.24  0 1 3 0  results/ic1ion.out
getvalue.pl -c 12.24  0 1 4 0  results/ic1ion.out
getvalue.pl -c 12.24  0 1 5 0  results/ic1ion.out
getvalue.pl -c 12.24  0 1 6 0  results/ic1ion.out
getvalue.pl -c 321.7  0 1 7 0  results/ic1ion.out
getvalue.pl -c 321.7  0 1 8 0  results/ic1ion.out
getvalue.pl -c 490  0 1 9 0  results/ic1ion.out
getvalue.pl -c 490  0 1 10 0  results/ic1ion.out
getvalue.pl -c 540  0 1 11 0  results/ic1ion.out
getvalue.pl -c 540  0 1 12 0  results/ic1ion.out
getvalue.pl -c 540  0 1 13 0  results/ic1ion.out
getvalue.pl -c 540  0 1 14 0  results/ic1ion.out



icf1ion pu3p.sipf
getvariable.pl -c -1122 E0 results/icf1ion.out
getvalue.pl -c 0  0 1 1 0  results/icf1ion.out
getvalue.pl -c 0  0 1 2 0  results/icf1ion.out
getvalue.pl -c 8.31  0 1 3 0  results/icf1ion.out
getvalue.pl -c 8.31  0 1 4 0  results/icf1ion.out
getvalue.pl -c 8.31  0 1 5 0  results/icf1ion.out
getvalue.pl -c 8.31  0 1 6 0  results/icf1ion.out
getvalue.pl -c 81.9  0 1 7 0  results/icf1ion.out
getvalue.pl -c 81.9  0 1 8 0  results/icf1ion.out
getvalue.pl -c 374  0 1 9 0  results/icf1ion.out
getvalue.pl -c 374  0 1 10 0  results/icf1ion.out
getvalue.pl -c 381  0 1 11 0  results/icf1ion.out
getvalue.pl -c 381  0 1 12 0  results/icf1ion.out
getvalue.pl -c 458  0 1 13 0  results/icf1ion.out
getvalue.pl -c 458  0 1 14 0  results/icf1ion.out

# 100sec
mcphasit -prefix test_

getvalue.pl -c -7497.1 2 8 1 0 results/test_mcphas.fum
getvalue.pl -c 0.025 2 10 21 0 results/test_mcphas.fum

#300sec
spins -c -prefix test_ 1 0 0 81 > dd
getvariable.pl -c 1.41050   "a(0,0)" dd 
getvariable.pl -c 0.000000   "a(2,-2)" dd 
getvariable.pl -c 0.000000   "a(2,-1)" dd 
getvariable.pl -c 0.00527   "a(2,0)" dd 
getvariable.pl -c 0.000000   "a(2,1)" dd 
getvariable.pl -c 0.000000   "a(2,2)" dd 
getvariable.pl -c 0.000000   "a(4,-4)" dd 
getvariable.pl -c 0.000000   "a(4,-3)" dd 
getvariable.pl -c 0.000000   "a(4,-2)" dd 
getvariable.pl -c 0.000000   "a(4,-1)" dd 
getvariable.pl -c 0.02813   "a(4,0)" dd 
getvariable.pl -c 0.000000   "a(4,1)" dd 
getvariable.pl -c 0.000000   "a(4,2)" dd 
getvariable.pl -c 0.000000   "a(4,3)" dd 
getvariable.pl -c 0.02584   "a(4,4)" dd 
getvariable.pl -c 0.000000   "a(6,-6)" dd 
getvariable.pl -c 0.000000   "a(6,-5)" dd 
getvariable.pl -c 0.000000   "a(6,-4)" dd 
getvariable.pl -c 0.000000   "a(6,-3)" dd 
getvariable.pl -c 0.000000   "a(6,-2)" dd 
getvariable.pl -c 0.000000   "a(6,-1)" dd 
getvariable.pl -c -0.08879   "a(6,0)" dd 
getvariable.pl -c 0.000000   "a(6,1)" dd 
getvariable.pl -c 0.000000   "a(6,2)" dd 
getvariable.pl -c 0.000000   "a(6,3)" dd 
getvariable.pl -c 0.24049   "a(6,4)" dd 
getvariable.pl -c 0.000000   "a(6,5)" dd 
getvariable.pl -c 0.000000   "a(6,6)" dd 

# 800sec
spins -s -M -prefix test_ 1 0 0 81 > dd

getvariable.pl -c 0.000000  "aS1(0,0)" dd 
getvariable.pl -c 0.000000  "aS1(1,-1)" dd 
getvariable.pl -c 0.000000  "aS1(1,0)" dd 
getvariable.pl -c 0.000000  "aS1(1,1)" dd 
getvariable.pl -c 0.000000  "aS1(2,-2)" dd 
getvariable.pl -c 0.000000  "aS1(2,-1)" dd 
getvariable.pl -c 0.000000  "aS1(2,0)" dd 
getvariable.pl -c 0.02540  "aS1(2,1)" dd 
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
getvariable.pl -c 0.07984  "aS1(4,1)" dd 
getvariable.pl -c 0.000000  "aS1(4,2)" dd 
getvariable.pl -c 0.00525  "aS1(4,3)" dd 
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
getvariable.pl -c 0.00393  "aS1(6,1)" dd 
getvariable.pl -c 0.000000  "aS1(6,2)" dd 
getvariable.pl -c 0.00703  "aS1(6,3)" dd 
getvariable.pl -c 0.000000  "aS1(6,4)" dd 
getvariable.pl -c 0.08258  "aS1(6,5)" dd 
getvariable.pl -c 0.000000  "aS1(6,6)" dd 
getvariable.pl -c 0.000000  "aS2(0,0)" dd 
getvariable.pl -c 0.000000  "aS2(1,-1)" dd 
getvariable.pl -c 0.000000  "aS2(1,0)" dd 
getvariable.pl -c 0.000000  "aS2(1,1)" dd 
getvariable.pl -c 0.000000  "aS2(2,-2)" dd 
getvariable.pl -c 0.02540  "aS2(2,-1)" dd 
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
getvariable.pl -c -0.00525  "aS2(4,-3)" dd 
getvariable.pl -c 0.000000  "aS2(4,-2)" dd 
getvariable.pl -c 0.07983  "aS2(4,-1)" dd 
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
getvariable.pl -c 0.08258  "aS2(6,-5)" dd 
getvariable.pl -c 0.000000  "aS2(6,-4)" dd 
getvariable.pl -c -0.00703  "aS2(6,-3)" dd 
getvariable.pl -c 0.000000  "aS2(6,-2)" dd 
getvariable.pl -c 0.00393  "aS2(6,-1)" dd 
getvariable.pl -c 0.000000  "aS2(6,0)" dd 
getvariable.pl -c 0.000000  "aS2(6,1)" dd 
getvariable.pl -c 0.000000  "aS2(6,2)" dd 
getvariable.pl -c 0.000000  "aS2(6,3)" dd 
getvariable.pl -c 0.000000  "aS2(6,4)" dd 
getvariable.pl -c 0.000000  "aS2(6,5)" dd 
getvariable.pl -c 0.000000  "aS2(6,6)" dd 
getvariable.pl -c 0.25685  "aS3(0,0)" dd 
getvariable.pl -c 0.000000  "aS3(1,-1)" dd 
getvariable.pl -c 0.000000  "aS3(1,0)" dd 
getvariable.pl -c 0.000000  "aS3(1,1)" dd 
getvariable.pl -c 0.000000  "aS3(2,-2)" dd 
getvariable.pl -c 0.000000  "aS3(2,-1)" dd 
getvariable.pl -c 0.01844  "aS3(2,0)" dd 
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
getvariable.pl -c 0.09719  "aS3(4,0)" dd 
getvariable.pl -c 0.000000  "aS3(4,1)" dd 
getvariable.pl -c 0.000000  "aS3(4,2)" dd 
getvariable.pl -c 0.000000  "aS3(4,3)" dd 
getvariable.pl -c  -0.0325  "aS3(4,4)" dd 
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
getvariable.pl -c -0.01255  "aS3(6,0)" dd 
getvariable.pl -c 0.000000  "aS3(6,1)" dd 
getvariable.pl -c 0.000000  "aS3(6,2)" dd 
getvariable.pl -c 0.000000  "aS3(6,3)" dd 
getvariable.pl -c 0.07758  "aS3(6,4)" dd 
getvariable.pl -c 0.000000  "aS3(6,5)" dd 
getvariable.pl -c 0.000000  "aS3(6,6)" dd 

# 100 sec
setup_mcdisp_mf -prefix test_ 1 0 0 81
mcdispit -prefix test_ -maxE 100

getvalue.pl -c   0.04852 9 10 6.33994 0 results/test_mcdisp.qei
getvalue.pl -c    0.04977 9 11 6.33994 0  results/test_mcdisp.qei

#display_densities  -c -M -prefix test_ 1 0 0 81 0 0 1 6.38398

#display_densities  -s -M -prefix test_ 1 0 0 81 0 0 1 6.38398

rm dd
cd ../../demo

