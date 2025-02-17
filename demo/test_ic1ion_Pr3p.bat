cd ../examples/testic1ion
call pointc testic1ion.sipf pointcharges_octahedron.dat > octahedronic1ion.sipf
getvariable.pl -c 91.118 L64 results/pointc.Llm
call substitute "MODULE=so1ion" "MODULE=ic1ion" octahedronic1ion.sipf

call densplt c octahedronic1ion.sipf 2 0 0 0 
call ic1ion octahedronic1ion.sipf
getvalue.pl -c 188.112 0 1 9 0 results/ic1ion.out
densplt c -M Co2p_atom1_rotated_pm_z.sipf 2 0 0 100
singleion -M -r Co2p_atom1_rotated_pm_z.sipf 2 10 0 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 0.06930 3 12 10 0 m.clc
singleion -M -r Co2p_atom1_rotated_pm_z.sipf 2 0 10 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 0.0698 4 13 10 0 m.clc
singleion -M -r Co2p_atom1_rotated_pm_z.sipf 2 0 0 10  0 0 0 0 0 0 > m.clc
getvalue.pl -c 5.998 5 14 10 0 m.clc

call substitute "MODULE=ic1ion" "MODULE=icf1ion" octahedronic1ion.sipf
call densplt c octahedronic1ion.sipf 2 0 0 0 
call icf1ion octahedronic1ion.sipf
getvalue.pl -c 173.306 0 1 9 0 results/icf1ion.out
singleion -M -r octahedronic1ion.sipf -Tsteps 20 400 2 10 0 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 0.2384 2 12  81.6 0 m.clc
singleion -M -r octahedronic1ion.sipf 2 0 10 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 1.9939 4 13 10 0 m.clc
singleion -M -r octahedronic1ion.sipf 2 0 0 10  0 0 0 0 0 0 > m.clc
getvalue.pl -c 1.9939 5 14 10 0 m.clc


rem compare singleion and dnsplt with ic1ion and icf1ion for the same very assymetrical problem
singleion -r Ce3p_xyz_ic1ion.sipf 2 20 10 30   0 0 0 0 0 0   0 0 0 > m.clc
getvalue.pl -c -0.20326 2 15 2 0 m.clc
getvalue.pl -c -0.1016  2 16 2 0 m.clc
getvalue.pl -c -0.3049  2 17 2 0 m.clc
getvalue.pl -c 1.5396  2 18 2 0 m.clc
getvalue.pl -c 0.76978  2 19 2 0 m.clc
getvalue.pl -c 2.3094  2 20 2 0 m.clc
getvalue.pl -c -0.10268  2 21 2 0 m.clc
getvalue.pl -c -0.15402  2 22 2 0 m.clc
getvalue.pl -c -0.13623  2 23 2 0 m.clc
getvalue.pl -c -1e-06  2 24 2 0 m.clc
getvalue.pl -c 9.1716  2 25 2 0 m.clc
getvalue.pl -c 35.337  2 26 2 0 m.clc
getvalue.pl -c 40.78  2 27 2 0 m.clc
getvalue.pl -c 56.154 2 28 2 0 m.clc

getvalue.pl -c 35.337       0 6 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 3.372      0 7 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.0910267  0 8 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0           0 9 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 2            0 10 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0158445    0 11 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.0205977   0 12 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.00792222  0 13 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   1.14772    0 14 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   1.49202    0 15 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.573855    0 16 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.02728    0 17 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0525    0 18 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0341042    0 19 4 0 results/Ce3p_xyz_ic1ion.sipf.trs

getvalue.pl -c 0.0       0 6 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 7 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 8 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 9 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 10 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 11 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 12 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 13 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 14 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 15 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 16 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 17 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 18 1 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0       0 19 1 0 results/Ce3p_xyz_ic1ion.sipf.trs



densplt c Ce3p_xyz_ic1ion.sipf 2 20 10 30 > m.clc
getvariable.pl -c 0.282 "a(0,0)" m.clc 
getvariable.pl -c -0.046 "a(2,-2)" m.clc 
getvariable.pl -c -0.0687 "a(2,-1)" m.clc 
getvariable.pl -c -0.0859 "a(2,0)" m.clc 
getvariable.pl -c -0.1374 "a(2,1)" m.clc 
getvariable.pl -c -0.0343 "a(2,2)" m.clc 
getvariable.pl -c 0.0042 "a(4,-4)" m.clc 
getvariable.pl -c 0.0163 "a(4,-3)" m.clc 
getvariable.pl -c 0.0258 "a(4,-2)" m.clc 
getvariable.pl -c 0.0117 "a(4,-1)" m.clc 
getvariable.pl -c -0.0105 "a(4,0)" m.clc 
getvariable.pl -c 0.0235 "a(4,1)" m.clc 
getvariable.pl -c 0.0194 "a(4,2)" m.clc 
getvariable.pl -c 0.00296 "a(4,3)" m.clc 
getvariable.pl -c -0.0012 "a(4,4)" m.clc 
getvariable.pl -c -0.00002 "a(6,-6)" m.clc 
getvariable.pl -c -0.00020 "a(6,-5)" m.clc 
getvariable.pl -c -0.0007 "a(6,-4)" m.clc 
getvariable.pl -c -0.00119 "a(6,-3)" m.clc 
getvariable.pl -c -0.00076 "a(6,-2)" m.clc 
getvariable.pl -c 0.00015 "a(6,-1)" m.clc 
getvariable.pl -c 0.00076 "a(6,0)" m.clc 
getvariable.pl -c 0.00031 "a(6,1)" m.clc 
getvariable.pl -c -0.00057 "a(6,2)" m.clc 
getvariable.pl -c -0.00022 "a(6,3)" m.clc 
getvariable.pl -c 0.00021 "a(6,4)" m.clc 
getvariable.pl -c 0.00019 "a(6,5)" m.clc 
getvariable.pl -c 0.00006 "a(6,6)" m.clc


singleion -r Ce3p_xyz_icf1ion.sipf 2 20 10 30   0 0 0 0 0 0   0 0 0 > m.clc
getvalue.pl -c -0.20326 2 15 2 0 m.clc
getvalue.pl -c -0.1016  2 16 2 0 m.clc
getvalue.pl -c -0.3049  2 17 2 0 m.clc
getvalue.pl -c 1.5396  2 18 2 0 m.clc
getvalue.pl -c 0.76978  2 19 2 0 m.clc
getvalue.pl -c 2.3094  2 20 2 0 m.clc
getvalue.pl -c -0.10268  2 21 2 0 m.clc
getvalue.pl -c -0.15402  2 22 2 0 m.clc
getvalue.pl -c -0.13623  2 23 2 0 m.clc
getvalue.pl -c -1e-06  2 24 2 0 m.clc
getvalue.pl -c 9.1716  2 25 2 0 m.clc
getvalue.pl -c 35.337  2 26 2 0 m.clc
getvalue.pl -c 40.78  2 27 2 0 m.clc
getvalue.pl -c 56.154 2 28 2 0 m.clc

getvalue.pl -c 35.337       0 6 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 3.372        0 7 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0.0910267  0 8 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0           0 9 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 2            0 10 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.0158445    0 11 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.0205977   0 12 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.0079222   0 13 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   1.14772    0 14 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   1.49202    0 15 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.573855    0 16 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.02728      0 17 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.0525       0 18 4 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.0341042    0 19 4 0 results/Ce3p_xyz_icf1ion.sipf.trs

getvalue.pl -c 0       0 6 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0       0 7 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0       0 8 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0      0 9 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0       0 10 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0       0 11 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0      0 12 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0      0 13 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0      0 14 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0      0 15 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0      0 16 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0       0 17 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0       0 18 1 0 results/Ce3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0       0 19 1 0 results/Ce3p_xyz_icf1ion.sipf.trs


densplt c Ce3p_xyz_icf1ion.sipf 2 20 10 30 > m.clc
getvariable.pl -c 0.282099 "a(0,0)" m.clc 
getvariable.pl -c -0.045797 "a(2,-2)" m.clc 
getvariable.pl -c -0.068695 "a(2,-1)" m.clc 
getvariable.pl -c -0.085934 "a(2,0)" m.clc 
getvariable.pl -c -0.137391 "a(2,1)" m.clc 
getvariable.pl -c -0.034348 "a(2,2)" m.clc 
getvariable.pl -c 0.004186 "a(4,-4)" m.clc 
getvariable.pl -c 0.016281 "a(4,-3)" m.clc 
getvariable.pl -c 0.025843 "a(4,-2)" m.clc 
getvariable.pl -c 0.011748 "a(4,-1)" m.clc 
getvariable.pl -c -0.010525 "a(4,0)" m.clc 
getvariable.pl -c 0.023495 "a(4,1)" m.clc 
getvariable.pl -c 0.019383 "a(4,2)" m.clc 
getvariable.pl -c 0.002960 "a(4,3)" m.clc 
getvariable.pl -c -0.001221 "a(4,4)" m.clc 
getvariable.pl -c -0.000021 "a(6,-6)" m.clc 
getvariable.pl -c -0.000200 "a(6,-5)" m.clc 
getvariable.pl -c -0.000708 "a(6,-4)" m.clc 
getvariable.pl -c -0.001191 "a(6,-3)" m.clc 
getvariable.pl -c -0.000761 "a(6,-2)" m.clc 
getvariable.pl -c 0.000153 "a(6,-1)" m.clc 
getvariable.pl -c 0.000757 "a(6,0)" m.clc 
getvariable.pl -c 0.000305 "a(6,1)" m.clc 
getvariable.pl -c -0.000571 "a(6,2)" m.clc 
getvariable.pl -c -0.000217 "a(6,3)" m.clc 
getvariable.pl -c 0.000206 "a(6,4)" m.clc 
getvariable.pl -c 0.000186 "a(6,5)" m.clc 
getvariable.pl -c 0.000055 "a(6,6)" m.clc

ic1ion Nd3p_xyz_ic1ion.sipf
getvalue.pl -c 0.7923  4 7 1.73205 0 results/icf1ion.mag 

singleion -r Nd3p_xyz_ic1ion.sipf 2 20 10 30   0 0 0 0 0 0   0 0 0 > m.clc

getvalue.pl -c -0.707389  2 15 2 0  m.clc
getvalue.pl -c  -0.353694  2 16 2 0  m.clc
getvalue.pl -c  -1.06109  2 17 2 0  m.clc
getvalue.pl -c  3.11274  2 18 2 0  m.clc
getvalue.pl -c  1.55637  2 19 2 0  m.clc
getvalue.pl -c  4.66912  2 20 2 0  m.clc
getvalue.pl -c  -0.0416526  2 21 2 0  m.clc
getvalue.pl -c  -0.062478  2 22 2 0  m.clc
getvalue.pl -c  -0.0552648  2 23 2 0  m.clc
getvalue.pl -c  -1e-06  2 24 2 0  m.clc
getvalue.pl -c  13.7361  2 25 2 0  m.clc
getvalue.pl -c  80.0232  2 26 2 0  m.clc
getvalue.pl -c  90.8995  2 27 2 0  m.clc
getvalue.pl -c  95.7337   2 28 2 0  m.clc

  
getvalue.pl -c 80.0232        0 6 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   6.63754       0 7 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c     0.137944   0 8 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    0           0 9 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    2           0 10 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    0.064627    0 11 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.0840141    0 12 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.0323135    0 13 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    2.30302     0 14 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    2.99391     0 15 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    1.15151     0 16 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.001951   0 17 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.003758  0 18 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.0024391    0 19 4 0 results/Nd3p_xyz_ic1ion.sipf.trs

icf1ion Nd3p_xyz_icf1ion.sipf
getvalue.pl -c 0.792305 4 7 1.73205 0 results/icf1ion.mag 

densplt c Nd3p_xyz_icf1ion.sipf 2 20 10 30 > m.clc

getvariable.pl -c 0.8463 "a(0,0)" m.clc 
getvariable.pl -c -0.01939 "a(2,-2)" m.clc 
getvariable.pl -c -0.02908 "a(2,-1)" m.clc 
getvariable.pl -c -0.03638 "a(2,0)" m.clc 
getvariable.pl -c -0.05816 "a(2,1)" m.clc 
getvariable.pl -c -0.01454 "a(2,2)" m.clc 
getvariable.pl -c -0.00549 "a(4,-4)" m.clc 
getvariable.pl -c -0.02135 "a(4,-3)" m.clc 
getvariable.pl -c -0.03389 "a(4,-2)" m.clc 
getvariable.pl -c -0.01541 "a(4,-1)" m.clc 
getvariable.pl -c 0.01380 "a(4,0)" m.clc 
getvariable.pl -c -0.03081 "a(4,1)" m.clc 
getvariable.pl -c -0.02542 "a(4,2)" m.clc 
getvariable.pl -c -0.00388 "a(4,3)" m.clc 
getvariable.pl -c 0.00160 "a(4,4)" m.clc 
getvariable.pl -c -0.00077 "a(6,-6)" m.clc 
getvariable.pl -c -0.00746 "a(6,-5)" m.clc 
getvariable.pl -c -0.02639 "a(6,-4)" m.clc 
getvariable.pl -c -0.04443 "a(6,-3)" m.clc 
getvariable.pl -c -0.02839 "a(6,-2)" m.clc 
getvariable.pl -c 0.00569 "a(6,-1)" m.clc 
getvariable.pl -c 0.02822 "a(6,0)" m.clc 
getvariable.pl -c 0.01138 "a(6,1)" m.clc 
getvariable.pl -c -0.0213 "a(6,2)" m.clc 
getvariable.pl -c -0.00808 "a(6,3)" m.clc 
getvariable.pl -c 0.0077 "a(6,4)" m.clc 
getvariable.pl -c 0.00692 "a(6,5)" m.clc 
getvariable.pl -c 0.00205 "a(6,6)" m.clc 


   

singleion -r Nd3p_xyz_icf1ion.sipf 2 20 10 30   0 0 0 0 0 0   0 0 0 > m.clc

getvalue.pl -c -0.72716  2 15 2 0  m.clc
getvalue.pl -c  -0.36358 2 16 2 0  m.clc
getvalue.pl -c   -1.0908 2 17 2 0  m.clc
getvalue.pl -c   3.1325  2 18 2 0  m.clc
getvalue.pl -c   1.5663  2 19 2 0  m.clc
getvalue.pl -c   4.6988  2 20 2 0  m.clc
getvalue.pl -c   -0.043466  2 21 2 0  m.clc
getvalue.pl -c   -0.0652  2 22 2 0  m.clc
getvalue.pl -c   -0.057672  2 23 2 0  m.clc
getvalue.pl -c   -1e-06   2 24 2 0  m.clc
getvalue.pl -c   13.576  2 25 2 0  m.clc
getvalue.pl -c   86.965  2 26 2 0  m.clc
getvalue.pl -c   97.606  2 27 2 0  m.clc
getvalue.pl -c   111.2  2 28 2 0  m.clc

getvalue.pl -c   86.9645    0 6 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   6.64096    0 7 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c    0.140107  0 8 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0          0 9 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   2          0 10 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0.0628451  0 11 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0.081698   0 12 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0.0314225  0 13 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c    2.30548   0 14 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   2.99711    0 15 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   1.15274    0 16 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0.0023 0 17 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0.0045 0 18 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.0029  0 19 4 0 results/Nd3p_xyz_icf1ion.sipf.trs



densplt c Nd3p_xyz_ic1ion.sipf 2 20 10 30 > m.clc

getvariable.pl -c 0.85 "a(0,0)" m.clc 
getvariable.pl -c -0.019 "a(2,-2)" m.clc 
getvariable.pl -c -0.028 "a(2,-1)" m.clc 
getvariable.pl -c -0.035 "a(2,0)" m.clc 
getvariable.pl -c -0.056 "a(2,1)" m.clc 
getvariable.pl -c -0.014 "a(2,2)" m.clc 
getvariable.pl -c -0.0052 "a(4,-4)" m.clc 
getvariable.pl -c -0.020 "a(4,-3)" m.clc 
getvariable.pl -c -0.032 "a(4,-2)" m.clc 
getvariable.pl -c -0.015 "a(4,-1)" m.clc 
getvariable.pl -c 0.013 "a(4,0)" m.clc 
getvariable.pl -c -0.03 "a(4,1)" m.clc 
getvariable.pl -c -0.024 "a(4,2)" m.clc 
getvariable.pl -c -0.004 "a(4,3)" m.clc 
getvariable.pl -c 0.0015 "a(4,4)" m.clc 
getvariable.pl -c -0.0007 "a(6,-6)" m.clc 
getvariable.pl -c -0.007 "a(6,-5)" m.clc 
getvariable.pl -c -0.024 "a(6,-4)" m.clc 
getvariable.pl -c -0.041 "a(6,-3)" m.clc 
getvariable.pl -c -0.026 "a(6,-2)" m.clc 
getvariable.pl -c 0.0053 "a(6,-1)" m.clc 
getvariable.pl -c 0.026 "a(6,0)" m.clc 
getvariable.pl -c 0.011 "a(6,1)" m.clc 
getvariable.pl -c -0.020 "a(6,2)" m.clc 
getvariable.pl -c -0.0075 "a(6,3)" m.clc 
getvariable.pl -c 0.0071 "a(6,4)" m.clc 
getvariable.pl -c 0.0064 "a(6,5)" m.clc 
getvariable.pl -c 0.0019 "a(6,6)" m.clc

singleion -r Pr3p_xyz.sipf -MQ 1 0 0 2 20 10 30 0 0 0  > m.clc 

getvalue.pl -c  1.68  2 9 2 0 m.clc
getvalue.pl -c  1.68  2 10 2 0 m.clc
getvalue.pl -c  1.68  2 12 2 0 m.clc
getvalue.pl -c  0.84  2 13 2 0 m.clc
getvalue.pl -c  0.84  2 14 2 0 m.clc
getvalue.pl -c  0.84  2 16 2 0 m.clc
getvalue.pl -c  2.52   2 17 2 0 m.clc

getvalue.pl -c 0 0 8 1 0    results/Pr3p_xyz.sipf.trs
getvalue.pl -c 0   0 11 1 0 results/Pr3p_xyz.sipf.trs
getvalue.pl -c 0  0 12 1 0  results/Pr3p_xyz.sipf.trs
getvalue.pl -c 0  0 13 1 0  results/Pr3p_xyz.sipf.trs

getvalue.pl -c 0.12  0 8 4 0 results/Pr3p_xyz.sipf.trs
getvalue.pl -c 0.88  0 11 4 0 results/Pr3p_xyz.sipf.trs
getvalue.pl -c 1.15  0 12 4 0 results/Pr3p_xyz.sipf.trs
getvalue.pl -c 0.45  0 13 4 0 results/Pr3p_xyz.sipf.trs

singleion -r Pr3p_xyz_icf1ion.sipf -MQ 1 0 0 2 20 10 30 0 0 0  > m.clc 

getvalue.pl -c  1.67  2 9 2 0 m.clc
getvalue.pl -c  1.67  2 10 2 0 m.clc
getvalue.pl -c  1.61  2 12 2 0 m.clc
getvalue.pl -c  0.83  2 13 2 0 m.clc
getvalue.pl -c  0.83  2 14 2 0 m.clc
getvalue.pl -c   0.81  2 16 2 0 m.clc
getvalue.pl -c  2.50   2 17 2 0 m.clc

getvalue.pl -c 0 0 8 1 0    results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0   0 11 1 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0  0 12 1 0  results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0  0 13 1 0  results/Pr3p_xyz_icf1ion.sipf.trs

getvalue.pl -c 0.13  0 8 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.90   0 11 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 1.18  0 12 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.46  0 13 4 0 results/Pr3p_xyz_icf1ion.sipf.trs

singleion -r Pr3p_xyz_ic1ion.sipf -MQ 1 0 0 2 20 10 30 0 0 0  > m.clc 

getvalue.pl -c  1.69  2 9 2 0 m.clc
getvalue.pl -c  1.69  2 10 2 0 m.clc
getvalue.pl -c  1.63   2 12 2 0 m.clc
getvalue.pl -c  0.84  2 13 2 0 m.clc
getvalue.pl -c  0.84  2 14 2 0 m.clc
getvalue.pl -c  0.81  2 16 2 0 m.clc
getvalue.pl -c  2.52  2 17 2 0 m.clc

getvalue.pl -c 0 0 8 1 0    results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0   0 11 1 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0  0 12 1 0  results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0  0 13 1 0  results/Pr3p_xyz_ic1ion.sipf.trs

getvalue.pl -c 0.13  0 8 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.92  0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 1.20  0 12 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.46  0 13 4 0 results/Pr3p_xyz_ic1ion.sipf.trs


densplt s Pr3p_xyz_ic1ion.sipf 2 20 10 30 > m.clc

getvariable.pl -c -0.24 "aS1(0,0)" m.clc 
getvariable.pl -c 0.000000 "aS1(1,-1)" m.clc 
getvariable.pl -c 0.000000 "aS1(1,0)" m.clc 
getvariable.pl -c 0.000000 "aS1(1,1)" m.clc 
getvariable.pl -c 0.040 "aS1(2,-2)" m.clc 
getvariable.pl -c 0.014 "aS1(2,-1)" m.clc 
getvariable.pl -c -0.019 "aS1(2,0)" m.clc 
getvariable.pl -c 0.12 "aS1(2,1)" m.clc 
getvariable.pl -c 0.069 "aS1(2,2)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,-3)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,-2)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,-1)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,0)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,1)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,2)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,3)" m.clc 
getvariable.pl -c 0.0097 "aS1(4,-4)" m.clc 
getvariable.pl -c 0.024 "aS1(4,-3)" m.clc 
getvariable.pl -c 0.019 "aS1(4,-2)" m.clc 
getvariable.pl -c -0.0042 "aS1(4,-1)" m.clc 
getvariable.pl -c -0.027 "aS1(4,0)" m.clc 
getvariable.pl -c 0.0047 "aS1(4,1)" m.clc 
getvariable.pl -c 0.032 "aS1(4,2)" m.clc 
getvariable.pl -c 0.016 "aS1(4,3)" m.clc 
getvariable.pl -c 0.0012 "aS1(4,4)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,-5)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,-4)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,-3)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,-2)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,-1)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,0)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,1)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,2)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,3)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,4)" m.clc 
getvariable.pl -c 0.000000 "aS1(5,5)" m.clc 
getvariable.pl -c -0.0023 "aS1(6,-6)" m.clc 
getvariable.pl -c -0.012 "aS1(6,-5)" m.clc 
getvariable.pl -c -0.025 "aS1(6,-4)" m.clc 
getvariable.pl -c -0.022 "aS1(6,-3)" m.clc 
getvariable.pl -c 0.0053 "aS1(6,-2)" m.clc 
getvariable.pl -c 0.012 "aS1(6,-1)" m.clc 
getvariable.pl -c 0.007 "aS1(6,0)" m.clc 
getvariable.pl -c 0.028 "aS1(6,1)" m.clc 
getvariable.pl -c -0.0054 "aS1(6,2)" m.clc 
getvariable.pl -c -0.018 "aS1(6,3)" m.clc 
getvariable.pl -c -0.0047 "aS1(6,4)" m.clc 
getvariable.pl -c 0.0035 "aS1(6,5)" m.clc 
getvariable.pl -c 0.0022 "aS1(6,6)" m.clc 
getvariable.pl -c -0.12 "aS2(0,0)" m.clc 
getvariable.pl -c 0.000000 "aS2(1,-1)" m.clc 
getvariable.pl -c 0.000000 "aS2(1,0)" m.clc 
getvariable.pl -c 0.000000 "aS2(1,1)" m.clc 
getvariable.pl -c 0.067 "aS2(2,-2)" m.clc 
getvariable.pl -c 0.10 "aS2(2,-1)" m.clc 
getvariable.pl -c -0.0096 "aS2(2,0)" m.clc 
getvariable.pl -c 0.014 "aS2(2,1)" m.clc 
getvariable.pl -c -0.028 "aS2(2,2)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,-3)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,-2)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,-1)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,0)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,1)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,2)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,3)" m.clc 
getvariable.pl -c 0.0021 "aS2(4,-4)" m.clc 
getvariable.pl -c 0.017 "aS2(4,-3)" m.clc 
getvariable.pl -c 0.031 "aS2(4,-2)" m.clc 
getvariable.pl -c 0.011 "aS2(4,-1)" m.clc 
getvariable.pl -c -0.014 "aS2(4,0)" m.clc 
getvariable.pl -c -0.004 "aS2(4,1)" m.clc 
getvariable.pl -c -0.012 "aS2(4,2)" m.clc 
getvariable.pl -c -0.019 "aS2(4,3)" m.clc 
getvariable.pl -c -0.0087 "aS2(4,4)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,-5)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,-4)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,-3)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,-2)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,-1)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,0)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,1)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,2)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,3)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,4)" m.clc 
getvariable.pl -c 0.000000 "aS2(5,5)" m.clc 
getvariable.pl -c 0.0021 "aS2(6,-6)" m.clc 
getvariable.pl -c 0.0033 "aS2(6,-5)" m.clc 
getvariable.pl -c -0.0046 "aS2(6,-4)" m.clc 
getvariable.pl -c -0.017 "aS2(6,-3)" m.clc 
getvariable.pl -c -0.0086 "aS2(6,-2)" m.clc 
getvariable.pl -c 0.0092 "aS2(6,-1)" m.clc 
getvariable.pl -c 0.0037 "aS2(6,0)" m.clc 
getvariable.pl -c 0.012 "aS2(6,1)" m.clc 
getvariable.pl -c 0.012 "aS2(6,2)" m.clc 
getvariable.pl -c 0.026 "aS2(6,3)" m.clc 
getvariable.pl -c 0.026 "aS2(6,4)" m.clc 
getvariable.pl -c 0.012 "aS2(6,5)" m.clc 
getvariable.pl -c 0.0023 "aS2(6,6)" m.clc 
getvariable.pl -c -0.36 "aS3(0,0)" m.clc 
getvariable.pl -c 0.000000 "aS3(1,-1)" m.clc 
getvariable.pl -c 0.000000 "aS3(1,0)" m.clc 
getvariable.pl -c 0.000000 "aS3(1,1)" m.clc 
getvariable.pl -c 0.014 "aS3(2,-2)" m.clc 
getvariable.pl -c 0.052 "aS3(2,-1)" m.clc 
getvariable.pl -c 0.13 "aS3(2,0)" m.clc 
getvariable.pl -c 0.103 "aS3(2,1)" m.clc 
getvariable.pl -c 0.01 "aS3(2,2)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,-3)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,-2)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,-1)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,0)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,1)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,2)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,3)" m.clc 
getvariable.pl -c 0.0018 "aS3(4,-4)" m.clc 
getvariable.pl -c 0.013 "aS3(4,-3)" m.clc 
getvariable.pl -c 0.032 "aS3(4,-2)" m.clc 
getvariable.pl -c 0.024 "aS3(4,-1)" m.clc 
getvariable.pl -c 0.0003 "aS3(4,0)" m.clc 
getvariable.pl -c 0.05 "aS3(4,1)" m.clc 
getvariable.pl -c 0.02 "aS3(4,2)" m.clc 
getvariable.pl -c 0.002 "aS3(4,3)" m.clc 
getvariable.pl -c -0.0005 "aS3(4,4)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,-5)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,-4)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,-3)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,-2)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,-1)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,0)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,1)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,2)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,3)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,4)" m.clc 
getvariable.pl -c 0.000000 "aS3(5,5)" m.clc 
getvariable.pl -c -0.00002 "aS3(6,-6)" m.clc 
getvariable.pl -c -0.0015 "aS3(6,-5)" m.clc 
getvariable.pl -c -0.011 "aS3(6,-4)" m.clc 
getvariable.pl -c -0.029 "aS3(6,-3)" m.clc 
getvariable.pl -c -0.032 "aS3(6,-2)" m.clc 
getvariable.pl -c -0.0051 "aS3(6,-1)" m.clc 
getvariable.pl -c 0.025 "aS3(6,0)" m.clc 
getvariable.pl -c -0.010 "aS3(6,1)" m.clc 
getvariable.pl -c -0.024 "aS3(6,2)" m.clc 
getvariable.pl -c -0.0053 "aS3(6,3)" m.clc 
getvariable.pl -c 0.0031 "aS3(6,4)" m.clc 
getvariable.pl -c 0.0014 "aS3(6,5)" m.clc 
getvariable.pl -c 0.00006 "aS3(6,6)" m.clc

densplt s Pr3p_xyz_icf1ion.sipf 2 20 10 30 > m.clc

getvariable.pl -c -0.25   "aS1(0,0)" m.clc 
getvariable.pl -c 0.000000   "aS1(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aS1(1,0)" m.clc 
getvariable.pl -c 0.000000   "aS1(1,1)" m.clc 
getvariable.pl -c 0.0308   "aS1(2,-2)" m.clc 
getvariable.pl -c 0.020   "aS1(2,-1)" m.clc 
getvariable.pl -c 0.0053   "aS1(2,0)" m.clc 
getvariable.pl -c 0.092   "aS1(2,1)" m.clc 
getvariable.pl -c 0.045   "aS1(2,2)" m.clc 
getvariable.pl -c 0.000000   "aS1(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aS1(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aS1(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aS1(3,0)" m.clc 
getvariable.pl -c 0.000000   "aS1(3,1)" m.clc 
getvariable.pl -c 0.000000   "aS1(3,2)" m.clc 
getvariable.pl -c 0.000000   "aS1(3,3)" m.clc 
getvariable.pl -c 0.011   "aS1(4,-4)" m.clc 
getvariable.pl -c 0.027   "aS1(4,-3)" m.clc 
getvariable.pl -c 0.021   "aS1(4,-2)" m.clc 
getvariable.pl -c -0.0056   "aS1(4,-1)" m.clc 
getvariable.pl -c -0.031   "aS1(4,0)" m.clc 
getvariable.pl -c 0.0042   "aS1(4,1)" m.clc 
getvariable.pl -c 0.037   "aS1(4,2)" m.clc 
getvariable.pl -c 0.018   "aS1(4,3)" m.clc 
getvariable.pl -c 0.0015   "aS1(4,4)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,-5)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,-4)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,-3)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,-2)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,-1)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,0)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,1)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,2)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,3)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,4)" m.clc 
getvariable.pl -c 0.000000   "aS1(5,5)" m.clc 
getvariable.pl -c -0.0017   "aS1(6,-6)" m.clc 
getvariable.pl -c -0.0086   "aS1(6,-5)" m.clc 
getvariable.pl -c -0.0188   "aS1(6,-4)" m.clc 
getvariable.pl -c -0.016   "aS1(6,-3)" m.clc 
getvariable.pl -c 0.0033   "aS1(6,-2)" m.clc 
getvariable.pl -c 0.0090   "aS1(6,-1)" m.clc 
getvariable.pl -c 0.0059   "aS1(6,0)" m.clc 
getvariable.pl -c 0.020   "aS1(6,1)" m.clc 
getvariable.pl -c -0.0043   "aS1(6,2)" m.clc 
getvariable.pl -c -0.013   "aS1(6,3)" m.clc 
getvariable.pl -c -0.0032   "aS1(6,4)" m.clc 
getvariable.pl -c 0.0027   "aS1(6,5)" m.clc 
getvariable.pl -c 0.0016   "aS1(6,6)" m.clc 

getvariable.pl -c -0.124   "aS2(0,0)" m.clc 
getvariable.pl -c 0.000000   "aS2(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aS2(1,0)" m.clc 
getvariable.pl -c 0.000000   "aS2(1,1)" m.clc 
getvariable.pl -c 0.041   "aS2(2,-2)" m.clc 
getvariable.pl -c 0.062   "aS2(2,-1)" m.clc 
getvariable.pl -c 0.0027   "aS2(2,0)" m.clc 
getvariable.pl -c 0.020   "aS2(2,1)" m.clc 
getvariable.pl -c -0.012   "aS2(2,2)" m.clc 
getvariable.pl -c 0.000000   "aS2(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aS2(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aS2(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aS2(3,0)" m.clc 
getvariable.pl -c 0.000000   "aS2(3,1)" m.clc 
getvariable.pl -c 0.000000   "aS2(3,2)" m.clc 
getvariable.pl -c 0.000000   "aS2(3,3)" m.clc 
getvariable.pl -c 0.0024   "aS2(4,-4)" m.clc 
getvariable.pl -c 0.020   "aS2(4,-3)" m.clc 
getvariable.pl -c 0.036   "aS2(4,-2)" m.clc 
getvariable.pl -c 0.012   "aS2(4,-1)" m.clc 
getvariable.pl -c -0.016   "aS2(4,0)" m.clc 
getvariable.pl -c -0.0056   "aS2(4,1)" m.clc 
getvariable.pl -c -0.015   "aS2(4,2)" m.clc 
getvariable.pl -c -0.023   "aS2(4,3)" m.clc 
getvariable.pl -c -0.010   "aS2(4,4)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,-5)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,-4)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,-3)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,-2)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,-1)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,0)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,1)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,2)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,3)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,4)" m.clc 
getvariable.pl -c 0.000000   "aS2(5,5)" m.clc 
getvariable.pl -c 0.0015   "aS2(6,-6)" m.clc 
getvariable.pl -c 0.0023   "aS2(6,-5)" m.clc 
getvariable.pl -c -0.0035   "aS2(6,-4)" m.clc 
getvariable.pl -c -0.013   "aS2(6,-3)" m.clc 
getvariable.pl -c -0.0065   "aS2(6,-2)" m.clc 
getvariable.pl -c 0.0067   "aS2(6,-1)" m.clc 
getvariable.pl -c 0.0029   "aS2(6,0)" m.clc 
getvariable.pl -c 0.0091   "aS2(6,1)" m.clc 
getvariable.pl -c 0.0087   "aS2(6,2)" m.clc 
getvariable.pl -c 0.019   "aS2(6,3)" m.clc 
getvariable.pl -c 0.019   "aS2(6,4)" m.clc 
getvariable.pl -c 0.0084   "aS2(6,5)" m.clc 
getvariable.pl -c 0.0017   "aS2(6,6)" m.clc 

getvariable.pl -c -0.37   "aS3(0,0)" m.clc 
getvariable.pl -c 0.000000   "aS3(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aS3(1,0)" m.clc 
getvariable.pl -c 0.000000   "aS3(1,1)" m.clc 
getvariable.pl -c 0.020   "aS3(2,-2)" m.clc 
getvariable.pl -c 0.048   "aS3(2,-1)" m.clc 
getvariable.pl -c 0.098   "aS3(2,0)" m.clc 
getvariable.pl -c 0.095   "aS3(2,1)" m.clc 
getvariable.pl -c 0.015   "aS3(2,2)" m.clc 
getvariable.pl -c 0.000000   "aS3(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aS3(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aS3(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aS3(3,0)" m.clc 
getvariable.pl -c 0.000000   "aS3(3,1)" m.clc 
getvariable.pl -c 0.000000   "aS3(3,2)" m.clc 
getvariable.pl -c 0.000000   "aS3(3,3)" m.clc 
getvariable.pl -c 0.0017   "aS3(4,-4)" m.clc 
getvariable.pl -c 0.014   "aS3(4,-3)" m.clc 
getvariable.pl -c 0.035   "aS3(4,-2)" m.clc 
getvariable.pl -c 0.027   "aS3(4,-1)" m.clc 
getvariable.pl -c 0.001   "aS3(4,0)" m.clc 
getvariable.pl -c 0.055   "aS3(4,1)" m.clc 
getvariable.pl -c 0.026   "aS3(4,2)" m.clc 
getvariable.pl -c 0.002   "aS3(4,3)" m.clc 
getvariable.pl -c -0.0005   "aS3(4,4)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,-5)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,-4)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,-3)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,-2)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,-1)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,0)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,1)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,2)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,3)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,4)" m.clc 
getvariable.pl -c 0.000000   "aS3(5,5)" m.clc 
getvariable.pl -c -0.000038   "aS3(6,-6)" m.clc 
getvariable.pl -c -0.0013   "aS3(6,-5)" m.clc 
getvariable.pl -c -0.0084   "aS3(6,-4)" m.clc 
getvariable.pl -c -0.0221   "aS3(6,-3)" m.clc 
getvariable.pl -c -0.024   "aS3(6,-2)" m.clc 
getvariable.pl -c -0.0035   "aS3(6,-1)" m.clc 
getvariable.pl -c 0.0188   "aS3(6,0)" m.clc 
getvariable.pl -c -0.0070   "aS3(6,1)" m.clc 
getvariable.pl -c -0.018   "aS3(6,2)" m.clc 
getvariable.pl -c -0.0040   "aS3(6,3)" m.clc 
getvariable.pl -c 0.0025   "aS3(6,4)" m.clc 
getvariable.pl -c 0.0012   "aS3(6,5)" m.clc 
getvariable.pl -c 0.0001   "aS3(6,6)" m.clc 

densplt c Pr3p_ic1ion_trunc.sipf 2 10 20 30 > m.clc

getvariable.pl -c 0.564  "a(0,0)" m.clc 
getvariable.pl -c -0.047   "a(2,-2)" m.clc 
getvariable.pl -c -0.074   "a(2,-1)" m.clc 
getvariable.pl -c -0.087   "a(2,0)" m.clc 
getvariable.pl -c -0.133   "a(2,1)" m.clc 
getvariable.pl -c -0.029   "a(2,2)" m.clc 
getvariable.pl -c -0.006   "a(4,-4)" m.clc 
getvariable.pl -c -0.025   "a(4,-3)" m.clc 
getvariable.pl -c -0.042   "a(4,-2)" m.clc 
getvariable.pl -c -0.020   "a(4,-1)" m.clc 
getvariable.pl -c 0.015   "a(4,0)" m.clc 
getvariable.pl -c -0.036412   "a(4,1)" m.clc 
getvariable.pl -c -0.026082   "a(4,2)" m.clc 
getvariable.pl -c -0.001080   "a(4,3)" m.clc 
getvariable.pl -c 0.002917   "a(4,4)" m.clc 
getvariable.pl -c 0.000064   "a(6,-6)" m.clc 
getvariable.pl -c 0.001616   "a(6,-5)" m.clc 
getvariable.pl -c 0.006933   "a(6,-4)" m.clc 
getvariable.pl -c 0.012890   "a(6,-3)" m.clc 
getvariable.pl -c 0.008912   "a(6,-2)" m.clc 
getvariable.pl -c -0.001560   "a(6,-1)" m.clc 
getvariable.pl -c -0.008213   "a(6,0)" m.clc 
getvariable.pl -c -0.002823   "a(6,1)" m.clc 
getvariable.pl -c 0.005568   "a(6,2)" m.clc 
getvariable.pl -c 0.000684   "a(6,3)" m.clc 
getvariable.pl -c -0.003370   "a(6,4)" m.clc 
getvariable.pl -c -0.002297   "a(6,5)" m.clc 
getvariable.pl -c -0.000592   "a(6,6)" m.clc 


densplt s Pr3p_ic1ion_trunc.sipf 2 10 20 30 > m.clc

getvariable.pl -c -0.234558  "aS1(0,0)" m.clc 
getvariable.pl -c 0.000000  "aS1(1,-1)" m.clc 
getvariable.pl -c 0.000000  "aS1(1,0)" m.clc 
getvariable.pl -c 0.000000  "aS1(1,1)" m.clc 
getvariable.pl -c 0.043282  "aS1(2,-2)" m.clc 
getvariable.pl -c 0.014483  "aS1(2,-1)" m.clc 
getvariable.pl -c -0.017849  "aS1(2,0)" m.clc 
getvariable.pl -c 0.120119  "aS1(2,1)" m.clc 
getvariable.pl -c 0.065969  "aS1(2,2)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,-3)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,-2)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,-1)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,0)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,1)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,2)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,3)" m.clc 
getvariable.pl -c 0.009419  "aS1(4,-4)" m.clc 
getvariable.pl -c 0.025135  "aS1(4,-3)" m.clc 
getvariable.pl -c 0.020879  "aS1(4,-2)" m.clc 
getvariable.pl -c -0.004203  "aS1(4,-1)" m.clc 
getvariable.pl -c -0.026569  "aS1(4,0)" m.clc 
getvariable.pl -c 0.006002  "aS1(4,1)" m.clc 
getvariable.pl -c 0.031211  "aS1(4,2)" m.clc 
getvariable.pl -c 0.012843  "aS1(4,3)" m.clc 
getvariable.pl -c -0.000199  "aS1(4,4)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,-5)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,-4)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,-3)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,-2)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,-1)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,0)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,1)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,2)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,3)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,4)" m.clc 
getvariable.pl -c 0.000000  "aS1(5,5)" m.clc 
getvariable.pl -c -0.001630  "aS1(6,-6)" m.clc 
getvariable.pl -c -0.010422  "aS1(6,-5)" m.clc 
getvariable.pl -c -0.025401  "aS1(6,-4)" m.clc 
getvariable.pl -c -0.023560  "aS1(6,-3)" m.clc 
getvariable.pl -c 0.004517  "aS1(6,-2)" m.clc 
getvariable.pl -c 0.013338  "aS1(6,-1)" m.clc 
getvariable.pl -c 0.008033  "aS1(6,0)" m.clc 
getvariable.pl -c 0.026382  "aS1(6,1)" m.clc 
getvariable.pl -c -0.007283  "aS1(6,2)" m.clc 
getvariable.pl -c -0.016585  "aS1(6,3)" m.clc 
getvariable.pl -c -0.000978  "aS1(6,4)" m.clc 
getvariable.pl -c 0.005463  "aS1(6,5)" m.clc 
getvariable.pl -c 0.002497  "aS1(6,6)" m.clc 

getvariable.pl -c -0.128455  "aS2(0,0)" m.clc 
getvariable.pl -c 0.000000  "aS2(1,-1)" m.clc 
getvariable.pl -c 0.000000  "aS2(1,0)" m.clc 
getvariable.pl -c 0.000000  "aS2(1,1)" m.clc 
getvariable.pl -c 0.065252  "aS2(2,-2)" m.clc 
getvariable.pl -c 0.102053  "aS2(2,-1)" m.clc 
getvariable.pl -c -0.010601  "aS2(2,0)" m.clc 
getvariable.pl -c 0.013718  "aS2(2,1)" m.clc 
getvariable.pl -c -0.031017  "aS2(2,2)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,-3)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,-2)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,-1)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,0)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,1)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,2)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,3)" m.clc 
getvariable.pl -c 0.000866  "aS2(4,-4)" m.clc 
getvariable.pl -c 0.015008  "aS2(4,-3)" m.clc 
getvariable.pl -c 0.030691  "aS2(4,-2)" m.clc 
getvariable.pl -c 0.011016  "aS2(4,-1)" m.clc 
getvariable.pl -c -0.014812  "aS2(4,0)" m.clc 
getvariable.pl -c -0.004391  "aS2(4,1)" m.clc 
getvariable.pl -c -0.014075  "aS2(4,2)" m.clc 
getvariable.pl -c -0.020596  "aS2(4,3)" m.clc 
getvariable.pl -c -0.008664  "aS2(4,4)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,-5)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,-4)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,-3)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,-2)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,-1)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,0)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,1)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,2)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,3)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,4)" m.clc 
getvariable.pl -c 0.000000  "aS2(5,5)" m.clc 
getvariable.pl -c 0.002455  "aS2(6,-6)" m.clc 
getvariable.pl -c 0.005253  "aS2(6,-5)" m.clc 
getvariable.pl -c -0.000834  "aS2(6,-4)" m.clc 
getvariable.pl -c -0.014495  "aS2(6,-3)" m.clc 
getvariable.pl -c -0.007866  "aS2(6,-2)" m.clc 
getvariable.pl -c 0.010133  "aS2(6,-1)" m.clc 
getvariable.pl -c 0.004500  "aS2(6,0)" m.clc 
getvariable.pl -c 0.013333  "aS2(6,1)" m.clc 
getvariable.pl -c 0.013052  "aS2(6,2)" m.clc 
getvariable.pl -c 0.027487  "aS2(6,3)" m.clc 
getvariable.pl -c 0.025446  "aS2(6,4)" m.clc 
getvariable.pl -c 0.010389  "aS2(6,5)" m.clc 
getvariable.pl -c 0.001645  "aS2(6,6)" m.clc 

getvariable.pl -c -0.363014  "aS3(0,0)" m.clc 
getvariable.pl -c 0.000000  "aS3(1,-1)" m.clc 
getvariable.pl -c 0.000000  "aS3(1,0)" m.clc 
getvariable.pl -c 0.000000  "aS3(1,1)" m.clc 
getvariable.pl -c 0.014202  "aS3(2,-2)" m.clc 
getvariable.pl -c 0.056218  "aS3(2,-1)" m.clc 
getvariable.pl -c 0.134939  "aS3(2,0)" m.clc 
getvariable.pl -c 0.099821  "aS3(2,1)" m.clc 
getvariable.pl -c 0.008650  "aS3(2,2)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,-3)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,-2)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,-1)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,0)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,1)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,2)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,3)" m.clc 
getvariable.pl -c 0.001577  "aS3(4,-4)" m.clc 
getvariable.pl -c 0.012730  "aS3(4,-3)" m.clc 
getvariable.pl -c 0.033637  "aS3(4,-2)" m.clc 
getvariable.pl -c 0.026663  "aS3(4,-1)" m.clc 
getvariable.pl -c 0.001260  "aS3(4,0)" m.clc 
getvariable.pl -c 0.047573  "aS3(4,1)" m.clc 
getvariable.pl -c 0.020609  "aS3(4,2)" m.clc 
getvariable.pl -c 0.000511  "aS3(4,3)" m.clc 
getvariable.pl -c -0.000796  "aS3(4,4)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,-5)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,-4)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,-3)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,-2)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,-1)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,0)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,1)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,2)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,3)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,4)" m.clc 
getvariable.pl -c 0.000000  "aS3(5,5)" m.clc 
getvariable.pl -c -0.000008  "aS3(6,-6)" m.clc 
getvariable.pl -c -0.001108  "aS3(6,-5)" m.clc 
getvariable.pl -c -0.009505  "aS3(6,-4)" m.clc 
getvariable.pl -c -0.028841  "aS3(6,-3)" m.clc 
getvariable.pl -c -0.034375  "aS3(6,-2)" m.clc 
getvariable.pl -c -0.006278  "aS3(6,-1)" m.clc 
getvariable.pl -c 0.024577  "aS3(6,0)" m.clc 
getvariable.pl -c -0.011144  "aS3(6,1)" m.clc 
getvariable.pl -c -0.020954  "aS3(6,2)" m.clc 
getvariable.pl -c -0.001065  "aS3(6,3)" m.clc 
getvariable.pl -c 0.004860  "aS3(6,4)" m.clc 
getvariable.pl -c 0.001649  "aS3(6,5)" m.clc 
getvariable.pl -c 0.000060  "aS3(6,6)" m.clc 

densplt o Pr3p_xyz_ic1ion.sipf 2 20 10 30 > m.clc

getvariable.pl -c 1.4472   "aL1(0,0)" m.clc 
getvariable.pl -c 0.000000   "aL1(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aL1(1,0)" m.clc 
getvariable.pl -c 0.000000   "aL1(1,1)" m.clc 
getvariable.pl -c -0.2791   "aL1(2,-2)" m.clc 
getvariable.pl -c 0.0059   "aL1(2,-1)" m.clc 
getvariable.pl -c 0.3342   "aL1(2,0)" m.clc 
getvariable.pl -c -0.8372   "aL1(2,1)" m.clc 
getvariable.pl -c -0.5631   "aL1(2,2)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,0)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,1)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,2)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,3)" m.clc 
getvariable.pl -c 0.002923   "aL1(4,-4)" m.clc 
getvariable.pl -c -0.01172   "aL1(4,-3)" m.clc 
getvariable.pl -c -0.05181   "aL1(4,-2)" m.clc 
getvariable.pl -c -0.04576   "aL1(4,-1)" m.clc 
getvariable.pl -c -0.01239   "aL1(4,0)" m.clc 
getvariable.pl -c -0.06908   "aL1(4,1)" m.clc 
getvariable.pl -c -0.00800   "aL1(4,2)" m.clc 
getvariable.pl -c 0.0171   "aL1(4,3)" m.clc 
getvariable.pl -c 0.0061   "aL1(4,4)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,-5)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,-4)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,-3)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,-2)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,-1)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,0)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,1)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,2)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,3)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,4)" m.clc 
getvariable.pl -c 0.000000   "aL1(5,5)" m.clc 
getvariable.pl -c 0.004209   "aL1(6,-6)" m.clc 
getvariable.pl -c 0.021339   "aL1(6,-5)" m.clc 
getvariable.pl -c 0.046051   "aL1(6,-4)" m.clc 
getvariable.pl -c 0.038209   "aL1(6,-3)" m.clc 
getvariable.pl -c -0.01087   "aL1(6,-2)" m.clc 
getvariable.pl -c -0.02276   "aL1(6,-1)" m.clc 
getvariable.pl -c -0.01271   "aL1(6,0)" m.clc 
getvariable.pl -c -0.0511   "aL1(6,1)" m.clc 
getvariable.pl -c 0.00918   "aL1(6,2)" m.clc 
getvariable.pl -c 0.03385   "aL1(6,3)" m.clc 
getvariable.pl -c 0.00895   "aL1(6,4)" m.clc 
getvariable.pl -c -0.0062   "aL1(6,5)" m.clc 
getvariable.pl -c -0.0039   "aL1(6,6)" m.clc 

getvariable.pl -c 0.7236   "aL2(0,0)" m.clc 
getvariable.pl -c 0.000000   "aL2(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aL2(1,0)" m.clc 
getvariable.pl -c 0.000000   "aL2(1,1)" m.clc 
getvariable.pl -c -0.5640   "aL2(2,-2)" m.clc 
getvariable.pl -c -0.8461   "aL2(2,-1)" m.clc 
getvariable.pl -c 0.1671   "aL2(2,0)" m.clc 
getvariable.pl -c 0.0059   "aL2(2,1)" m.clc 
getvariable.pl -c 0.2845   "aL2(2,2)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,0)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,1)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,2)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,3)" m.clc 
getvariable.pl -c -0.0032   "aL2(4,-4)" m.clc 
getvariable.pl -c 0.00262   "aL2(4,-3)" m.clc 
getvariable.pl -c 0.01112   "aL2(4,-2)" m.clc 
getvariable.pl -c -0.00044   "aL2(4,-1)" m.clc 
getvariable.pl -c -0.0062   "aL2(4,0)" m.clc 
getvariable.pl -c -0.0458   "aL2(4,1)" m.clc 
getvariable.pl -c -0.0534   "aL2(4,2)" m.clc 
getvariable.pl -c -0.0381   "aL2(4,3)" m.clc 
getvariable.pl -c -0.0129   "aL2(4,4)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,-5)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,-4)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,-3)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,-2)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,-1)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,0)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,1)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,2)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,3)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,4)" m.clc 
getvariable.pl -c 0.000000   "aL2(5,5)" m.clc 
getvariable.pl -c -0.0039   "aL2(6,-6)" m.clc 
getvariable.pl -c -0.0062   "aL2(6,-5)" m.clc 
getvariable.pl -c 0.00799   "aL2(6,-4)" m.clc 
getvariable.pl -c 0.03094   "aL2(6,-3)" m.clc 
getvariable.pl -c 0.01537   "aL2(6,-2)" m.clc 
getvariable.pl -c -0.0169   "aL2(6,-1)" m.clc 
getvariable.pl -c -0.0064   "aL2(6,0)" m.clc 
getvariable.pl -c -0.02276   "aL2(6,1)" m.clc 
getvariable.pl -c -0.02314   "aL2(6,2)" m.clc 
getvariable.pl -c -0.04818   "aL2(6,3)" m.clc 
getvariable.pl -c -0.04709   "aL2(6,4)" m.clc 
getvariable.pl -c -0.0213   "aL2(6,5)" m.clc 
getvariable.pl -c -0.0042   "aL2(6,6)" m.clc 

getvariable.pl -c 2.17   "aL3(0,0)" m.clc 
getvariable.pl -c 0.000000   "aL3(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aL3(1,0)" m.clc 
getvariable.pl -c 0.000000   "aL3(1,1)" m.clc 
getvariable.pl -c 0.006   "aL3(2,-2)" m.clc 
getvariable.pl -c -0.27   "aL3(2,-1)" m.clc 
getvariable.pl -c -0.97   "aL3(2,0)" m.clc 
getvariable.pl -c -0.55   "aL3(2,1)" m.clc 
getvariable.pl -c 0.004   "aL3(2,2)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,0)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,1)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,2)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,3)" m.clc 
getvariable.pl -c -0.018   "aL3(4,-4)" m.clc 
getvariable.pl -c -0.058   "aL3(4,-3)" m.clc 
getvariable.pl -c -0.072   "aL3(4,-2)" m.clc 
getvariable.pl -c -0.016   "aL3(4,-1)" m.clc 
getvariable.pl -c 0.052   "aL3(4,0)" m.clc 
getvariable.pl -c -0.033   "aL3(4,1)" m.clc 
getvariable.pl -c -0.054   "aL3(4,2)" m.clc 
getvariable.pl -c -0.011   "aL3(4,3)" m.clc 
getvariable.pl -c 0.0051   "aL3(4,4)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,-5)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,-4)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,-3)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,-2)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,-1)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,0)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,1)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,2)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,3)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,4)" m.clc 
getvariable.pl -c 0.000000   "aL3(5,5)" m.clc 
getvariable.pl -c 0.000000   "aL3(6,-6)" m.clc 
getvariable.pl -c 0.0024   "aL3(6,-5)" m.clc 
getvariable.pl -c 0.018   "aL3(6,-4)" m.clc 
getvariable.pl -c 0.051   "aL3(6,-3)" m.clc 
getvariable.pl -c 0.058   "aL3(6,-2)" m.clc 
getvariable.pl -c 0.01   "aL3(6,-1)" m.clc 
getvariable.pl -c -0.05   "aL3(6,0)" m.clc 
getvariable.pl -c 0.019   "aL3(6,1)" m.clc 
getvariable.pl -c 0.043   "aL3(6,2)" m.clc 
getvariable.pl -c 0.009   "aL3(6,3)" m.clc 
getvariable.pl -c -0.005   "aL3(6,4)" m.clc 
getvariable.pl -c -0.002   "aL3(6,5)" m.clc 
getvariable.pl -c 0.000000   "aL3(6,6)" m.clc 


singleion -sx -r Pr3p_xyz_ic1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.241 2 9 2 0 m.clc
getvalue.pl -c 0.040    2 13 2 0 m.clc
getvalue.pl -c 0.0135   2 14 2 0 m.clc
getvalue.pl -c -0.01926 2 15 2 0 m.clc
getvalue.pl -c 0.12094  2 16 2 0 m.clc
getvalue.pl -c 0.0694 2 17 2 0 m.clc
getvalue.pl -c 0.0154   0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.00267  0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0002   0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

singleion -sx -r Pr3p_ic1ion_trunc.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.24097 2 9 2 0 m.clc
getvalue.pl -c 0.0403  2 13 2 0 m.clc
getvalue.pl -c 0.01351   2 14 2 0 m.clc
getvalue.pl -c -0.01921 2 15 2 0 m.clc
getvalue.pl -c 0.1209  2 16 2 0 m.clc
getvalue.pl -c 0.0693 2 17 2 0 m.clc
getvalue.pl -c 0.0154   0 11 4 0 results/Pr3p_ic1ion_trunc.sipf.trs
getvalue.pl -c 0.00267  0 15 4 0 results/Pr3p_ic1ion_trunc.sipf.trs
getvalue.pl -c 0.0002   0 16 4 0 results/Pr3p_ic1ion_trunc.sipf.trs


singleion -sx -r Pr3p_xyz_icf1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.24891 2 9 2 0 m.clc
getvalue.pl -c 0.031 2 13 2 0 m.clc
getvalue.pl -c 0.02023  2 14 2 0 m.clc
getvalue.pl -c 0.005308  2 15 2 0 m.clc 
getvalue.pl -c 0.0924  2 16 2 0 m.clc
getvalue.pl -c 0.0447  2 17 2 0 m.clc
getvalue.pl -c 0.0162     0 11 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.00155    0 15 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.00049   0 16 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
  
 
singleion -lx -r Pr3p_xyz_ic1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c 1.44     2 9 2 0 m.clc
getvalue.pl -c -0.28    2 13 2 0 m.clc
getvalue.pl -c 0.006  2 14 2 0 m.clc
getvalue.pl -c 0.34 2 15 2 0 m.clc
getvalue.pl -c -0.84  2 16 2 0 m.clc
getvalue.pl -c -0.56 2 17 2 0 m.clc
getvalue.pl -c 0.64    0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.12  0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 6e-05  0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

singleion -lx -r Pr3p_xyz_icf1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c 1.4552  2 9 2 0 m.clc
getvalue.pl -c -0.282 2 13 2 0 m.clc
getvalue.pl -c 0  2 14 2 0 m.clc
getvalue.pl -c 0.3254  2 15 2 0 m.clc 
getvalue.pl -c -0.8454   2 16 2 0 m.clc
getvalue.pl -c -0.5636  2 17 2 0 m.clc
getvalue.pl -c 0.6430      0 11 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.1254     0 15 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0             0 16 4 0 results/Pr3p_xyz_icf1ion.sipf.trs

singleion -sy -r Pr3p_xyz_ic1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.12    2 9 2 0 m.clc
getvalue.pl -c 0.067    2 13 2 0 m.clc
getvalue.pl -c 0.10     2 14 2 0 m.clc
getvalue.pl -c -0.01  2 15 2 0 m.clc 
getvalue.pl -c 0.01    2 16 2 0 m.clc
getvalue.pl -c -0.028   2 17 2 0 m.clc
    
getvalue.pl -c 0.01998      0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.001468    0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.000677   0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

singleion -sy -r Pr3p_xyz_icf1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.124457   2 9 2 0 m.clc
getvalue.pl -c 0.04137 2 13 2 0 m.clc
getvalue.pl -c 0.06206   2 14 2 0 m.clc
getvalue.pl -c 0.002654  2 15 2 0 m.clc 
getvalue.pl -c 0.02023   2 16 2 0 m.clc
getvalue.pl -c -0.01226  2 17 2 0 m.clc
getvalue.pl -c 0.02106      0 11 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.000671     0 15 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.0005508    0 16 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
    

singleion -ly -r Pr3p_xyz_ic1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c  0.72    2 9 2 0 m.clc
getvalue.pl -c -0.56    2 13 2 0 m.clc
getvalue.pl -c -0.84     2 14 2 0 m.clc
getvalue.pl -c 0.17  2 15 2 0 m.clc 
getvalue.pl -c 0.006    2 16 2 0 m.clc
getvalue.pl -c 0.28   2 17 2 0 m.clc
        
getvalue.pl -c 0.828884      0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.097736      0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0500397     0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

singleion -ly -r Pr3p_xyz_icf1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c 0.727599    2 9 2 0 m.clc
getvalue.pl -c -0.563596 2 13 2 0 m.clc
getvalue.pl -c -0.845397   2 14 2 0 m.clc
getvalue.pl -c 0.162696  2 15 2 0 m.clc 
getvalue.pl -c 0   2 16 2 0 m.clc
getvalue.pl -c 0.281798  2 17 2 0 m.clc
getvalue.pl -c 0.835901       0 11 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.0964507     0 15 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.0482251     0 16 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
 
singleion -sz -r Pr3p_xyz_ic1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.36   2 9 2 0 m.clc
getvalue.pl -c 0.014        2 13 2 0 m.clc
getvalue.pl -c 0.052    2 14 2 0 m.clc
getvalue.pl -c 0.133 2 15 2 0 m.clc 
getvalue.pl -c 0.103    2 16 2 0 m.clc
getvalue.pl -c 0.010   2 17 2 0 m.clc
    
getvalue.pl -c 0.00768446       0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.000218828  0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.00400595   0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

singleion -sz -r Pr3p_xyz_icf1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.373372     2 9 2 0 m.clc
getvalue.pl -c 0.0202265   2 13 2 0 m.clc
getvalue.pl -c 0.0476551   2 14 2 0 m.clc
getvalue.pl -c 0.0979351  2 15 2 0 m.clc 
getvalue.pl -c 0.0953101 2 16 2 0 m.clc
getvalue.pl -c  0.0151699   2 17 2 0 m.clc
getvalue.pl -c 0.00809825       0 11 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.000506804   0 15 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0.00319181     0 16 4 0 results/Pr3p_xyz_icf1ion.sipf.trs

singleion -lz -r Pr3p_xyz_ic1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c 2.17     2 9 2 0 m.clc
getvalue.pl -c 0.0059      2 13 2 0 m.clc
getvalue.pl -c -0.274    2 14 2 0 m.clc
getvalue.pl -c -0.969  2 15 2 0 m.clc 
getvalue.pl -c -0.55    2 16 2 0 m.clc
getvalue.pl -c 0.0045   2 17 2 0 m.clc
    
getvalue.pl -c 0.318801     0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 5.88991e-05     0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.119468   0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

singleion -lz -r Pr3p_xyz_icf1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c 2.18281      2 9 2 0 m.clc
getvalue.pl -c 0   2 13 2 0 m.clc
getvalue.pl -c -0.281798     2 14 2 0 m.clc
getvalue.pl -c  -0.976181  2 15 2 0 m.clc 
getvalue.pl -c -0.563596 2 16 2 0 m.clc
getvalue.pl -c  0   2 17 2 0 m.clc
getvalue.pl -c 0.3215       0 11 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0   0 15 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.125385    0 16 4 0 results/Pr3p_xyz_icf1ion.sipf.trs

singleion -opmat 2 -r Pr3p_xyz_icf1ion.sipf 2 20 10 30 0 0 0

getvalue.pl -c +0.282843 0 2 34 0 results/Pr3p_xyz_icf1ion.sipf.opmat 

singleion -opmat 2 -r Pr3p_xyz_ic1ion.sipf 2 20 10 30 0 0 0

getvalue.pl -c +0.577350 0 2 92 0 results/Pr3p_xyz_ic1ion.sipf.opmat 



rm m.clc
cd ../../demo

# commands to write above checks for densplt using output (m.clc)
# cp m.clc dd ;  substitute "#\!" "getvariable.pl -c" dd ; swapcol 3 5 dd; substitute "=" " " dd ; 
# substitute 'aS' '"aS' dd ; substitute ')' ')" m.clc' dd
