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
getvalue.pl -c 3.37204      0 7 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.0910267  0 8 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0           0 9 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 2            0 10 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0158445    0 11 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.0205977   0 12 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.00792222  0 13 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   1.14772    0 14 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   1.49202    0 15 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.573855    0 16 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0275589    0 17 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0524248    0 18 4 0 results/Ce3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0341042    0 19 4 0 results/Ce3p_xyz_ic1ion.sipf.trs



densplt c Ce3p_xyz_ic1ion.sipf 2 20 10 30 > m.clc
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
getvalue.pl -c   6.6393       0 7 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c     0.137944   0 8 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    0           0 9 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    2           0 10 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    0.064627    0 11 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.0840141    0 12 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.0323135    0 13 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    2.30302     0 14 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    2.99391     0 15 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    1.15151     0 16 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.00270587   0 17 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c    0.00475756  0 18 4 0 results/Nd3p_xyz_ic1ion.sipf.trs
getvalue.pl -c   0.0024391    0 19 4 0 results/Nd3p_xyz_ic1ion.sipf.trs


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
getvalue.pl -c   0.00231595 0 17 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c   0.00446039 0 18 4 0 results/Nd3p_xyz_icf1ion.sipf.trs
getvalue.pl -c  0.00289497  0 19 4 0 results/Nd3p_xyz_icf1ion.sipf.trs


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

densplt c Nd3p_xyz_ic1ion.sipf 2 20 10 30 > m.clc

getvariable.pl -c 0.8463 "a(0,0)" m.clc 
getvariable.pl -c -0.01858 "a(2,-2)" m.clc 
getvariable.pl -c -0.02787 "a(2,-1)" m.clc 
getvariable.pl -c -0.03486 "a(2,0)" m.clc 
getvariable.pl -c -0.05573 "a(2,1)" m.clc 
getvariable.pl -c -0.01393 "a(2,2)" m.clc 
getvariable.pl -c -0.00522 "a(4,-4)" m.clc 
getvariable.pl -c -0.02029 "a(4,-3)" m.clc 
getvariable.pl -c -0.03220 "a(4,-2)" m.clc 
getvariable.pl -c -0.01464 "a(4,-1)" m.clc 
getvariable.pl -c 0.01312 "a(4,0)" m.clc 
getvariable.pl -c -0.02928 "a(4,1)" m.clc 
getvariable.pl -c -0.02415 "a(4,2)" m.clc 
getvariable.pl -c -0.0037 "a(4,3)" m.clc 
getvariable.pl -c 0.00152 "a(4,4)" m.clc 
getvariable.pl -c -0.00071 "a(6,-6)" m.clc 
getvariable.pl -c -0.00692 "a(6,-5)" m.clc 
getvariable.pl -c -0.02445 "a(6,-4)" m.clc 
getvariable.pl -c -0.04116 "a(6,-3)" m.clc 
getvariable.pl -c -0.02630 "a(6,-2)" m.clc 
getvariable.pl -c 0.00527 "a(6,-1)" m.clc 
getvariable.pl -c 0.02614 "a(6,0)" m.clc 
getvariable.pl -c 0.01055 "a(6,1)" m.clc 
getvariable.pl -c -0.01973 "a(6,2)" m.clc 
getvariable.pl -c -0.00748 "a(6,3)" m.clc 
getvariable.pl -c 0.00713 "a(6,4)" m.clc 
getvariable.pl -c 0.00641 "a(6,5)" m.clc 
getvariable.pl -c 0.0019 "a(6,6)" m.clc


densplt s Pr3p_xyz_ic1ion.sipf 2 20 10 30 > m.clc

getvariable.pl -c -0.240938 "aS1(0,0)" m.clc 
getvariable.pl -c 0.000000 "aS1(1,-1)" m.clc 
getvariable.pl -c 0.000000 "aS1(1,0)" m.clc 
getvariable.pl -c 0.000000 "aS1(1,1)" m.clc 
getvariable.pl -c 0.040311 "aS1(2,-2)" m.clc 
getvariable.pl -c 0.013504 "aS1(2,-1)" m.clc 
getvariable.pl -c -0.019260 "aS1(2,0)" m.clc 
getvariable.pl -c 0.120935 "aS1(2,1)" m.clc 
getvariable.pl -c 0.069370 "aS1(2,2)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,-3)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,-2)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,-1)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,0)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,1)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,2)" m.clc 
getvariable.pl -c 0.000000 "aS1(3,3)" m.clc 
getvariable.pl -c 0.009669 "aS1(4,-4)" m.clc 
getvariable.pl -c 0.024219 "aS1(4,-3)" m.clc 
getvariable.pl -c 0.019196 "aS1(4,-2)" m.clc 
getvariable.pl -c -0.004151 "aS1(4,-1)" m.clc 
getvariable.pl -c -0.027234 "aS1(4,0)" m.clc 
getvariable.pl -c 0.004707 "aS1(4,1)" m.clc 
getvariable.pl -c 0.032285 "aS1(4,2)" m.clc 
getvariable.pl -c 0.015579 "aS1(4,3)" m.clc 
getvariable.pl -c 0.001204 "aS1(4,4)" m.clc 
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
getvariable.pl -c -0.002290 "aS1(6,-6)" m.clc 
getvariable.pl -c -0.011678 "aS1(6,-5)" m.clc 
getvariable.pl -c -0.025401 "aS1(6,-4)" m.clc 
getvariable.pl -c -0.021510 "aS1(6,-3)" m.clc 
getvariable.pl -c 0.005325 "aS1(6,-2)" m.clc 
getvariable.pl -c 0.012414 "aS1(6,-1)" m.clc 
getvariable.pl -c 0.007413 "aS1(6,0)" m.clc 
getvariable.pl -c 0.027833 "aS1(6,1)" m.clc 
getvariable.pl -c -0.005375 "aS1(6,2)" m.clc 
getvariable.pl -c -0.018451 "aS1(6,3)" m.clc 
getvariable.pl -c -0.004687 "aS1(6,4)" m.clc 
getvariable.pl -c 0.003498 "aS1(6,5)" m.clc 
getvariable.pl -c 0.002148 "aS1(6,6)" m.clc 
getvariable.pl -c -0.120469 "aS2(0,0)" m.clc 
getvariable.pl -c 0.000000 "aS2(1,-1)" m.clc 
getvariable.pl -c 0.000000 "aS2(1,0)" m.clc 
getvariable.pl -c 0.000000 "aS2(1,1)" m.clc 
getvariable.pl -c 0.067119 "aS2(2,-2)" m.clc 
getvariable.pl -c 0.100679 "aS2(2,-1)" m.clc 
getvariable.pl -c -0.009630 "aS2(2,0)" m.clc 
getvariable.pl -c 0.013504 "aS2(2,1)" m.clc 
getvariable.pl -c -0.027933 "aS2(2,2)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,-3)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,-2)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,-1)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,0)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,1)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,2)" m.clc 
getvariable.pl -c 0.000000 "aS2(3,3)" m.clc 
getvariable.pl -c 0.002130 "aS2(4,-4)" m.clc 
getvariable.pl -c 0.017027 "aS2(4,-3)" m.clc 
getvariable.pl -c 0.031063 "aS2(4,-2)" m.clc 
getvariable.pl -c 0.010935 "aS2(4,-1)" m.clc 
getvariable.pl -c -0.013617 "aS2(4,0)" m.clc 
getvariable.pl -c -0.004151 "aS2(4,1)" m.clc 
getvariable.pl -c -0.012478 "aS2(4,2)" m.clc 
getvariable.pl -c -0.019256 "aS2(4,3)" m.clc 
getvariable.pl -c -0.008670 "aS2(4,4)" m.clc 
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
getvariable.pl -c 0.002101 "aS2(6,-6)" m.clc 
getvariable.pl -c 0.003292 "aS2(6,-5)" m.clc 
getvariable.pl -c -0.004572 "aS2(6,-4)" m.clc 
getvariable.pl -c -0.017153 "aS2(6,-3)" m.clc 
getvariable.pl -c -0.008580 "aS2(6,-2)" m.clc 
getvariable.pl -c 0.009212 "aS2(6,-1)" m.clc 
getvariable.pl -c 0.003706 "aS2(6,0)" m.clc 
getvariable.pl -c 0.012414 "aS2(6,1)" m.clc 
getvariable.pl -c 0.012303 "aS2(6,2)" m.clc 
getvariable.pl -c 0.025962 "aS2(6,3)" m.clc 
getvariable.pl -c 0.025525 "aS2(6,4)" m.clc 
getvariable.pl -c 0.011601 "aS2(6,5)" m.clc 
getvariable.pl -c 0.002295 "aS2(6,6)" m.clc 
getvariable.pl -c -0.361409 "aS3(0,0)" m.clc 
getvariable.pl -c 0.000000 "aS3(1,-1)" m.clc 
getvariable.pl -c 0.000000 "aS3(1,0)" m.clc 
getvariable.pl -c 0.000000 "aS3(1,1)" m.clc 
getvariable.pl -c 0.013504 "aS3(2,-2)" m.clc 
getvariable.pl -c 0.051565 "aS3(2,-1)" m.clc 
getvariable.pl -c 0.133796 "aS3(2,0)" m.clc 
getvariable.pl -c 0.103129 "aS3(2,1)" m.clc 
getvariable.pl -c 0.010128 "aS3(2,2)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,-3)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,-2)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,-1)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,0)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,1)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,2)" m.clc 
getvariable.pl -c 0.000000 "aS3(3,3)" m.clc 
getvariable.pl -c 0.001755 "aS3(4,-4)" m.clc 
getvariable.pl -c 0.012835 "aS3(4,-3)" m.clc 
getvariable.pl -c 0.031860 "aS3(4,-2)" m.clc 
getvariable.pl -c 0.024130 "aS3(4,-1)" m.clc 
getvariable.pl -c 0.000290 "aS3(4,0)" m.clc 
getvariable.pl -c 0.048259 "aS3(4,1)" m.clc 
getvariable.pl -c 0.023895 "aS3(4,2)" m.clc 
getvariable.pl -c 0.002334 "aS3(4,3)" m.clc 
getvariable.pl -c -0.000512 "aS3(4,4)" m.clc 
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
getvariable.pl -c -0.000022 "aS3(6,-6)" m.clc 
getvariable.pl -c -0.001530 "aS3(6,-5)" m.clc 
getvariable.pl -c -0.010601 "aS3(6,-4)" m.clc 
getvariable.pl -c -0.028861 "aS3(6,-3)" m.clc 
getvariable.pl -c -0.031951 "aS3(6,-2)" m.clc 
getvariable.pl -c -0.005080 "aS3(6,-1)" m.clc 
getvariable.pl -c 0.024890 "aS3(6,0)" m.clc 
getvariable.pl -c -0.010161 "aS3(6,1)" m.clc 
getvariable.pl -c -0.023963 "aS3(6,2)" m.clc 
getvariable.pl -c -0.005247 "aS3(6,3)" m.clc 
getvariable.pl -c 0.003092 "aS3(6,4)" m.clc 
getvariable.pl -c 0.001418 "aS3(6,5)" m.clc 
getvariable.pl -c 0.000059 "aS3(6,6)" m.clc



rm m.clc
cd ../../demo

