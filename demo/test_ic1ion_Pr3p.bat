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


singleion -XM 0 0 -muBT -nt 10000 -r Co2p_atom1_rotated_pm_z.sipf 2 0 0 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 12.08 2 20 2 0 m.clc




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
getvariable.pl -c 0.3 "a(0,0)" m.clc 
getvariable.pl -c -0.05 "a(2,-2)" m.clc 
getvariable.pl -c -0.07 "a(2,-1)" m.clc 
getvariable.pl -c -0.09 "a(2,0)" m.clc 
getvariable.pl -c -0.1 "a(2,1)" m.clc 
getvariable.pl -c -0.03 "a(2,2)" m.clc 
getvariable.pl -c 0.004 "a(4,-4)" m.clc 
getvariable.pl -c 0.02 "a(4,-3)" m.clc 
getvariable.pl -c 0.03 "a(4,-2)" m.clc 
getvariable.pl -c 0.01 "a(4,-1)" m.clc 
getvariable.pl -c -0.01 "a(4,0)" m.clc 
getvariable.pl -c 0.02 "a(4,1)" m.clc 
getvariable.pl -c 0.02 "a(4,2)" m.clc 
getvariable.pl -c 0.003 "a(4,3)" m.clc 
getvariable.pl -c -0.001 "a(4,4)" m.clc 
getvariable.pl -c -0.00002 "a(6,-6)" m.clc 
getvariable.pl -c -0.0002 "a(6,-5)" m.clc 
getvariable.pl -c -0.0007 "a(6,-4)" m.clc 
getvariable.pl -c -0.001 "a(6,-3)" m.clc 
getvariable.pl -c -0.0008 "a(6,-2)" m.clc 
getvariable.pl -c 0.0002 "a(6,-1)" m.clc 
getvariable.pl -c 0.0008 "a(6,0)" m.clc 
getvariable.pl -c 0.0003 "a(6,1)" m.clc 
getvariable.pl -c -0.0006 "a(6,2)" m.clc 
getvariable.pl -c -0.0002 "a(6,3)" m.clc 
getvariable.pl -c 0.0002 "a(6,4)" m.clc 
getvariable.pl -c 0.0002 "a(6,5)" m.clc 
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
getvalue.pl -c 0.8035 4 7 1.73205 0 results/ic1ion.mag 

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

getvariable.pl -c 0.9 "a(0,0)" m.clc 
getvariable.pl -c -0.02 "a(2,-2)" m.clc 
getvariable.pl -c -0.03 "a(2,-1)" m.clc 
getvariable.pl -c -0.04 "a(2,0)" m.clc 
getvariable.pl -c -0.06 "a(2,1)" m.clc 
getvariable.pl -c -0.01 "a(2,2)" m.clc 
getvariable.pl -c -0.005 "a(4,-4)" m.clc 
getvariable.pl -c -0.02 "a(4,-3)" m.clc 
getvariable.pl -c -0.03 "a(4,-2)" m.clc 
getvariable.pl -c -0.02 "a(4,-1)" m.clc 
getvariable.pl -c 0.01 "a(4,0)" m.clc 
getvariable.pl -c -0.03 "a(4,1)" m.clc 
getvariable.pl -c -0.02 "a(4,2)" m.clc 
getvariable.pl -c -0.004 "a(4,3)" m.clc 
getvariable.pl -c 0.002 "a(4,4)" m.clc 
getvariable.pl -c -0.0007 "a(6,-6)" m.clc 
getvariable.pl -c -0.007 "a(6,-5)" m.clc 
getvariable.pl -c -0.02 "a(6,-4)" m.clc 
getvariable.pl -c -0.04 "a(6,-3)" m.clc 
getvariable.pl -c -0.03 "a(6,-2)" m.clc 
getvariable.pl -c 0.005 "a(6,-1)" m.clc 
getvariable.pl -c 0.03 "a(6,0)" m.clc 
getvariable.pl -c 0.01 "a(6,1)" m.clc 
getvariable.pl -c -0.02 "a(6,2)" m.clc 
getvariable.pl -c -0.008 "a(6,3)" m.clc 
getvariable.pl -c 0.007 "a(6,4)" m.clc 
getvariable.pl -c 0.006 "a(6,5)" m.clc 
getvariable.pl -c 0.002 "a(6,6)" m.clc

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

getvariable.pl -c 0.6  "a(0,0)" m.clc 
getvariable.pl -c -0.05   "a(2,-2)" m.clc 
getvariable.pl -c -0.07   "a(2,-1)" m.clc 
getvariable.pl -c -0.09   "a(2,0)" m.clc 
getvariable.pl -c -0.1   "a(2,1)" m.clc 
getvariable.pl -c -0.03   "a(2,2)" m.clc 
getvariable.pl -c -0.006   "a(4,-4)" m.clc 
getvariable.pl -c -0.03   "a(4,-3)" m.clc 
getvariable.pl -c -0.04   "a(4,-2)" m.clc 
getvariable.pl -c -0.02   "a(4,-1)" m.clc 
getvariable.pl -c 0.02   "a(4,0)" m.clc 
getvariable.pl -c -0.04   "a(4,1)" m.clc 
getvariable.pl -c -0.03   "a(4,2)" m.clc 
getvariable.pl -c -0.001   "a(4,3)" m.clc 
getvariable.pl -c 0.003   "a(4,4)" m.clc 
getvariable.pl -c 0.00006   "a(6,-6)" m.clc 
getvariable.pl -c 0.002   "a(6,-5)" m.clc 
getvariable.pl -c 0.007   "a(6,-4)" m.clc 
getvariable.pl -c 0.01   "a(6,-3)" m.clc 
getvariable.pl -c 0.009   "a(6,-2)" m.clc 
getvariable.pl -c -0.002  "a(6,-1)" m.clc 
getvariable.pl -c -0.008   "a(6,0)" m.clc 
getvariable.pl -c -0.003   "a(6,1)" m.clc 
getvariable.pl -c 0.006   "a(6,2)" m.clc 
getvariable.pl -c 0.0007   "a(6,3)" m.clc 
getvariable.pl -c -0.003   "a(6,4)" m.clc 
getvariable.pl -c -0.002   "a(6,5)" m.clc 
getvariable.pl -c -0.0006   "a(6,6)" m.clc 


densplt s Pr3p_ic1ion_trunc.sipf 2 10 20 30 > m.clc

getvariable.pl -c -0.2  "aS1(0,0)" m.clc 
getvariable.pl -c 0.000000  "aS1(1,-1)" m.clc 
getvariable.pl -c 0.000000  "aS1(1,0)" m.clc 
getvariable.pl -c 0.000000  "aS1(1,1)" m.clc 
getvariable.pl -c 0.04  "aS1(2,-2)" m.clc 
getvariable.pl -c 0.01  "aS1(2,-1)" m.clc 
getvariable.pl -c -0.02  "aS1(2,0)" m.clc 
getvariable.pl -c 0.1  "aS1(2,1)" m.clc 
getvariable.pl -c 0.07  "aS1(2,2)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,-3)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,-2)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,-1)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,0)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,1)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,2)" m.clc 
getvariable.pl -c 0.000000  "aS1(3,3)" m.clc 
getvariable.pl -c 0.009  "aS1(4,-4)" m.clc 
getvariable.pl -c 0.03  "aS1(4,-3)" m.clc 
getvariable.pl -c 0.02  "aS1(4,-2)" m.clc 
getvariable.pl -c -0.004  "aS1(4,-1)" m.clc 
getvariable.pl -c -0.03  "aS1(4,0)" m.clc 
getvariable.pl -c 0.006  "aS1(4,1)" m.clc 
getvariable.pl -c 0.03  "aS1(4,2)" m.clc 
getvariable.pl -c 0.01  "aS1(4,3)" m.clc 
getvariable.pl -c -0.0002  "aS1(4,4)" m.clc 
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
getvariable.pl -c -0.002  "aS1(6,-6)" m.clc 
getvariable.pl -c -0.01  "aS1(6,-5)" m.clc 
getvariable.pl -c -0.03  "aS1(6,-4)" m.clc 
getvariable.pl -c -0.02  "aS1(6,-3)" m.clc 
getvariable.pl -c 0.005  "aS1(6,-2)" m.clc 
getvariable.pl -c 0.01  "aS1(6,-1)" m.clc 
getvariable.pl -c 0.008  "aS1(6,0)" m.clc 
getvariable.pl -c 0.03  "aS1(6,1)" m.clc 
getvariable.pl -c -0.007  "aS1(6,2)" m.clc 
getvariable.pl -c -0.02  "aS1(6,3)" m.clc 
getvariable.pl -c -0.001  "aS1(6,4)" m.clc 
getvariable.pl -c 0.005  "aS1(6,5)" m.clc 
getvariable.pl -c 0.002  "aS1(6,6)" m.clc 

getvariable.pl -c -0.1  "aS2(0,0)" m.clc 
getvariable.pl -c 0.000000  "aS2(1,-1)" m.clc 
getvariable.pl -c 0.000000  "aS2(1,0)" m.clc 
getvariable.pl -c 0.000000  "aS2(1,1)" m.clc 
getvariable.pl -c 0.07  "aS2(2,-2)" m.clc 
getvariable.pl -c 0.1  "aS2(2,-1)" m.clc 
getvariable.pl -c -0.01  "aS2(2,0)" m.clc 
getvariable.pl -c 0.01  "aS2(2,1)" m.clc 
getvariable.pl -c -0.03  "aS2(2,2)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,-3)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,-2)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,-1)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,0)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,1)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,2)" m.clc 
getvariable.pl -c 0.000000  "aS2(3,3)" m.clc 
getvariable.pl -c 0.0009  "aS2(4,-4)" m.clc 
getvariable.pl -c 0.02  "aS2(4,-3)" m.clc 
getvariable.pl -c 0.03  "aS2(4,-2)" m.clc 
getvariable.pl -c 0.01  "aS2(4,-1)" m.clc 
getvariable.pl -c -0.01  "aS2(4,0)" m.clc 
getvariable.pl -c -0.004  "aS2(4,1)" m.clc 
getvariable.pl -c -0.01  "aS2(4,2)" m.clc 
getvariable.pl -c -0.02  "aS2(4,3)" m.clc 
getvariable.pl -c -0.009  "aS2(4,4)" m.clc 
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
getvariable.pl -c 0.002  "aS2(6,-6)" m.clc 
getvariable.pl -c 0.005  "aS2(6,-5)" m.clc 
getvariable.pl -c -0.0008  "aS2(6,-4)" m.clc 
getvariable.pl -c -0.01  "aS2(6,-3)" m.clc 
getvariable.pl -c -0.008  "aS2(6,-2)" m.clc 
getvariable.pl -c 0.01  "aS2(6,-1)" m.clc 
getvariable.pl -c 0.005  "aS2(6,0)" m.clc 
getvariable.pl -c 0.01  "aS2(6,1)" m.clc 
getvariable.pl -c 0.01  "aS2(6,2)" m.clc 
getvariable.pl -c 0.03  "aS2(6,3)" m.clc 
getvariable.pl -c 0.03  "aS2(6,4)" m.clc 
getvariable.pl -c 0.01  "aS2(6,5)" m.clc 
getvariable.pl -c 0.002  "aS2(6,6)" m.clc 

getvariable.pl -c -0.4  "aS3(0,0)" m.clc 
getvariable.pl -c 0.000000  "aS3(1,-1)" m.clc 
getvariable.pl -c 0.000000  "aS3(1,0)" m.clc 
getvariable.pl -c 0.000000  "aS3(1,1)" m.clc 
getvariable.pl -c 0.01  "aS3(2,-2)" m.clc 
getvariable.pl -c 0.06  "aS3(2,-1)" m.clc 
getvariable.pl -c 0.1  "aS3(2,0)" m.clc 
getvariable.pl -c 0.1  "aS3(2,1)" m.clc 
getvariable.pl -c 0.009  "aS3(2,2)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,-3)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,-2)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,-1)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,0)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,1)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,2)" m.clc 
getvariable.pl -c 0.000000  "aS3(3,3)" m.clc 
getvariable.pl -c 0.002  "aS3(4,-4)" m.clc 
getvariable.pl -c 0.01  "aS3(4,-3)" m.clc 
getvariable.pl -c 0.03  "aS3(4,-2)" m.clc 
getvariable.pl -c 0.03  "aS3(4,-1)" m.clc 
getvariable.pl -c 0.001  "aS3(4,0)" m.clc 
getvariable.pl -c 0.05  "aS3(4,1)" m.clc 
getvariable.pl -c 0.02  "aS3(4,2)" m.clc 
getvariable.pl -c 0.0005  "aS3(4,3)" m.clc 
getvariable.pl -c -0.0008  "aS3(4,4)" m.clc 
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
getvariable.pl -c -0.001  "aS3(6,-5)" m.clc 
getvariable.pl -c -0.01  "aS3(6,-4)" m.clc 
getvariable.pl -c -0.03  "aS3(6,-3)" m.clc 
getvariable.pl -c -0.03  "aS3(6,-2)" m.clc 
getvariable.pl -c -0.006  "aS3(6,-1)" m.clc 
getvariable.pl -c 0.02  "aS3(6,0)" m.clc 
getvariable.pl -c -0.01  "aS3(6,1)" m.clc 
getvariable.pl -c -0.02  "aS3(6,2)" m.clc 
getvariable.pl -c -0.001  "aS3(6,3)" m.clc 
getvariable.pl -c 0.005  "aS3(6,4)" m.clc 
getvariable.pl -c 0.002  "aS3(6,5)" m.clc 
getvariable.pl -c 0.00006  "aS3(6,6)" m.clc 

densplt o Pr3p_xyz_ic1ion.sipf 2 20 10 30 > m.clc

getvariable.pl -c 1.45   "aL1(0,0)" m.clc 
getvariable.pl -c 0.000000   "aL1(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aL1(1,0)" m.clc 
getvariable.pl -c 0.000000   "aL1(1,1)" m.clc 
getvariable.pl -c -0.3   "aL1(2,-2)" m.clc 
getvariable.pl -c 0.01   "aL1(2,-1)" m.clc 
getvariable.pl -c 0.3   "aL1(2,0)" m.clc 
getvariable.pl -c -0.8   "aL1(2,1)" m.clc 
getvariable.pl -c -0.6   "aL1(2,2)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,0)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,1)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,2)" m.clc 
getvariable.pl -c 0.000000   "aL1(3,3)" m.clc 
getvariable.pl -c 0.00   "aL1(4,-4)" m.clc 
getvariable.pl -c -0.01   "aL1(4,-3)" m.clc 
getvariable.pl -c -0.05   "aL1(4,-2)" m.clc 
getvariable.pl -c -0.05   "aL1(4,-1)" m.clc 
getvariable.pl -c -0.01   "aL1(4,0)" m.clc 
getvariable.pl -c -0.07   "aL1(4,1)" m.clc 
getvariable.pl -c -0.01   "aL1(4,2)" m.clc 
getvariable.pl -c 0.02   "aL1(4,3)" m.clc 
getvariable.pl -c 0.01   "aL1(4,4)" m.clc 
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
getvariable.pl -c 0.00   "aL1(6,-6)" m.clc 
getvariable.pl -c 0.02   "aL1(6,-5)" m.clc 
getvariable.pl -c 0.05   "aL1(6,-4)" m.clc 
getvariable.pl -c 0.04   "aL1(6,-3)" m.clc 
getvariable.pl -c -0.01   "aL1(6,-2)" m.clc 
getvariable.pl -c -0.02   "aL1(6,-1)" m.clc 
getvariable.pl -c -0.01   "aL1(6,0)" m.clc 
getvariable.pl -c -0.05   "aL1(6,1)" m.clc 
getvariable.pl -c 0.01   "aL1(6,2)" m.clc 
getvariable.pl -c 0.03   "aL1(6,3)" m.clc 
getvariable.pl -c 0.01   "aL1(6,4)" m.clc 
getvariable.pl -c -0.01   "aL1(6,5)" m.clc 
getvariable.pl -c -0.00   "aL1(6,6)" m.clc 

getvariable.pl -c 0.7   "aL2(0,0)" m.clc 
getvariable.pl -c 0.000000   "aL2(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aL2(1,0)" m.clc 
getvariable.pl -c 0.000000   "aL2(1,1)" m.clc 
getvariable.pl -c -0.6   "aL2(2,-2)" m.clc 
getvariable.pl -c -0.8   "aL2(2,-1)" m.clc 
getvariable.pl -c 0.2   "aL2(2,0)" m.clc 
getvariable.pl -c 0.01   "aL2(2,1)" m.clc 
getvariable.pl -c 0.3   "aL2(2,2)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,0)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,1)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,2)" m.clc 
getvariable.pl -c 0.000000   "aL2(3,3)" m.clc 
getvariable.pl -c -0.00   "aL2(4,-4)" m.clc 
getvariable.pl -c 0.00   "aL2(4,-3)" m.clc 
getvariable.pl -c 0.01   "aL2(4,-2)" m.clc 
getvariable.pl -c -0.00   "aL2(4,-1)" m.clc 
getvariable.pl -c -0.01   "aL2(4,0)" m.clc 
getvariable.pl -c -0.05   "aL2(4,1)" m.clc 
getvariable.pl -c -0.05   "aL2(4,2)" m.clc 
getvariable.pl -c -0.04   "aL2(4,3)" m.clc 
getvariable.pl -c -0.01   "aL2(4,4)" m.clc 
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
getvariable.pl -c -0.01   "aL2(6,-6)" m.clc 
getvariable.pl -c -0.01   "aL2(6,-5)" m.clc 
getvariable.pl -c 0.01   "aL2(6,-4)" m.clc 
getvariable.pl -c 0.03   "aL2(6,-3)" m.clc 
getvariable.pl -c 0.02   "aL2(6,-2)" m.clc 
getvariable.pl -c -0.02   "aL2(6,-1)" m.clc 
getvariable.pl -c -0.01   "aL2(6,0)" m.clc 
getvariable.pl -c -0.02   "aL2(6,1)" m.clc 
getvariable.pl -c -0.02   "aL2(6,2)" m.clc 
getvariable.pl -c -0.05   "aL2(6,3)" m.clc 
getvariable.pl -c -0.05   "aL2(6,4)" m.clc 
getvariable.pl -c -0.02   "aL2(6,5)" m.clc 
getvariable.pl -c -0.00   "aL2(6,6)" m.clc 

getvariable.pl -c 2.17   "aL3(0,0)" m.clc 
getvariable.pl -c 0.000000   "aL3(1,-1)" m.clc 
getvariable.pl -c 0.000000   "aL3(1,0)" m.clc 
getvariable.pl -c 0.000000   "aL3(1,1)" m.clc 
getvariable.pl -c 0.01   "aL3(2,-2)" m.clc 
getvariable.pl -c -0.3   "aL3(2,-1)" m.clc 
getvariable.pl -c -1   "aL3(2,0)" m.clc 
getvariable.pl -c -0.6   "aL3(2,1)" m.clc 
getvariable.pl -c 0.00   "aL3(2,2)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,-3)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,-2)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,-1)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,0)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,1)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,2)" m.clc 
getvariable.pl -c 0.000000   "aL3(3,3)" m.clc 
getvariable.pl -c -0.02   "aL3(4,-4)" m.clc 
getvariable.pl -c -0.06   "aL3(4,-3)" m.clc 
getvariable.pl -c -0.07   "aL3(4,-2)" m.clc 
getvariable.pl -c -0.02   "aL3(4,-1)" m.clc 
getvariable.pl -c 0.05   "aL3(4,0)" m.clc 
getvariable.pl -c -0.03   "aL3(4,1)" m.clc 
getvariable.pl -c -0.05   "aL3(4,2)" m.clc 
getvariable.pl -c -0.01   "aL3(4,3)" m.clc 
getvariable.pl -c 0.01   "aL3(4,4)" m.clc 
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
getvariable.pl -c 0.00   "aL3(6,-5)" m.clc 
getvariable.pl -c 0.02   "aL3(6,-4)" m.clc 
getvariable.pl -c 0.05   "aL3(6,-3)" m.clc 
getvariable.pl -c 0.06   "aL3(6,-2)" m.clc 
getvariable.pl -c 0.01   "aL3(6,-1)" m.clc 
getvariable.pl -c -0.05   "aL3(6,0)" m.clc 
getvariable.pl -c 0.02   "aL3(6,1)" m.clc 
getvariable.pl -c 0.04   "aL3(6,2)" m.clc 
getvariable.pl -c 0.01   "aL3(6,3)" m.clc 
getvariable.pl -c -0.01   "aL3(6,4)" m.clc 
getvariable.pl -c -0.00   "aL3(6,5)" m.clc 
getvariable.pl -c 0.000000   "aL3(6,6)" m.clc 


singleion -sx -r Pr3p_xyz_ic1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.2 2 9 2 0 m.clc
getvalue.pl -c 0.04    2 13 2 0 m.clc
getvalue.pl -c 0.01   2 14 2 0 m.clc
getvalue.pl -c -0.02 2 15 2 0 m.clc
getvalue.pl -c 0.1  2 16 2 0 m.clc
getvalue.pl -c 0.07 2 17 2 0 m.clc
getvalue.pl -c 0.01   0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.00  0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.00   0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

singleion -sx -r Pr3p_ic1ion_trunc.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.2 2 9 2 0 m.clc
getvalue.pl -c 0.04  2 13 2 0 m.clc
getvalue.pl -c 0.01   2 14 2 0 m.clc
getvalue.pl -c -0.02 2 15 2 0 m.clc
getvalue.pl -c 0.1  2 16 2 0 m.clc
getvalue.pl -c 0.07 2 17 2 0 m.clc
getvalue.pl -c 0.02   0 11 4 0 results/Pr3p_ic1ion_trunc.sipf.trs
getvalue.pl -c 0.003  0 15 4 0 results/Pr3p_ic1ion_trunc.sipf.trs
getvalue.pl -c 0.0002   0 16 4 0 results/Pr3p_ic1ion_trunc.sipf.trs


singleion -sx -r Pr3p_xyz_icf1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c -0.2 2 9 2 0 m.clc
getvalue.pl -c 0.03 2 13 2 0 m.clc
getvalue.pl -c 0.02  2 14 2 0 m.clc
getvalue.pl -c 0.005  2 15 2 0 m.clc 
getvalue.pl -c 0.0924  2 16 2 0 m.clc
getvalue.pl -c 0.0447  2 17 2 0 m.clc
getvalue.pl -c 0.0162     0 11 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.00155    0 15 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
getvalue.pl -c 0.00049   0 16 4 0 results/Pr3p_xyz_icf1ion.sipf.trs
  
 
singleion -lx -r Pr3p_xyz_ic1ion.sipf  2 20 10 30 0 0 0 > m.clc
getvalue.pl -c 1     2 9 2 0 m.clc
getvalue.pl -c -0.3    2 13 2 0 m.clc
getvalue.pl -c 0.006  2 14 2 0 m.clc
getvalue.pl -c 0.3 2 15 2 0 m.clc
getvalue.pl -c -0.8  2 16 2 0 m.clc
getvalue.pl -c -0.6 2 17 2 0 m.clc
getvalue.pl -c 0.6    0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.1  0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
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
getvalue.pl -c -0.1    2 9 2 0 m.clc
getvalue.pl -c 0.07    2 13 2 0 m.clc
getvalue.pl -c 0.1     2 14 2 0 m.clc
getvalue.pl -c -0.01  2 15 2 0 m.clc 
getvalue.pl -c 0.01    2 16 2 0 m.clc
getvalue.pl -c -0.03   2 17 2 0 m.clc
    
getvalue.pl -c 0.02      0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.001    0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0007   0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

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
getvalue.pl -c  0.7    2 9 2 0 m.clc
getvalue.pl -c -0.6    2 13 2 0 m.clc
getvalue.pl -c -0.8     2 14 2 0 m.clc
getvalue.pl -c 0.2  2 15 2 0 m.clc 
getvalue.pl -c 0.01    2 16 2 0 m.clc
getvalue.pl -c 0.3   2 17 2 0 m.clc
        
getvalue.pl -c 0.8     0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.1      0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.05     0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

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
getvalue.pl -c -0.4   2 9 2 0 m.clc
getvalue.pl -c 0.01        2 13 2 0 m.clc
getvalue.pl -c 0.05    2 14 2 0 m.clc
getvalue.pl -c 0.1 2 15 2 0 m.clc 
getvalue.pl -c 0.1    2 16 2 0 m.clc
getvalue.pl -c 0.01   2 17 2 0 m.clc
    
getvalue.pl -c 0.008       0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 0.0002  0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.004   0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

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
getvalue.pl -c 2.2     2 9 2 0 m.clc
getvalue.pl -c 0.01      2 13 2 0 m.clc
getvalue.pl -c -0.3    2 14 2 0 m.clc
getvalue.pl -c -1  2 15 2 0 m.clc 
getvalue.pl -c -0.6    2 16 2 0 m.clc
getvalue.pl -c 0.01   2 17 2 0 m.clc
    
getvalue.pl -c 0.3     0 11 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c 6e-05     0 15 4 0 results/Pr3p_xyz_ic1ion.sipf.trs
getvalue.pl -c  0.1   0 16 4 0 results/Pr3p_xyz_ic1ion.sipf.trs

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
