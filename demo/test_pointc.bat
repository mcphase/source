# Test of pointc
pointc Pr3+ 2.4 0 0 3  >  Pr3p.sipf
getvariable.pl -c 0.8 GJ Pr3p.sipf
getvariable.pl -c 4.505 B20 Pr3p.sipf
getvariable.pl -c -37.1829 L40  Pr3p.sipf
getvariable.pl -c -6.37011 L60  Pr3p.sipf
# singleion
singleion -r Pr3p.sipf 5 1000 0 0  0 0 0 > Ma.clc
getvalue.pl -c 3.91584 3 9 1000 0 Ma.clc
getvalue.pl -c 78.8805 0 6 2    0 results/Pr3p.sipf.trs
getvalue.pl -c 4.47757  0 7 2    0 results/Pr3p.sipf.trs
getvalue.pl -c 0.138796 0 8 2    0 results/Pr3p.sipf.trs
getvariable.pl -c -253.893 Eigenvalues results/Pr3p.sipf.levels.cef
# gauss
gauss 2 0.2 -4 4 > res.dat
# convolute (must be .pl because of eval in convolute command batch)
convolute.pl  5 7 results/Pr3p.sipf.trs 1 2 res.dat 
# test of display_density
densplt c -M Pr3p.sipf 5 1000 0 0 > dd
getvariable.pl -c 0.564198   "a(0,0)" dd 
getvariable.pl -c 0.000000   "a(2,-2)" dd 
getvariable.pl -c 0.000000   "a(2,-1)" dd 
getvariable.pl -c 0.109704   "a(2,0)" dd 
getvariable.pl -c 0.000000   "a(2,1)" dd 
getvariable.pl -c -0.139594   "a(2,2)" dd 
getvariable.pl -c 0.000000   "a(4,-4)" dd 
getvariable.pl -c 0.000000   "a(4,-3)" dd 
getvariable.pl -c 0.000000   "a(4,-2)" dd 
getvariable.pl -c 0.000000   "a(4,-1)" dd 
getvariable.pl -c -0.043045   "a(4,0)" dd 
getvariable.pl -c 0.000000   "a(4,1)" dd 
getvariable.pl -c 0.047341   "a(4,2)" dd 
getvariable.pl -c 0.000000   "a(4,3)" dd 
getvariable.pl -c -0.024525   "a(4,4)" dd 
getvariable.pl -c 0.000000   "a(6,-6)" dd 
getvariable.pl -c 0.000000   "a(6,-5)" dd 
getvariable.pl -c 0.000000   "a(6,-4)" dd 
getvariable.pl -c 0.000000   "a(6,-3)" dd 
getvariable.pl -c 0.000000   "a(6,-2)" dd 
getvariable.pl -c 0.000000   "a(6,-1)" dd 
getvariable.pl -c -0.020352   "a(6,0)" dd 
getvariable.pl -c 0.000000   "a(6,1)" dd 
getvariable.pl -c 0.022013   "a(6,2)" dd 
getvariable.pl -c 0.000000   "a(6,3)" dd 
getvariable.pl -c -0.009557   "a(6,4)" dd 
getvariable.pl -c 0.000000   "a(6,5)" dd 
getvariable.pl -c 0.002600   "a(6,6)" dd 

singleion -r Pr3p.sipf -MQ 1 0 1 5 1000 0 0  0 0 0 > Ma.clc

getvalue.pl -c 3.13267 3 12 1000 0 Ma.clc



rem rm Pr3p.sipf Ma.clc res.dat dd


