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
densplt c -M Pr3p.sipf 5 1000 0 0 
rm Pr3p.sipf Ma.clc res.dat


