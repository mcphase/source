cd ../examples/bfk

call singleion -r Pr3p.sipf -nt 100 10 0 0 0 0 0 0
call formfactor Pr3p.sipf

call bfk 10 0.03 2 1 results/Pr3p.sipf.levels.cef bfk.par

getvalue.pl -c 0.30186  2 3 3.0000 0 results/bfk2.res


cd ../../demo
