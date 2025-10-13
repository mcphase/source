cd ../examples/Ru3p_create_sipf

formfactor Ru3p.sipf
getvalue.pl -c 0.90033 1 2 1.0 0 results/formfactor.out
getvalue.pl -c 0.81059 1 3 1.0 0 results/formfactor.out
getvalue.pl -c 0.90033 1 4 1.0 0 results/formfactor.out
getvalue.pl -c 0.03620 1 5 1.0 0 results/formfactor.out

cd ../coo

mcdispit -c -max 200  -minE -0.1 -maxE 3000 -prefix 002
call range 6 -0.1 3000 results/002mcdisp.trs
getvalue.pl -c 113.691 5 6 10 0  results/002mcdisp.trs
 mcdispit -t -prefix 002
range 5 13 14 results/002mcdisp.qei
getvalue.pl -c 4e-06  9 10 1.33376 0 results/002mcdisp.qei
getvalue.pl -c 0.0001 9 11 1.33376 0 results/002mcdisp.qei
   
 

call mcdiff coo_exp.msf
call factcol 10 4 results/mcdiff.out
call factcol 11 4 results/mcdiff.out

getvalue.pl -c  0.028 4 10 5.3263E-01 0 results/mcdiff.out
getvalue.pl -c  0.008 4 11 5.3263E-01 0 results/mcdiff.out



cd ../../demo
