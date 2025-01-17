cd ../examples/cepd2si2
cp rkky.j mcphas.j
makenn 30

addj -s -1 results/makenn.j rkky.j > mcphas.j
getvalue.pl -c -5.98781e-02 0 4 6 0 mcphas.j
mcphasit
getvalue.pl -c -14.6456 1 8 1.3 0 results/mcphas.fum
setup_mcdiff_in 1.3 0 0 0 

copy exp_int_mag.dat scaled_exp_int_mag.dat
call delcol 4 scaled_exp_int_mag.dat
call delcol 4 scaled_exp_int_mag.dat
call delcol 4 scaled_exp_int_mag.dat
call delcol 4 scaled_exp_int_mag.dat

call factcol 4 0.0000064 scaled_exp_int_mag.dat
call factcol 5 0.0000064 scaled_exp_int_mag.dat

call mcdiff scaled_exp_int_mag.dat
getvalue.pl -c 1.7066e-03 0 5 30 0 results/mcdiff.out
getvalue.pl -c 1.2511e-03 0 7 30 0 results/mcdiff.out

getvariable.pl -c 9.0  rpvaluecol5 results/mcdiff.out
cd ../../demo
