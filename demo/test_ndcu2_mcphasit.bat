# Test of mcphasit
cd ../examples/ndcu2b_new
mcphasit
getvalue.pl -c -6.23770 2  8 0   0 results/mcphas.fum
getvalue.pl -c 0.7261   2 10 2.2 0 results/mcphas.fum
getvalue.pl -c 1.081    2 10 3.0 0 results/mcphas.fum
getvalue.pl -c 2.073    2 10 3.2 0 results/mcphas.fum
setup_mcdiff_in 1 0 0 0 
mcdiff
getvalue.pl -c 1.9671E-02 4 8 2.1922E+01 0 results/mcdiff.out
formfactor Nd3p.sipf
getvalue.pl -c 0.9205278 1 2 2.4 0 results/formfactor.out
radwavfunc Nd3p.sipf
getvalue.pl -c 0.09074942 1 2 0.017958563  0 results/radwavfunc.out
 
setup_mcdisp_mf 1.5 0 0 0
cp mcdispall.par mcdisp.par
mcdispit -max 2 
getvalue.pl -c  -1.787 5 10 0.02 0 results/mcdisp.qom
range -d 5 0.039 0.041 results/mcdisp.qei
getvalue.pl -c  0.0408399 9 10 1.24943 0 results/mcdisp.qei
cd ../../demo




