cd ../examples/tbcu2a

mcphasit  -prefix test_ -v

getvalue.pl -c -346.0 2 8 3 0 results/test_mcphas.fum
getvalue.pl -c 2.94 2 10 3 0 results/test_mcphas.fum

setup_mcdiff_in -prefix test_ 5 0 0 0

setvariable thetamax 15 test_mcdiff.in

mcdiff -prefix test_

range 1 -0.001 1000 results/test_mcdiff.out
range 2 -0.001 1000 results/test_mcdiff.out
range 3 -0.001 1000 results/test_mcdiff.out

getvalue.pl -c 5.6861E-02 6 8 2.7730E+01 0 results/test_mcdiff.out
getvalue.pl -c 5.6878E-02 6 12 2.7730E+01 0 results/test_mcdiff.out

setup_mcdisp_mf -prefix test_ 5 6 0 0

mcdispit -maxE 30 -prefix test_

getvalue.pl -c 0.6567 9 10  6.07745  0 results/test_mcdisp.qei   
getvalue.pl -c 0.6572 9 11  6.07745  0 results/test_mcdisp.qei   

# add some spin charge moment densities

# spins -c -M -prefix test_ 5 6 0 0  0 0 1 5.76932 
# spins -s -M -prefix test_ 5 6 0 0  0 0 1 5.76932 

cd ../../demo
