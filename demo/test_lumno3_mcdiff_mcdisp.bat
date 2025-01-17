cd ../examples/lumno3

mcphasit -v

echo "***********************************************"
echo "LuMnO3 MCDIFF TO CALCULATE MAGNETIC DIFFRACTION"
echo "***********************************************"

 call setup_mcdiff_in 2 0 0 0
REM ... edit mcdiff.in to put in epxerimental setup (wavelength, ...)
 call setvariable thetamax 30 mcdiff.in
 mcdiff
average -dmin=0.1 -av 5 results/mcdiff.out
getvalue.pl -c 0.1196 4 8 5.35285  0  results/mcdiff.out
echo "***************************************"
echo "LuMnO3 MCDISP TO CALCULATE EXCITATIONS"
echo "***************************************"


  call setup_mcdisp_mf 2 0 0 0

REM initialize mcdisp
 mcdispit -max 3 -c
REM select only strong transitions (gamma large)
call range 7 0.1 3000 results/mcdisp.trs
 mcdispit -t -d
range 5 1.019 1.021 results/mcdisp.qei
range 6 -0.00001 0.0001  results/mcdisp.qei
range 7 -0.00001 0.0001 results/mcdisp.qei
getvalue.pl -c 0.654 9 10  6.52607  0 results/mcdisp.qei
  
cd ../../demo

