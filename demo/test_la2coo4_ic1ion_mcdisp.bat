cd ../examples/la2coo4

mcphasit

getvalue.pl -c -469460e1 1 8 0.3 0 results/mcphas.fum

call setup_mcdisp_mf 0.3 0 0 0

REM "editing mcdisp.mf and insert additional anisotropy field Ha ..."

 call factcol 1 1.0122453 mcdisp.mf
 call factcol 2 1.0122453 mcdisp.mf

REM initialize mcdisp
 mcdispit -c -maxE 400
REM select only strong transitions (gamma large)
call range 7 0.0001 10000 results/mcdisp.trs
REM restart mcdisp to calculate intensities
call mcdispit -t -d -maxE 400

range 6	0.34 0.36 results/mcdisp.qei
getvalue.pl -c 0.629 9 10 33.6062 0  results/mcdisp.qei
getvalue.pl -c 0.265  9 13 33.6062 0  results/mcdisp.qei
getvalue.pl -c -0.378 9 15 33.6062 0  results/mcdisp.qei


cd ../../demo
