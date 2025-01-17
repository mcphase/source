cd ../examples/prni2b2c



REM do a mcphas run to get stable structure
call mcphasit

getvalue.pl -c -16.616 1 8 15 0 results/mcphas.fum

call powdermagnon 0.2 0.2 0.05 3 0.2 2

getvalue.pl -c 0.00357 0 1 7 0 powdermagnon.hkl

REM Ei=3.15 meV ... ki=1.238 A^-1

call substitute "kf=2" "ki=1.238" mcdisp.par

call setup_mcdisp_mf 2 0 0 0

call mcdispit -minE 0.1 -max 2

getvalue.pl -c 0.5907 0 9 9 0 results/mcdisp.qei
getvalue.pl -c 0.00287 0 10 9 0 results/mcdisp.qei
getvalue.pl -c 0.00287 0 11 9 0 results/mcdisp.qei

cd ../../demo
