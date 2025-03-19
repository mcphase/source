cd ../examples/prni2si2

mcphasit

getvalue.pl -c -12.61782  2 8 8.1 0 results/mcphas.fum
getvalue.pl -c 1.425  2 10 8.1 0 results/mcphas.fum

setup_mcdisp_mf 4 0 0 0.1
 mcdispit -c
call range -d 6 -7 7 results/mcdisp.trs
call range -d 7 0.1 100 results/mcdisp.trs
call mcdispit -t

range 7 1.059 1.061 results/mcdisp.qei

getvalue.pl -c 0.4 9 10 1.175 0 results/mcdisp.qei
getvalue.pl -c 0.4 9 11 1.175 0 results/mcdisp.qei
getvalue.pl -c 0.02 9 13 1.175 0 results/mcdisp.qei
getvalue.pl -c 0.02 9 15 1.175 0 results/mcdisp.qei
getvalue.pl -c -0.1 9 17 1.175 0 results/mcdisp.qei


cd ../../demo
