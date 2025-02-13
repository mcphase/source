cd ../examples/tbcu2a

mcphasit -prefix test_

getvalue.pl -c -345.96 2 8 3 0 results/test_mcphas.fum
getvalue.pl -c 2.915 2 10 3 0 results/test_mcphas.fum

setup_mcdisp_mf -prefix test_ 5 1 0 0


mcdispit -maxE 100 -prefix test_

cd ../../demo
