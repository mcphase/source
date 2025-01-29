cd ../examples/ndcu2b_1dim
mcphasit -v
range 1 1.9 2.1 results/mcphas.fum
getvalue.pl -c 0.85 4 10 2.2 0 results/mcphas.fum
getvalue.pl -c 1.35 4 10 2.6 0 results/mcphas.fum
getvalue.pl -c 2.533 4 10 2.8 0 results/mcphas.fum
cd ../cecu2a/
mcdispit
getvalue.pl -c  1.143 0 15 1 0 results/mcdisp.qom
getvalue.pl -c  1.835 0 16 1 0 results/mcdisp.qom
getvalue.pl -c 1.40466 0 6 2 0 results/mcdisp.trs
getvalue.pl -c 0.055165 0 8 2 0  results/mcdisp.trs
cd ../../demo

