cd ../examples/gdru2si2

mcphasit mcphas.xyt
range 5 0.09 0.11 results/mcphas.fum
getvalue.pl -c -4.564 3 8 15 0 results/mcphas.fum

cd ../../demo
