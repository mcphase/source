cd ../examples/tungsten_phonons

echo "create mcphas.j with interactions and sipf file K matrices from "
echo "Born van Karman springs in bvk_springs.dat"
copy mcphas0.j mcphas.j
call makenn 4.7 -bvk bvk_springs.dat 
copy results/makenn.j mcphas.j
getvalue.pl -c +2.7442835e+02 0 4 1 0 mcphas.j

call mcdispit
range 1 9.971 9.973 results/mcdisp.qei
getvalue.pl -c 0.002343 9 12 19.458       0   results/mcdisp.qei


cd ../../demo
