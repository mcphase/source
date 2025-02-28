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

singleion -r W2p1.sipf -HE -P 2 0 0 0 10 0 0   0 0 0  > d.clc
getvalue.pl -c 0.000115 6 12 10 0 d.clc

singleion -r W2p1.sipf -HE -pel 2 0 0 0 10 0 0   0 0 0 > d.clc
getvalue.pl -c 0.02307 6 12 10 0 d.clc
getvalue -c 22.902 0 8 1 0 results/W2p1.sipf.trs
getvalue -c 17.036 0 11 1 0 results/W2p1.sipf.trs
getvalue -c 5.866 0 12 1 0 results/W2p1.sipf.trs


rm d.clc


cd ../../demo
