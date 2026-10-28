cd ../examples/tungsten_phonons

echo "create mcphas.j with interactions and sipf file K matrices from "
echo "Born van Karman springs in bvk_springs.dat"
copy mcphas0.j mcphas.j
call makenn 4.7 -bvk bvk_springs.dat 
copy results/makenn.j mcphas.j
getvalue.pl -c 274.40469 0 4 1 0 mcphas.j

call mcdispit
range -d 1 9.971 9.973 results/mcdisp.qei
getvalue.pl -c 0.002397 9 12 19.1897       0   results/mcdisp.qei

singleion -r W2p1.sipf -HE -P 2 0 0 0 10 0 0   0 0 0  > d.clc
getvalue.pl -c 0.000119 6 12 10 0 d.clc

singleion -r W2p1.sipf -HE -pel 2 0 0 0 10 0 0   0 0 0 > d.clc
getvalue.pl -c 0.02372 6 12 10 0 d.clc
getvalue.pl -c 23.222 0 13 3 0 results/W2p1.sipf.trs
getvalue.pl -c 17.812 0 11 1 0 results/W2p1.sipf.trs
getvalue.pl -c 5.411 0 12 1 0 results/W2p1.sipf.trs


singleion -r W2p1.sipf -pel -HE  2 0 0 0 0 0 1  0 0 0 > d.clc
getvalue.pl -c 0.0023 2 14 2 0 d.clc

singleion -r W2p1.sipf -Xpel 0 0  -HE   2 0 0 0 0 0 0  0 0 0 > d.clc
getvalue.pl -c 2.372	2 13 2 0 d.clc




rm d.clc


cd ../../demo
