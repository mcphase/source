cd ../examples/upd3

call mcphasit -v
getvalue.pl -c -21.491 0 9 1 0 results/mcphas.fum
call spins -f results/mcphas.mf 3 0 0 0 > mcdisp.mf
call mcdispit -c -minE 0 -maxE 20
call range 7 1e-10 1e10 results/mcdisp.trs
call mcdispit -t -d
getvalue.pl -c 2.614 1 11 0.008139 0  results/mcdisp.qom

cd ../../demo


