cd ../examples/testic1ion
call pointc testic1ion.sipf pointcharges_octahedron.dat > octahedronic1ion.sipf
getvariable.pl -c 91.118 L64 results/pointc.Llm
call substitute "MODULE=so1ion" "MODULE=ic1ion" octahedronic1ion.sipf
call densplt c octahedronic1ion.sipf 2 0 0 0 
call ic1ion octahedronic1ion.sipf
getvalue.pl -c 188.112 0 1 9 0 results/ic1ion.out
densplt c -M Co2p_atom1_rotated_pm_z.sipf 2 0 0 100
singleion -M -r Co2p_atom1_rotated_pm_z.sipf 2 10 0 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 0.06930 3 12 10 0 m.clc
singleion -M -r Co2p_atom1_rotated_pm_z.sipf 2 0 10 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 0.0698 4 13 10 0 m.clc
singleion -M -r Co2p_atom1_rotated_pm_z.sipf 2 0 0 10  0 0 0 0 0 0 > m.clc
getvalue.pl -c 5.998 5 14 10 0 m.clc

call substitute "MODULE=ic1ion" "MODULE=icf1ion" octahedronic1ion.sipf
call densplt c octahedronic1ion.sipf 2 0 0 0 
call icf1ion octahedronic1ion.sipf
getvalue.pl -c 173.306 0 1 9 0 results/icf1ion.out
singleion -M -r octahedronic1ion.sipf 2 10 0 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 1.9939 3 12 10 0 m.clc
singleion -M -r octahedronic1ion.sipf 2 0 10 0  0 0 0 0 0 0 > m.clc
getvalue.pl -c 1.9939 4 13 10 0 m.clc
singleion -M -r octahedronic1ion.sipf 2 0 0 10  0 0 0 0 0 0 > m.clc
getvalue.pl -c 1.9939 5 14 10 0 m.clc
rm m.clc

cd ../../demo
