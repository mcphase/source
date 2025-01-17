cd ../examples/NiO

echo "First we set up the model by Hutchings"

echo "insert anisotropy parameters to sipf file Ni2p_rot.sipf"
echo " ... D1=1.13,D2=0.06 K ...convert to meV"
call setvariable B20 "0.0863x(1.13+0.06)/(-6)"  Ni2p_rot.sipf
call setvariable B22 "0.0862x(1.13-0.06)/2" Ni2p_rot.sipf

echo "rotate sipf file coordinate system to mcphase convention xyz||abc"
echo "first rotate around (1-10) so that 11-2 becomes z axis ..."

call rotateBlm -th +139.1066 -i Ni2p_rot.sipf -o results/dummy
echo "then rotate around z such that 1-10 becomes y axis"
call rotateBlm -fi 135 -i results/dummy -o Ni2p.sipf
getvariable.pl -c 0.02523 B22S Ni2p.sipf
getvariable.pl -c 0.0682 B21S Ni2p.sipf
getvariable.pl -c 0.003769 B20 Ni2p.sipf
getvariable.pl -c 0.06821 B21 Ni2p.sipf

echo "check anisotropy around 1-10 direction - easy axis should be at 311deg"
call anisotropyit 2 1 1 -1 0 60 -r Ni2p.sipf 0 0 0
getvalue.pl -c 0.89 8 13 312 0 results/anisotropy.out

echo "crystal structure is stored in mcphas0.j"
echo " ... here we generate the neighbour list ... "
copy mcphas0.j mcphas.j
call makenn 5 -e > exchng.table
call sortf 1 exchng.table

getvalue.pl -c 2.95227 0 1 2 0 exchng.table

echo "Now we insert the magnetic coupling"
echo "insert J1+ into exchng.table"
call setvalue 1 2 15.7 exchng.table
echo "insert J1- into exchng.table"
call setvalue 2 2 16.1 exchng.table
echo "insert J2 into exchng.table"
call setvalue 3 2 -221 exchng.table
echo "convert from K to meV"
call factcol 2 0.086277 exchng.table

call makenn 5 -e exchng.table
getvalue.pl -c +1.3545 0 4 4 0 results/makenn.j
copy results/makenn.j mcphas.j


echo "... calculate selfconsistently the magnetic structure"
call mcphasit -v
getvalue.pl -c -57.3223 2 8 0 0 results/mcphas.fum

echo "for spinwave calculate take magnetic structure which was calculated for T=2 K H=0"
call setup_mcdisp_mf 2 0 0 1

echo "calculate the spinwaves and compare to literature"
call mcdispit -max 3 -prefix 001
range 9 23 23.6 results/001mcdisp.qei
getvalue.pl -c 23.4588  1 9 1.376 0 results/001mcdisp.qei
getvalue.pl -c 0.6807  1 10 1.376 0 results/001mcdisp.qei
getvalue.pl -c 0.4784   1 13 1.376 0 results/001mcdisp.qei
getvalue.pl -c -0.01514 1 17 1.376 0 results/001mcdisp.qei
getvalue.pl -c 0.5428  1 21 1.376 0 results/001mcdisp.qei

rem <h3> now do the rixs cross section at Q=(0.75 0.75 0.75) </h3>

rem <h3> (i) sigma1r=1=sigmat1u (haverkort PRL 105 -167404 equ (10)), all other sigma=0</h3>
call setvariable SIGMA0r 0 Ni2p.sipf
call setvariable SIGMA1r 1 Ni2p.sipf
call setvariable SIGMA2r 0 Ni2p.sipf
call setvariable SIGMA0i 0 Ni2p.sipf
call setvariable SIGMA1i 0 Ni2p.sipf
call setvariable SIGMA2i 0 Ni2p.sipf
call mcdispit -max 3 -prefix 002 -xa 3
call range -d 9 0.1 1000 results/002mcdisp.qex

rem <h3> (ii) sigma2r=1=sigmat2eg=sigma2t2g (haverkort PRL 105 -167404 equ (10)), all other sigma=0</h3>
call setvariable SIGMA0r 0 Ni2p.sipf
call setvariable SIGMA1r 0 Ni2p.sipf
call setvariable SIGMA2r 1 Ni2p.sipf
call setvariable SIGMA0i 0 Ni2p.sipf
call setvariable SIGMA1i 0 Ni2p.sipf
call setvariable SIGMA2i 0 Ni2p.sipf
call mcdispit -max 3 -prefix 002 -xa 3
call range -d 9 0.1 1000 results/002mcdisp.qex
delline 100 1000000  results/002mcdisp.qex
getvalue.pl -c 9.4926E-01 11 10 30 0 results/002mcdisp.qex
getvalue.pl -c 5.4766E-02 13 12 30 0 results/002mcdisp.qex
getvalue.pl -c 1.8716E-01 15 14 30 0 results/002mcdisp.qex
getvalue.pl -c 2.2818E-02 17 16 30 0 results/002mcdisp.qex
getvalue.pl -c 2.9025E-01 19 18 30 0 results/002mcdisp.qex
getvalue.pl -c 7.5980E-01 21 20 30 0 results/002mcdisp.qex
getvalue.pl -c 6.6911E-02 23 22 30 0 results/002mcdisp.qex
getvalue.pl -c 9.7045E-02 25 24 30 0 results/002mcdisp.qex

cd ../../demo
