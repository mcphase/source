cd ../examples/ndba2cu3o7

REM this calculates the <h2> anisotropic formfactor for Nd3p in NBCO </h2>

copy mcdiff_ff.in mcdiff.in
mcdiff 


echo "*********************************************************************"
echo "col10 contains 2|MSF.P|/sin^2(angle between Q and P)  ~ M(Q)[mb]"
echo "           here we assume that M(Q)||P !!!"
echo "*********************************************************************"

REM REMove peaks which are far out of scattering plane
call range -d 3 -5 5 results/mcdiff.out
average -dmin=0.0001 -av 5  results/mcdiff.out
getvalue.pl -c 1.329 4 10 3.912 0 results/mcdiff.out
getvalue.pl -c 1.328 4 10 2.7449 0 results/mcdiff.out


cd ../../demo