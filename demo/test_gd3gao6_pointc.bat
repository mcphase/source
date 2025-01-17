cd ../examples/gd3gao6

call makenn 3 -kaneyoshi 0 0.1 1 -d

call pointc Gd3p.sipf results/makenn.a1.pc > Gd3p1.sipf
call pointc Gd3p.sipf results/makenn.a2.pc > Gd3p2.sipf
call pointc Gd3p.sipf results/makenn.a3.pc > Gd3p3.sipf
call pointc Gd3p.sipf results/makenn.a4.pc > Gd3p4.sipf
call pointc Gd3p.sipf results/makenn.a5.pc > Gd3p5.sipf
call pointc Gd3p.sipf results/makenn.a6.pc > Gd3p6.sipf


getvariable.pl -c -46.407 L22S Gd3p1.sipf 
getvariable.pl -c 1.5627 L21S Gd3p1.sipf 
getvariable.pl -c 5.4692 L20 Gd3p1.sipf 
getvariable.pl -c 5.0149 L21 Gd3p1.sipf 
getvariable.pl -c 9.3304 L22 Gd3p1.sipf 
getvariable.pl -c -1.0071 L44S Gd3p1.sipf 
getvariable.pl -c 15.227 L43S Gd3p1.sipf 
getvariable.pl -c -8.5108 L42S Gd3p1.sipf 
getvariable.pl -c -2.6764 L41S Gd3p1.sipf 
getvariable.pl -c -23.191 L40 Gd3p1.sipf 
getvariable.pl -c 2.9253 L41 Gd3p1.sipf 
getvariable.pl -c -0.63880 L42 Gd3p1.sipf 
getvariable.pl -c -1.6193 L43 Gd3p1.sipf 
getvariable.pl -c -7.6234 L44 Gd3p1.sipf 
getvariable.pl -c 0.21048 L66S Gd3p1.sipf 
getvariable.pl -c 0.058662 L65S Gd3p1.sipf 
getvariable.pl -c 0.21712 L64S Gd3p1.sipf 
getvariable.pl -c 1.5212 L63S Gd3p1.sipf 
getvariable.pl -c 1.275 L62S Gd3p1.sipf 
getvariable.pl -c 0.1922 L61S Gd3p1.sipf 
getvariable.pl -c 2.3180 L60 Gd3p1.sipf 
getvariable.pl -c -0.16464 L61 Gd3p1.sipf 
getvariable.pl -c -0.30133 L62 Gd3p1.sipf 
getvariable.pl -c -0.3411 L63 Gd3p1.sipf 
getvariable.pl -c -2.4409 L64 Gd3p1.sipf 
getvariable.pl -c 2.288 L65 Gd3p1.sipf 
getvariable.pl -c -0.71918 L66 Gd3p1.sipf 



cd ../../demo
