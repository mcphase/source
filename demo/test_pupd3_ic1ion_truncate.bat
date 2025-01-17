cd ../examples/pupd3

REM in mcphas 5.5 there was definitely something buggy in
REM output - see emails with Duc autumn 2024 and some bugs were removed
REM (ic1ion and icf1ion were not consistent) 
REM 16.1.2025 tested with respect to XTLS on pu3p.sipf 
REM with alpha=beta=gamma=0  see 

ic1ion pu3p.sipf
getvariable.pl -c -7496.1 E0 results/ic1ion.out
getvalue.pl -c 0  0 1 1 0  results/ic1ion.out
getvalue.pl -c 0  0 1 2 0  results/ic1ion.out
getvalue.pl -c 12.2  0 1 3 0  results/ic1ion.out
getvalue.pl -c 12.2  0 1 4 0  results/ic1ion.out
getvalue.pl -c 12.2  0 1 5 0  results/ic1ion.out
getvalue.pl -c 12.2  0 1 6 0  results/ic1ion.out
getvalue.pl -c 322  0 1 7 0  results/ic1ion.out
getvalue.pl -c 322  0 1 8 0  results/ic1ion.out
getvalue.pl -c 490  0 1 9 0  results/ic1ion.out
getvalue.pl -c 490  0 1 10 0  results/ic1ion.out
getvalue.pl -c 540  0 1 11 0  results/ic1ion.out
getvalue.pl -c 540  0 1 12 0  results/ic1ion.out
getvalue.pl -c 540  0 1 13 0  results/ic1ion.out
getvalue.pl -c 540  0 1 14 0  results/ic1ion.out



icf1ion pu3p.sipf
getvariable.pl -c -1121.8 E0 results/icf1ion.out
getvalue.pl -c 0  0 1 1 0  results/icf1ion.out
getvalue.pl -c 0  0 1 2 0  results/icf1ion.out
getvalue.pl -c 8.31  0 1 3 0  results/icf1ion.out
getvalue.pl -c 8.31  0 1 4 0  results/icf1ion.out
getvalue.pl -c 8.31  0 1 5 0  results/icf1ion.out
getvalue.pl -c 8.31  0 1 6 0  results/icf1ion.out
getvalue.pl -c 81.9  0 1 7 0  results/icf1ion.out
getvalue.pl -c 81.9  0 1 8 0  results/icf1ion.out
getvalue.pl -c 374  0 1 9 0  results/icf1ion.out
getvalue.pl -c 374  0 1 10 0  results/icf1ion.out
getvalue.pl -c 381  0 1 11 0  results/icf1ion.out
getvalue.pl -c 381  0 1 12 0  results/icf1ion.out
getvalue.pl -c 458  0 1 13 0  results/icf1ion.out
getvalue.pl -c 458  0 1 14 0  results/icf1ion.out


cp mcphas.ini test_mcphas.ini
setvariable xmax 1 test_mcphas.ini
mcphasit -prefix test_

getvalue.pl -c -7497.1447 1 8 1 0 results/test_mcphas.fum
getvalue.pl -c 0.00116 1 10 1 0 results/test_mcphas.fum


cd ../../demo
