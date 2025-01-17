cd ../examples/tbcu2a

cp mcphas.ini test_mcphas.ini
setvariable xmax 5 test_mcphas.ini
mcphasit -prefix test_

getvalue.pl -c -22.0805 2 8 3 0 results/test_mcphas.fum
getvalue.pl -c 2.859 2 10 3 0 results/test_mcphas.fum



cd ../../demo
