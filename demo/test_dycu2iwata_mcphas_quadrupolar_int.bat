cd ../examples/dycu2iwata


copy mcphas.ini test_mcphas.ini
setvariable xmax 2 test_mcphas.ini


call mcphasit -prefix test_
getvalue.pl -c 9.892 2 10 10 0 results/test_mcphas.fum
getvalue.pl -c -12.4180  2 9 10 0 results/test_mcphas.fum
     

cd ../../demo
