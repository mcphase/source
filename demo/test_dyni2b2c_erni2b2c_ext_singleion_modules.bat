cd ../examples/dyni2b2c

copy mcphas.ini test_mcphas.ini
setvariable xmax 4  test_mcphas.ini
mcphasit -prefix test_
getvalue.pl -c -0.504516 2 8 0 0 results/test_mcphas.fum

cd ../erni2b2c

copy mcphas.ini test_mcphas.ini
setvariable xmax 2  test_mcphas.ini
setvariable ymax 0  test_mcphas.ini
setvariable maxnofspins 30  test_mcphas.ini
mcphasit -prefix test_
getvalue.pl -c -8.63607 2 8 0 0 results/test_mcphas.fum


cd ../../demo
