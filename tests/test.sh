echo "Please execute this file from the pcasuite main dir"

echo "\n Testing pcazip:"
./pcazip -i tests/2hk5.netcdf -p tests/2hk5.pdb -o tests/2hk5.pcz -v -M @C,CA

echo "\n Testing pcadump:"
./pczdump -i tests/2hk5.pcz --info

echo "\n Testing pcaunzip:"
./pcaunzip -i tests/2hk5.pcz -o tests/2hk5.x -v
