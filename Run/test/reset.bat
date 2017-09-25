rm ESC fort.* logfile
cp fort10 fort.10
cp ../../Ncode/nbody6mp .
./nbody6mp < input > logfile &
