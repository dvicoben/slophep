## Hammer Checks
This is the source code for HAMMER checks. The `CMakeLists.txt` is directly lifted from the demo/Examples directory HAMMER produces upon install. Include paths and the like will need to be change to run this in another machine.

Scripts are for each mode, FF parameterisation. Accessing evaluating and accessing the FFs requires a bit of a work-around since `evalAtPSPoint` is a protected member function for the HAMMER FFs. This is achieved in the scripts with a simple wrapper class.

## Instructions
These should be run such that they output to a `.txt` file that will be read in for cross-checks. Compile first, then run appropiate decay mode. From this directory:
```
cmake .
make
./BdToDstFFBGL > ../Hammer_BdToDstFFBGL.txt # or some other executable
```