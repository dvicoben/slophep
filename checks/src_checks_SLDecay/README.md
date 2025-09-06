## SL_Decay Checks

This is the source code for SL_decay checks (SL_Decay scripts lifted from https://gitlab.cern.ch/scali/SL_Decay, added here for ease of reproducibility). Note that a bug in the BGL outer functions was found and is fixed in the scripts here (See [MR!2](https://gitlab.cern.ch/scali/SL_Decay/-/merge_requests/2)).

## Instructions
These should be run such that they output to a `.txt` file that will be read in for cross-checks. Compile first, from this directory, and then run appropiate decay mode. For example
```
make BdToDstFFBGL
make clean
./BdToDstFFBGL > ../checks_SLDecay_BdToDstFFBGL.txt 
```
