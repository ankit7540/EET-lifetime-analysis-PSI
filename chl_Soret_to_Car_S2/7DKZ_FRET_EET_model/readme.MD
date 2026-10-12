# Chromophore FRET analysis


A set of programs in python for analyzing coordinates from pdb file to look at distance and orientation factor for FRET analysis. This is designed as an example case and is not generalized. 

The associated dataset and code is for crystal structure in 7DKZ.

See `note.txt` for description of specific contents (very useful to understand how the data is arranged in this modelling)

-----------


### Usage 
Download zip, and use under a python environment

### Required packages
numpy, scipy, matplotlib, re, pandas

### software versions
matplotlib                    3.5.3
matplotlib-inline             0.1.6
numpy                         1.21.5
numpydoc                      1.4.0
pandas                        1.3.5
scipy                         1.7.3 
python   	                  3.7.1

## output
Resulting computed output (for distance and orientation factors) are returned as numpy arrays with sorting as per chain ID.

FRET output uses pandas dataframes, and uses distance and orientation information as dataframes. See associated jupyter notebook.


