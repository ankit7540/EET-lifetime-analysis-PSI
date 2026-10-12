# Chromophore TDC analysis


A set of programs in python for obtaining EET rates and lifetimes from V(coulombic) terms computed from prior ab initio calculation.

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
Resulting computed data is in a pandas dataframe with row and col indices already arranged as per the residues.

The `example_*.ipynb` files are supposed to serve as exmaples for EET rates/lifetimes from earlier V(coulombic) calc, and output saved as a specific arrangement, where donors are across rows, and acceptors are across cols. 


