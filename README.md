# Scripts for GEDDs paper

These scripts have been written by Edvinas Pakanavicius for the statistical analysis of GEDDs.


## System requirements
These Python scripts are dependent on the following Python packages:
matplotlib==2.2.3
numpy==1.15.1
statsmodels==0.6.1



These scripts have been tested on the following versions:
Centos 7
Ubuntu 16.04 LTS
Python version 2.7.12

There is no non-standard hardware required to run these scripts.


## Installation guide
There is no special installation required for the analysis scripts - they can be saved in the directory with the data. 
Python packages can be installed using pip.


## Demo
Run the script by typing the following on command line:
python flips.py 1 test_data.txt 
The number '1' indicates the chromosome number. The test data only has data from chr1.
Note that the script is dependent on the order of the columns and the data cannot contain column names.
The script takes a few minutes to run on a standard desktop computer. It outputs a .png file.


## Instructions to use
python flips.py <CHRNUM> <INPUTFILENAME>
python energy.py <CHRNUM> <INPUTFILENAME>
