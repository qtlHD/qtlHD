Xgap Binary Format
==================
Version: 0.0.1

Xgap binary is a format to store genotype phenotype and other information. The format saves matrices of different types into a single binary file.
A block structure separated by footprints is used as easy consistency check of the file, and provides a way to load specific partys of the data from 
file to memory.

Header
------
The header describes the content of the rest of the file, and starts and ends with a footprint:

byte[2] b_footprint = [ 0, 5 ];

The footprint is followed by the version of the binary file written:

byte[3] b_version = [ 0, 0, 1 ];

Then 4 bytes are reserved for the number of matrices inside the file (in little endian notation), followed by a closing footprint. The size of 
the entire header is thus: 2+3+4+2 = 11 bytes.

Matrices
--------