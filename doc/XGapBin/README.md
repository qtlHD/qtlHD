Xgap Binary Format
==================
Version: 0.0.1

Xgap binary is a format to store genotype phenotype and other information. The format saves matrices of different types into a single binary file.
A block structure separated by footprints is used as easy consistency check of the file, and provides a way to load specific partys of the data from 
file to memory.

Header
------
The header describes the content of the rest of the file, and starts and ends with a footprint, followed by the version of the binary file 
supplied as a byte[3] array (Major, Minor, Build):

 - byte[2] footprint = [ 0, 5 ]
 - byte[3] version = [ 0, 0, 1 ]
 - byte[4] Number of matrices
 - byte[2] footprint = [ 0, 5 ]

Then 4 bytes are reserved for the number of matrices inside the file (in little endian notation), followed by a closing footprint. The size of 
the entire header is thus: 2+3+4+2 = 11 bytes.

Matrices
--------
Xgap bin stores matrices using the XGAP defined data type, the possible types to store data are:

 - EMPTY - No data is stored in this matrix
 - INTMATRIX - Matrix of integers
 - DOUBLEMATRIX - Matrix of doubles
 - FIXEDCHARMATRIX - Matrix of characters with a fixed length (e.g. "AA","AB","BB" Genotype matrices)
 - VARCHARMATRIX - Matrix of any size strings
 
In the XGAP binary format all matrices are surrounded by footprints (see above), followed by the type encoded by 4 bytes. Note: The first 
footprint is not printed, this is because the previous matrix/header already printed a footprint:

 - byte[4] type
 - byte[4] nrow - number of row
 - byte[4] ncol - number of columns
 - byte[4] size || byte[4]*(nrow*ncol) sizes
 - byte[ ] DATA
 - byte[4] footprint = [ 0, 5 ]

The only weird thing here is size or sizes. Size if the length of a single element when we store a INTMATRIX, DOUBLEMATRIX or 
FIXEDCHARMATRIX. When storing a VARCHARMATRIX  we need to know the length of each individual entry in the matrix, and size is 
not a single entry but a list of entries.

Binary Cross Object
-------------------
Converting an existsing cross object (stored as CSV or CSVr) can be done using the csv2xgap utility. This will create an XGAP 
binary file with 3 matrices: the first one storing: phenotype information in a DOUBLEMATRIX, next the genotypes using a FIXEDCHARMATRIX.
Last is the marker map, stored as a VARCHARMATRIX each row containing markername, chromosome, location.
