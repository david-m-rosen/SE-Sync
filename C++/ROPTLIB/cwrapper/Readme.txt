This directory includes C wrappers of BLAS and LAPACK from Fortran.
The code is witten by the group of David Sovboda (http://www.fi.muni.cz/~xsvobod2/).
It can be downloaded from: 
http://www.fi.muni.cz/~xsvobod2/misc/lapack/download/headers.tar.gz

We make slightly modifications in the wrappers:
In blas/f2c.h and lapack/f2c.h files
*) line 10: long int has been changed to int, since in MAC and Ubuntu, the integer of fortain 
	is int in C++, not long int in C++.
*) lines 159 and 160 about max and min are comment out since they
	have been defined in cmath.


Wen Huang