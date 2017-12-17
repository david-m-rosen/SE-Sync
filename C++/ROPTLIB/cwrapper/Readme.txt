This directory includes C wrappers of BLAS and LAPACK from Fortran.
The code is witten by the group of David Sovboda (http://www.fi.muni.cz/~xsvobod2/).
It can be downloaded from: 
http://www.fi.muni.cz/~xsvobod2/misc/lapack/download/headers.tar.gz

We make slightly modifications in the wrappers:
We use only one f2c.h header file for both BLAS and LAPACK
*) line 10: long int has been changed to int, since in MAC and Ubuntu, the integer of fortain 
	is int in C++, not long int in C++.
*) line 14 (and other cases of real): real is renamed to f2c_real in order to
  avoid conflicts with std::real or real(...) in Eigen library
*) lines 158 to 166 are commented out since their definitions
	are not used in this context and provoke conflicts.

Wen Huang
Jesus Briales
