integralbinfile
===============

psi4 plugin to generate a binary file containing all the integrals

run "psi4 input.dat" and get "integral.bin"

read "integral.bin" using "read_integral.m" in Matlab

binary file format: 

1 int: number of atoms 

3 ints: nirrep, nrow, ncol for the overlap matrix (S),

% both nrow and ncol are equal to the number of basis functions (nbasis) 

nbasis*nbasis doubles: integral values in S 

3 ints: nirrep, nrow, ncol for the kinetic energy matrix (KE) 

nbasis*nbasis doubles: integral values in KE

% potential energy matrices; there are (number of atoms) matricies, each in the form: 

3 ints: nirrep, nrow, ncol for the ith potential energy matrix EN(:, :, i) 

nbasis*nbasis doubles: integral values in EN(:, :, i) 

% 2-electron integrals; there are (number of unique 2-e integrals) of them, each in the form: 

4 ints: indices i, j, k, l 

1 double: integral value H2(i, j, k, l) 

1 double: nuclear repulsion energy (e_nuc) 

% end of file

