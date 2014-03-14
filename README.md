integralbinfile
===============

psi4 plugin to generate a binary file containing all the integrals

run "psi4 input.dat" and get "integral.bin"

read "integral.bin" using "read_integral.m" in Matlab

binary file format: 

1 int: number of atoms 

3 ints: nirrep, nrow, ncol for overlap matrix (S) 

nbasis*nbasis doubles: S 

3 ints: nirrep, nrow, ncol for KE matrix (KE) 

nbasis*nbasis doubles: KE

potential matrices; there are number of atoms matricies, each has: 

3 ints: nirrep, nrow, ncol for EN(:, :, i) 

nbasis*nbasis doubles: EN(:, :, i) 

2-electron integrals; there are n_unique of them, each has: 

4 ints: indices i, j, k, l 

1 double: integral value H2(i, j, k, l) 

#EOF
