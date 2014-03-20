integralbinfile
===============

psi4 plugin to generate a binary file containing all the integrals

run "psi4 input.dat" and get "integral.bin"

read "integral.bin" using "read_integral.m" in Matlab

binary file format: 

1 int: number of atoms 

3 ints: nirrep, nrow, ncol for the overlap matrix (S), both nrow and ncol are equal to the number of basis functions (nbasis) 

nbasis*nbasis doubles: integral values in S 

3 ints: nirrep, nrow, ncol for the kinetic energy integral matrix (KE) 

nbasis*nbasis doubles: integral values in KE

% potential energy integral matrices; there are (number of atoms) matricies, each in the form: 

3 ints: nirrep, nrow, ncol for the ith potential energy matrix EN(:, :, i) 

nbasis*nbasis doubles: integral values in EN(:, :, i) 

% 2-electron integrals; there are (number of unique 2-e integrals) of them, each in the form: 

4 doubles: indices i, j, k, l 

1 double: integral value H2(i, j, k, l) 

1 int: variable do_env_in controlling whether we compute internal environment potential integrals 

1 int: number of internal environment point charges (num_ptq_in) 

% internal environment potential energy integral matrices; there are (num_ptq_in) matricies, each in the form: 

3 ints: nirrep, nrow, ncol for the ith internal environment potential energy matrix ENV_IN(:, :, i) 

1 int: variable do_env_ex controlling whether we compute external environment potential integrals 

1 int: number of external environment point charges (num_ptq_ex) 

% external environment potential energy integral matrices; there are (num_ptq_ex) matricies, each in the form: 

3 ints: nirrep, nrow, ncol for the ith external environment potential energy matrix ENV_EX(:, :, i) 

% end of file

