# Code for simulating the Z_2 Gross-Neveu model

The code currently simulate the Z_2 globally invariant Gross-Neveu model in (2+1) dimensions.  
To compile the code on a local machine, execute the command

> gcc -O3 *.c -o GNModel -lm

Some compilers may have troubles with the cronotime.* files. You may delete those files without problems.
When you run the code, 7 parameters are needed as input parameters

1. double : beta -> plaquette coupling
2. double : beta_scalar -> variance of the sigma field coupling
3. double : mass -> mass of the fermion field
4. double : epsilon_link -> infinitesimal time step for the link HMC
5. double : epsilon_sigma -> infinitesimal time step for the sigma HMC
6. double : epsilon_pi -> infinitesimal time step for the pi HMC
7. int : step_link -> number of steps of the HMC trajectory for the link variables
8. int : step_sigma -> number of steps of the HMC trajectory for the sigma field
9. int : step_pi -> number of steps of the HMC trajectory for the pi field
10. int : seed -> random number seed

For example, to run the GN_Z2_model with 0 mass, beta_scalar=0.85, 2 steps in the HMC and epsilon_sigma=0.2 (if **#define AUTOREGULATION** has been activated, note that this value will be changed) run

> ./GNModel 0.0 0.85 0.0 0.0 0.2 0.0 0 2 0 1231

We strongly suggest to run always the code with **#define START_X0_NULL** instead to **#define START_X0_PREVIOUS** to ensure detailed balance of the HMC algorithm.
The lattice hamiltonian and observables are discussed [here](pdf/Gross_Neveu_algorithm.pdf).
The main parameters required for the simulation are collected in *macro.h*, and they are discussed [here](pdf/discussionmacro.md).