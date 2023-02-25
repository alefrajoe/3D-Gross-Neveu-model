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

To run the code with $\beta_s=0.85$, and a HMC update with two-step update whose length is **step_sigma=0.2** use 

> ./GNModel 0.0 0.85 0.0 0.0 0.2 0.0 0 2 0 1231

We strongly suggest to run always the code with **#define START_X0_NULL** instead to **#define START_X0_PREVIOUS** to ensure detailed balance of the HMC algorithm.
The lattice hamiltonian reads as follows

$$ H = \sum_{x, y} \bar{\psi}_x M_{x, y} \psi_y + \frac{N_f \beta_s}{4}\sum_{\widetilde{x}}\sigma^2_{\widetilde{x}}$$

where the Dirac matrix $M_{x,y}$ is defined as follows

$$M_{x,y}=\sum_{\mu=1}^3\frac{\eta_\mu(x)}{2}\[\delta_{y,x+\mu}-\delta_{y,x-\mu}\] + \frac{1}{8}\sum_{\langle x,\widetilde{x}\rangle}\sigma_{\widetilde{x}}$$