#ifndef SIMULATION_H
#define SIMULATION_H

#include "lattice.h"
#include "moleculardynamics.h"
#include "condensate.h"
#include "info.h"
#include "inputoutput.h"

/**=================================================================
 * Struct Simulation.
 * The structure contains all the elements necessary for a simulation.
 * ================================================================*/
typedef struct Simulation
{
    struct Lattice * lattice;
    struct Info info;
    struct Fermioncoor * fcor;
    struct Fermionkit * kit;
    struct Force * force;
    struct Momenta * momenta;
    struct Momentasigma * momentasigma;
    struct Forcesigma * forcesigma;
    struct Momentapi * momentapi;
    struct Forcepi * forcepi;

    char observables_file[200];
    char acceptance_file[200];
    char terma_file[200];
    char terma_file_copy[200];
}Simulation;

// simulation functions
void AllocateSimulation(struct Simulation * simulation);
void StartSimulation(struct Simulation * simulation, double beta, double beta_scalar, double mass, double epsilon_link_HMC, double epsilon_sigma_HMC, double epsilon_pi_HMC, int step_link_HMC, int step_sigma_HMC, int step_pi_HMC, int sosia);
void DestroySimulation(struct Simulation * simulation);

// read and write lattice
double TestWriteReadLattice(const struct Simulation const * simulation);
void WriteBinaryLattice(const struct Simulation const * simulation);
void ReadBinaryLattice(struct Simulation * simulation, const char * terma_file, bool try_copy);
#endif