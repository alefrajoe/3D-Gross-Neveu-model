#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "simulation.h"
#include "monopole.h"

/**==================================================
 * Structure Observables.
 * This contains the observables to be printed out in the appropriate file.
 * =================================================*/
typedef struct Observables
{
    #ifdef OBS_MONOPOLE
    double * monopole;
    #endif

    #ifdef OBS_CONDENSATE
    double * condensate_tr1;
    double * condensate_tr1_abs;
    double * condensate_tr2;
    double * condensate_tr4;
    double * condensate_mu2;
    double * condensate_mu4;
    double * condensate_chi;
    double * condensate_chi_abs;
    double * condensate_chip;
    double * condensate12_nosigma;
    double * condensate12_withsigma;
    #endif
    
    #ifdef OBS_CONDENSATE_SIGMA
    double * condensate_sigma;
    double * condensate_sigma_abs;
    double * condensate_sigma_mu2;
    double * condensate_sigma_mu4;
    double * condensate_sigma_chi;
    double * condensate_sigma_chip;
    double * condensate_sigma_chip_pbc1;
    double * condensate_sigma_chip_pbc2;
    #endif

    #ifdef OBS_PLAQUETTE
    double * plaq1;
    double * plaq2;
    double * plaq3;
    #endif
}Observables;

// observables function
void AllocateObservables(struct Observables * obs);
void DestroyObservables(struct Observables * obs);

// take measures
void UpdateLatticeTermalization(struct Simulation * simulation);
void UpdateLattice(struct Simulation * simulation);
void TakeMeasures(struct Simulation * simulation, struct Observables * obs);

// observables functions
void CondensateSigmaPi(const struct Lattice * lat, const struct Info const * info, double * condensate_sigma, double * condensate_sigma_abs, double * mu2, double * mu4, double * chi, double * chip, double * chip_pbc1, double * chip_pbc2);


#endif