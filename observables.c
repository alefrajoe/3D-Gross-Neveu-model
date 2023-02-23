#include "observables.h"

//////////////////////////////////////////////////////////////////
//                   observables function
//////////////////////////////////////////////////////////////////

/**
 * Allocate the memory for all the observables taken into account in the simulation.
 */
void AllocateObservables(struct Observables * obs)
{
    #ifdef OBS_MONOPOLE
    obs->monopole = malloc(sizeof(double) * IMEAS);
    #endif

    #ifdef OBS_CONDENSATE
    obs->condensate_tr1 = malloc(sizeof(double) * IMEAS);
    obs->condensate_tr1_abs = malloc(sizeof(double) * IMEAS);
    obs->condensate_tr2 = malloc(sizeof(double) * IMEAS);
    obs->condensate_tr4 = malloc(sizeof(double) * IMEAS);
    obs->condensate_mu2 = malloc(sizeof(double) * IMEAS);
    obs->condensate_mu4 = malloc(sizeof(double) * IMEAS);
    obs->condensate_chi = malloc(sizeof(double) * IMEAS);
    obs->condensate_chi_abs = malloc(sizeof(double) * IMEAS);
    obs->condensate_chip = malloc(sizeof(double) * IMEAS);
    obs->condensate12_nosigma = malloc(sizeof(double) * IMEAS);
    obs->condensate12_withsigma = malloc(sizeof(double) * IMEAS);
    #endif

    #ifdef OBS_CONDENSATE_SIGMA
    obs->condensate_sigma = malloc(sizeof(double) * IMEAS);
    obs->condensate_sigma_abs = malloc(sizeof(double) * IMEAS);
    obs->condensate_sigma_mu2 = malloc(sizeof(double) * IMEAS);
    obs->condensate_sigma_mu4 = malloc(sizeof(double) * IMEAS);
    obs->condensate_sigma_chi = malloc(sizeof(double) * IMEAS);
    obs->condensate_sigma_chip = malloc(sizeof(double) * IMEAS);
    obs->condensate_sigma_chip_pbc1 = malloc(sizeof(double) * IMEAS);
    obs->condensate_sigma_chip_pbc2 = malloc(sizeof(double) * IMEAS);
    #endif

    #ifdef OBS_PLAQUETTE
    obs->plaq1 = malloc(sizeof(double) * IMEAS);
    obs->plaq2 = malloc(sizeof(double) * IMEAS);
    obs->plaq3 = malloc(sizeof(double) * IMEAS);
    #endif
}

/**
 * Destroy the allocated memory for the observables.
 */
void DestroyObservables(struct Observables * obs)
{
    #ifdef OBS_MONOPOLE
    free(obs->monopole);
    obs->monopole = NULL;
    #endif

    #ifdef OBS_CONDENSATE
    free(obs->condensate_tr1);
    obs->condensate_tr1 = NULL;
    free(obs->condensate_tr1_abs);
    obs->condensate_tr1_abs = NULL;
    free(obs->condensate_tr2);
    obs->condensate_tr2 = NULL;
    free(obs->condensate_tr4);
    obs->condensate_tr4 = NULL;
    free(obs->condensate_mu2);
    obs->condensate_mu2 = NULL;
    free(obs->condensate_mu4);
    obs->condensate_mu4 = NULL;
    free(obs->condensate_chi);
    obs->condensate_chi = NULL;
    free(obs->condensate_chi_abs);
    obs->condensate_chi_abs = NULL;
    free(obs->condensate_chip);
    obs->condensate_chip = NULL;
    free(obs->condensate12_nosigma);
    obs->condensate12_nosigma = NULL;
    free(obs->condensate12_withsigma);
    obs->condensate12_withsigma = NULL;
    #endif

    #ifdef OBS_CONDENSATE_SIGMA
    free(obs->condensate_sigma);
    obs->condensate_sigma = NULL;
    free(obs->condensate_sigma_abs);
    obs->condensate_sigma_abs = NULL;
    free(obs->condensate_sigma_mu2);
    obs->condensate_sigma_mu2 = NULL;
    free(obs->condensate_sigma_mu4);
    obs->condensate_sigma_mu4 = NULL;
    free(obs->condensate_sigma_chi);
    obs->condensate_sigma_chi = NULL;
    free(obs->condensate_sigma_chip);
    obs->condensate_sigma_chip = NULL;
    free(obs->condensate_sigma_chip_pbc1);
    obs->condensate_sigma_chip_pbc1 = NULL;
    free(obs->condensate_sigma_chip_pbc2);
    obs->condensate_sigma_chip_pbc2 = NULL;
    #endif

    #ifdef OBS_PLAQUETTE
    free(obs->plaq1);
    free(obs->plaq2);
    free(obs->plaq3);
    obs->plaq1 = NULL;
    obs->plaq2 = NULL;
    obs->plaq3 = NULL;
    #endif
}

//////////////////////////////////////////////////////////////////
//                    take measures
//////////////////////////////////////////////////////////////////

/**
 * Update all the sites of the lattice.
 * The update tipology depends on the updates that are active.
 */
void UpdateLatticeTermalization(struct Simulation * simulation)
{    
    for(int i=0; i<TERMALIZATION_STEP; i++)
    {
        #ifdef HYBRID_MC
        // ---------------------------------- Pure gauge -----------------------------------------
        #ifdef PURE_GAUGE
        //while(!LeapFrog2Link(simulation->lattice, &(simulation->info), simulation->force, simulation->momenta, simulation->fcor, simulation->kit, simulation->info.epsilon_link_HMC, simulation->info.steps_link_HMC, simulation->acceptance_file)){printf("Link\n");};
        #endif
        
        
        // ------------------------------ Sigma interaction -----------------------------------------
        #ifdef SIGMA_INTERACTION
        // Use LF2
        #ifdef LF2SIGMA
        while(!LeapFrog2Sigma(simulation->lattice, &(simulation->info), simulation->forcesigma, simulation->momentasigma, simulation->fcor, simulation->kit, simulation->info.epsilon_sigma_HMC, simulation->info.steps_sigma_HMC, simulation->acceptance_file)){printf("Sigma\n");};
        #endif
        // Use LF2MN2
        #ifdef LFMN2SIGMA
        while(!LeapFrog2MinimumNormSigma(simulation->lattice, &(simulation->info), simulation->forcesigma, simulation->momentasigma, simulation->fcor, simulation->kit, simulation->info.epsilon_sigma_HMC, simulation->info.steps_sigma_HMC, simulation->acceptance_file)){printf("Sigma\n");};
        #endif
        // Use LFMN2SIGMAVF
        #ifdef LFMN2SIGMAVF
        while(!LeapFrog2MinimumNormSigmaVelocityFirst(simulation->lattice, &(simulation->info), simulation->forcesigma, simulation->momentasigma, simulation->fcor, simulation->kit, simulation->info.epsilon_sigma_HMC, simulation->info.steps_sigma_HMC, simulation->acceptance_file)){printf("Sigma\n");};
        #endif
        // Use 4MN4FP
        #ifdef LF4MN4FPSIGMA
        printf("Integrator to be implemented!\n");
        exit(1);
        #endif
        #endif


        #ifdef PI_INTERACTION
        while(!LeapFrog2Pi(simulation->lattice, &(simulation->info), simulation->forcepi, simulation->momentapi, simulation->fcor, simulation->kit, simulation->info.epsilon_pi_HMC, simulation->info.steps_pi_HMC, simulation->acceptance_file)){printf("Pi\n");};
        #endif
        #endif

        // update HMC parameters
        UpdateHMCParameters(&(simulation->info));
    }

    // iterations printed out (both stdout to data file).. set to zero again
    simulation->lattice->done_sigma = 0;
    simulation->lattice->distance_sigma = 0;
    simulation->kit->iterations = 0.0;
    simulation->kit->CG_call = 0;

    // update the measure: it will start from run 1 at the next run
    simulation->info.measure++;
    // write the lattice to file
    WriteBinaryLattice(simulation);
}

/**
 * Update all the sites of the lattice.
 * The update tipology depends on the updates that are active.
 */
void UpdateLattice(struct Simulation * simulation)
{    
    #ifdef HYBRID_MC
    // ---------------------------------- Pure gauge -----------------------------------------
    #ifdef PURE_GAUGE
    //while(!LeapFrog2Link(simulation->lattice, &(simulation->info), simulation->force, simulation->momenta, simulation->fcor, simulation->kit, simulation->info.epsilon_link_HMC, simulation->info.steps_link_HMC, simulation->acceptance_file)){printf("Link\n");};
    #endif
    
    
    // ------------------------------ Sigma interaction -----------------------------------------
    #ifdef SIGMA_INTERACTION
    // Use LF2
    #ifdef LF2SIGMA
    while(!LeapFrog2Sigma(simulation->lattice, &(simulation->info), simulation->forcesigma, simulation->momentasigma, simulation->fcor, simulation->kit, simulation->info.epsilon_sigma_HMC, simulation->info.steps_sigma_HMC, simulation->acceptance_file)){fprintf(stderr, "Inversion problem sigma!!!\n"); exit(2);};
    #endif
    // Use LF2MN2
    #ifdef LFMN2SIGMA
    while(!LeapFrog2MinimumNormSigma(simulation->lattice, &(simulation->info), simulation->forcesigma, simulation->momentasigma, simulation->fcor, simulation->kit, simulation->info.epsilon_sigma_HMC, simulation->info.steps_sigma_HMC, simulation->acceptance_file)){fprintf(stderr, "Inversion problem sigma!!!\n"); exit(2);};
    #endif
    // Use LFMN2SIGMAVF
    #ifdef LFMN2SIGMAVF
    while(!LeapFrog2MinimumNormSigmaVelocityFirst(simulation->lattice, &(simulation->info), simulation->forcesigma, simulation->momentasigma, simulation->fcor, simulation->kit, simulation->info.epsilon_sigma_HMC, simulation->info.steps_sigma_HMC, simulation->acceptance_file)){fprintf(stderr, "Inversion problem sigma!!!\n"); exit(2);};
    #endif
    // Use 4MN4FP
    #ifdef LF4MN4FPSIGMA
    printf("Integrator to be implemented!\n");
    exit(1);
    #endif
    #endif


    #ifdef PI_INTERACTION
    while(!LeapFrog2Pi(simulation->lattice, &(simulation->info), simulation->forcepi, simulation->momentapi, simulation->fcor, simulation->kit, simulation->info.epsilon_pi_HMC, simulation->info.steps_pi_HMC, simulation->acceptance_file)){printf("Pi\n");};
    #endif
    #endif
}

/**
 * The function computes and take the observables defined into the macro.h file.
 * The observables are stored into the appropriate array "* obs".
 */
void TakeMeasures(struct Simulation * simulation, struct Observables * obs)
{
    int path = 0;

    #ifdef OBS_CONDENSATE
    // allocate the memory for a random vector field if required
    struct RandomVec randomvec;
    AllocateRandomVec(&randomvec);
    #endif

    // along the path
    for(path=0; path<(MEASURE * IMEAS); path++)
    {     
        // update lattice
        for(int i=0; i<LATTICE_UPDATE; i++)
        UpdateLattice(simulation);

        // update acceptance if a metropolis update has been used
        #ifdef METROPOLIS_LINK
        UpdateMetropolisAcceptance(&(simulation->info), simulation->acceptance_file);
        #endif

        // compute observables
        #ifdef OBS_MONOPOLE
        obs->monopole[path % IMEAS] = MonopoleDensity(simulation->lattice);
        #endif

        #ifdef OBS_CONDENSATE
        ComputeCondensateDirac(simulation->lattice, &(obs->condensate_tr1[path % IMEAS]), &(obs->condensate_tr1_abs[path % IMEAS]),&(obs->condensate_tr2[path % IMEAS]), &(obs->condensate_tr4[path % IMEAS]), &(obs->condensate_mu2[path % IMEAS]), &(obs->condensate_mu4[path % IMEAS]), &(obs->condensate_chi[path % IMEAS]), &(obs->condensate_chi_abs[path % IMEAS]), &(obs->condensate_chip[path % IMEAS]), &(obs->condensate12_nosigma[path % IMEAS]), &(obs->condensate12_withsigma[path % IMEAS]), &randomvec, &(simulation->kit->inverted_phi_obs.flav[0]), &(simulation->info), simulation->fcor, HIGH_PRECISION_INVERSION, simulation->kit, 0);
        #endif

        #ifdef OBS_CONDENSATE_SIGMA
        CondensateSigmaPi(simulation->lattice, &(simulation->info), &(obs->condensate_sigma[path % IMEAS]), &(obs->condensate_sigma_abs[path % IMEAS]), &(obs->condensate_sigma_mu2[path % IMEAS]), &(obs->condensate_sigma_mu4[path % IMEAS]), &(obs->condensate_sigma_chi[path % IMEAS]), &(obs->condensate_sigma_chip[path % IMEAS]), &(obs->condensate_sigma_chip_pbc1[path % IMEAS]), &(obs->condensate_sigma_chip_pbc2[path % IMEAS]));
        #endif

        #ifdef OBS_PLAQUETTE
        ObsPlaquette(simulation->lattice, &(obs->plaq1[path % IMEAS]), &(obs->plaq2[path % IMEAS]), &(obs->plaq3[path % IMEAS]));
        #endif

        // perform the average
        if((path % IMEAS) == (IMEAS - 1))
        {            
            // print the measure number to stdout
            fprintf(stdout, "%d\n", simulation->info.measure);

            FILE * opf=fopen(simulation->observables_file, "a");
            // every time fprintf beta and pointer
            fprintf(opf, "%d\t%d\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f", simulation->info.measure, LSPACE, simulation->info.beta, simulation->info.beta_scalar, simulation->info.mass, simulation->kit->iterations / simulation->kit->CG_call, simulation->lattice->distance_sigma / simulation->lattice->done_sigma);

            #ifdef DEBUG_PRINT_ITERATIONS
            printf("Average iteration: %.16f\n", simulation->kit->iterations / simulation->kit->CG_call);
            #endif

            // iterations printed out (both stdout to data file).. set to zero again
            simulation->kit->iterations = 0.0;
            simulation->kit->CG_call = 0;
            // reset also the distance done in fcor
            simulation->lattice->distance_sigma = 0.0;
            simulation->lattice->done_sigma = 0;

            // write to file all the observables taken
            #ifdef OBS_MONOPOLE
            fprintf(opf, "\t%.16f", Mean(obs->monopole, IMEAS));
            #endif

            #ifdef OBS_CONDENSATE
            fprintf(opf, "\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f", Mean(obs->condensate_tr1, IMEAS), Mean(obs->condensate_tr1_abs, IMEAS), Mean(obs->condensate_tr2, IMEAS), Mean(obs->condensate_tr4, IMEAS), Mean(obs->condensate_mu2, IMEAS), Mean(obs->condensate_mu4, IMEAS), Mean(obs->condensate_chi, IMEAS), Mean(obs->condensate_chi_abs, IMEAS), Mean(obs->condensate_chip, IMEAS), Mean(obs->condensate12_nosigma, IMEAS), Mean(obs->condensate12_withsigma, IMEAS));
            #endif

            #ifdef OBS_CONDENSATE_SIGMA
            fprintf(opf, "\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f\t%.16f", Mean(obs->condensate_sigma, IMEAS), Mean(obs->condensate_sigma_abs, IMEAS), Mean(obs->condensate_sigma_mu2, IMEAS), Mean(obs->condensate_sigma_mu4, IMEAS), Mean(obs->condensate_sigma_chi, IMEAS), Mean(obs->condensate_sigma_chip, IMEAS), Mean(obs->condensate_sigma_chip_pbc1, IMEAS), Mean(obs->condensate_sigma_chip_pbc2, IMEAS));
            #endif

            #ifdef OBS_PLAQUETTE
            fprintf(opf, "\t%.16f\t%.16f\t%.16f", Mean(obs->plaq1, IMEAS), Mean(obs->plaq2, IMEAS), Mean(obs->plaq3, IMEAS));
            #endif
            fprintf(opf, "\n");

            // measure++
            simulation->info.measure++;

            // always write to a binary file
            WriteBinaryLattice(simulation);

            fclose(opf);
        }
    }

    #ifdef OBS_CONDENSATE
    DestroyRandomVec(&randomvec);
    #endif

    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                observables functions
//////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Compute the observables that is related to the sigma-pi condensate.
 * In particular, the function saves in the proper variables the value of the sigma-pi condensate, the susceptibility chi,
 * its fourier transform chip and the useful Binder cumulants mu2 and mu4.
*/
void CondensateSigmaPi(const struct Lattice * lat, const struct Info const * info, double * condensate_sigma, double * condensate_sigma_abs, double * mu2, double * mu4, double * chi, double * chip, double * chip_pbc1, double * chip_pbc2)
{
    double total = 0.0, mu2_box = 0.0;
    complex box[DIM];

    // initialize all variables
    (*condensate_sigma) = 0.0; (*condensate_sigma_abs) = 0.0; (*mu2) = 0.0; (*mu4) = 0.0; (*chi) = 0.0; (*chip) = 0.0;
    // initialize all DIM instances of the box even if just the first DIM-1 are used
    for(int d=0; d<DIM; d++) box[d] = 0.0;

    // for all lattice sites
    // x here is in lessicographic coordinates
    for(int x=0; x<VOLUME; x++)
    {
        // Z2 Gross-Neveu model
        #ifdef SIGMA_INTERACTION
        #ifndef PI_INTERACTION        
        double temp = lat->sigma[lat->fermionsite[x]] / VOLUME;
        #endif
        #endif
        // U(1) Gross-Neveu model
        #ifdef SIGMA_INTERACTION
        #ifdef PI_INTERACTION
        double temp = pow((pow(lat->sigma[lat->fermionsite[x]], 2) + pow(lat->pi[lat->fermionsite[x]], 2)), 0.5) / VOLUME;
        #endif
        #endif

        // compute the condensate_sigma by adding the temp value
        (*condensate_sigma) += temp;
        // to compute the Fourier transform of the order parameter, save the value into the DIM-1 box variables
        // the quantity is averaged only over the temporal direction
        for(int d=0; d<DIM; d++)
        box[d] += cexp(I_UNIT * P_MIN * lat->lattice_site[COOR(x, d)]) * temp;
    }

    (*condensate_sigma_abs) = fabs((*condensate_sigma));
    (*mu2) = pow((*condensate_sigma), 2);
    (*mu4) = pow((*mu2), 2);
    (*chi) = VOLUME * (*mu2);

    // chip is averaged over all the lattice directions (PBC needed)
    // at least 2 directions in the lattice are needed
    //for(int d=0; d<DIM; d++)
    (*chip) = (VOLUME * pow(cabs(box[2]), 2));
    (*chip_pbc1) = (VOLUME * pow(cabs(box[0]), 2));
    (*chip_pbc2) = (VOLUME * pow(cabs(box[1]), 2));
}
