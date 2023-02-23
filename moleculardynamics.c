#include "moleculardynamics.h"

///////////////////////////////////////////////////////////
//                Momenta functions
///////////////////////////////////////////////////////////
/**
 * The function allocates memory for DIM * VOLUME momenta.
 */
void AllocateMomenta(struct Momenta * momenta)
{
    // allocate memory
    momenta->momenta = malloc(sizeof(double) * DIM * VOLUME);
    // initialize trivially all the momenta to zero
    for(int x=0; x<DIM*VOLUME; x++) momenta->momenta[x] = 0.0;
}

/**
 * Print all the momenta to stdout.
 */
void PrintMomenta(const struct Momenta const * momenta)
{
    printf("************* Momenta ***************\n");
    // for all DIM * VOLUME elements
    for(int r=0; r<DIM * VOLUME; r++)
    {
        printf("Momenta [%d] is %.16f\n", r, momenta->momenta[r]);
    }
    printf("***********************************\n");
}

/**
 * The function draws DIM * VOLUME momenta (VOLUME is always an EVEN number), so that
 * these are DIM * VOLUME elements drawn from a gaussia distribution with zero average and unit variance.
 * **************************************************************************
 * P ( momenta ) \propto  exp ( - ( momenta^2 ) / 2 )
 * **************************************************************************
 */
void DrawMomentaFromGaussianDistribution(struct Momenta * momenta)
{
    // for DIM * VOLUME / 2 steps
    for(int r=0; r<(DIM * VOLUME / 2); r++)
    {
        // draw 2 gaussian numbers with mean zero and unit variance
        Gauss2MeanVariance(&(momenta->momenta[(2*r)]), &(momenta->momenta[(2*r)+1]), 0.0, 1.0);
    }
}

/**
 * Deallocate the memory of the passed momenta.
 * Set the pointer of the passed momenta to NULL.
 */
void DestroyMomenta(struct Momenta * momenta)
{
    free(momenta->momenta);
    momenta->momenta = NULL;
}

/**
 * Return the norm of the force, according to the L_2 norm.
 */
double NormForce(const struct Force const * force)
{
    double norm = 0.0;

    for(int p=0; p<(DIM * VOLUME); p++)
    norm += cabs(force->force[p])*cabs(force->force[p]);

    return pow(norm, 0.5);
}

//////////////////////////////////////////////////////////////////////////////////
//                             force functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * Allocate the memory for the force associated with the gauge fields.
 */
void AllocateForce(struct Force * force)
{
    // allocate memory
    force->force = malloc(sizeof(double) * DIM * VOLUME);
    // initialize trivially all the forces to zero
    for(int x=0; x<DIM*VOLUME; x++) force->force[x] = 0.0;
}

/**
 * Print the force to stdout.
 */
void PrintForce(const struct Force const * force)
{
    printf("************* Force ***************\n");
    // for all DIM * VOLUME elements
    for(int r=0; r<DIM * VOLUME; r++)
    {
        printf("Force [%d] is %.16f\n", r, force->force[r]);
    }
    printf("***********************************\n");
}

/**
 * Set all the elements of the force to zero.
 */
void SetForceToZero(struct Force * force);

/**
 * Deallocate the memory of the passed force.
 * Set the pointer of the passed force to NULL.
 */
void DestroyForce(struct Force * force)
{
    free(force->force);
    force->force = NULL;
}

///////////////////////////////////////////////////////////////////////
//                    force computation
///////////////////////////////////////////////////////////////////////
/**
 * Add the forces of all the DIM * VOLUME coordinates, associated with the pure gauge sector.
 * The force is not initialized to zero.
 * Remember that the coordinates for "struct Momenta" and "struct Force" are given by MOM_COOR(x, d).
 * If Q_\mu ( x ) defines the coordinates conjugated to the respective momenta, we find for the force associated
 * with the pure gauge sector F_G
 * ====================================================================================
 * F_G [ Q_\mu ( x ) ] = -   \beta  \sum_{\nu=1}^{DIM}  Imag { \Pi_{\mu\nu}(x) }
 * ====================================================================================
 */
void AddForceFromPureGaugeSector(const struct Lattice const * lattice, const struct Info const * info, struct Force * force)
{
    // for all DIM * VOLUME elements in force
    // this x running over the whole lattice is in Lessicographic coordinates:
    // save into the force in fermion-Lessicographic coordinate
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        // pass to the function the site in fermionic coordinates
        force->force[COOR(lattice->fermionsite[x], d)] = info->beta * ComputeImagAllPlaquettesIncludingLink(lattice, x, d);
    }
}

/**
 * Auxiliary function required to compute the force associated with the pseudo-fermion component of the action.
 * The function computes the product of the partial derivative of the Dirac operator and a fermion (with both upper and lower components).
 * Note that the target fermion is not initialized, only the relevant vector elements (x and x+mu) are changed.
 * - fermion_x : is the site associated with the derivative in Lessicographic fermion coordinates
 * - mu = is the direction associated with the derivative
 * .. such that the link that has been derived is U_mu(x)
 * --------------------------------------------------------------------------------------------------------------
 * \partial D^\dagger / \partial \theta_mu(x)   =  - (0.5) i \eta_mu(x) [ U^*_mu(x) \delta_{x+mu, x} + U_mu(x) \delta_{x, x+mu}]
 * --------------------------------------------------------------------------------------------------------------
 * */ 
void PartialAdjDiracFermion(const struct Lattice const * lat, const struct Fermion const * start, struct Fermion * target, int fermion_x, int mu, const struct Fermioncoor const * fermioncoor)
{
    // set the first component of target --> U^*_mu(x)
    target->fermion[fermioncoor->fermionpp[COOR(fermion_x, mu)]] = - (0.5) * I_UNIT * fermioncoor->eta_pp[COOR(fermion_x, mu)] * conj(lat->dlink[COOR(fermioncoor->site[fermion_x], mu)].link) * start->fermion[fermion_x];

    // set the second non vanishing component --> U_mu(x)
    target->fermion[fermion_x] = - (0.5) * I_UNIT * fermioncoor->eta_pp[COOR(fermion_x, mu)] * lat->dlink[COOR(fermioncoor->site[fermion_x], mu)].link * start->fermion[fermioncoor->fermionpp[COOR(fermion_x, mu)]];
}

/**
 * The function add the forces associated with the pure fermionic sector to "force" (that hopefully has already been initialized to zero).
 * We remember that the "x" passed to the function should be in FERMIONIC Lessicographic coordinates.
 */
bool AddForceFromFermionSectorRespectToLink(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Force * force, struct Fermionkit * kit)
{
    // first invert the Phi field
    // both the Phi field and the inverted_phi field are located into the fermionkit
    // invert all flavors
    for(int f=0; f<HALFNFLAV; f++)
    {
        if(!ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fcor, HIGH_PRECISION_INVERSION, kit, f)) return false;
    }

    // for all flavors
    // compute (D^\dagger)_{ij}( F_inv )_j
    for(int f=0; f<HALFNFLAV; f++)
    AdjDiracFermion(lat, &(kit->inverted_phi.flav[f]), &(kit->inverted_phi_obs.flav[f]), info, fcor);

    // at this stage we have computed (D D^\dagger)^{-1} (\phi)
    // compute all components of the force
    for(int f=0; f<HALFNFLAV; f++)
    for(int x=0; x<VOLUME; x++)
    for(int mu=0; mu<DIM; mu++)
    {
        //
        force->force[COOR(x, mu)] += creal(conj(kit->inverted_phi.flav[f].fermion[x]) * I_UNIT * fcor->eta_pp[COOR(x, mu)] * lat->dlink[COOR(fcor->site[x], mu)].link * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(x, mu)]]);
        //
        force->force[COOR(x, mu)] += creal(conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(x, mu)]]) * I_UNIT * fcor->eta_mm[COOR(fcor->fermionpp[COOR(x, mu)], mu)] * conj(lat->dlink[COOR(fcor->site[x], mu)].link) * kit->inverted_phi_obs.flav[f].fermion[x]);
    }

    // return if all the inversions has been carried out correctly
    return true;
}

/**
 * Compute the global force required for the updates of the momenta.
 * The force is initially set to zero.
 * The elements of the force vectors are saved in Lessicographic fermionic coordinates.
 * Compute force is composed of two main terms:
 *                                                - F_G  :  from pure gauge sector
 *                                                - F_PF :  from the pseudo-fermionic gauge-sector
 */
bool ComputeForce(const struct Lattice const * lattice, const struct Info const * info, struct Force * force, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit)
{
    // first initialize the forces
    SetForceToZero(force);

    // -------------------- PURE GAUGE SECTOR --------------------------
    // add the forces of 
    #ifdef PURE_GAUGE
    AddForceFromPureGaugeSector(lattice, info, force);
    #endif

    // ------------------- PSEUDO FERMIONIC SECTOR -------------------
    // before evolving in time, compute phi and update the Dirac matrix
    #ifdef PSEUDO_FERMION
    // Add the force from fermionic sector
    if(!AddForceFromFermionSectorRespectToLink(lattice, info, fermioncoor, force, kit)) return false;
    #endif

    return true;
}

////////////////////////////////////////////////////////////////////////
//                   energy functions
////////////////////////////////////////////////////////////////////////
/**
 * Compute the total fictitious energy associated with the model.
 * Recall that the lattice model fictitious hamiltonian is made of three different terms.
 * The fermionic component of the ficitious energy is computed by means of the field \Chi.
 *                                                      
 * ===============================================================
 * H =  (P^2) / 2     +       S_G        +        S_PF
 * ===============================================================
 */
double ComputeTotalFictitiousEnergyWithChi(const struct Lattice const * lattice, const struct Info const * info, const struct Momenta const * momenta, const struct Fermionkit const * kit)
{
    // initialize the variable that stores the energy
    double total_energy = 0.0;

    // for all momenta elements
    // add      (P^2) / 2
    for(int r=0; r<(DIM*VOLUME); r++) total_energy += pow(momenta->momenta[r], 2.0) / 2.0;

    // add the pure gauge contribution
    #ifdef PURE_GAUGE
    total_energy -= (info->beta * TotalRealTracePlaquette(lattice));
    #endif

    #ifdef PSEUDO_FERMION
    // add the pure fermionic contribution
    // in the case of the chi field the calcualation is extremely simple
    // sum over the non-null HALFVOLUME components
    // S_PF  =  ( \chi^* )_i  ( \chi )_i
    for(int f=0; f<HALFNFLAV; f++)
    for(int x=0; x<VOLUME; x++)
    total_energy += pow(cabs(kit->chi.flav[f].fermion[x]), 2.0);
    #endif

    return total_energy;
}

/**
 * Compute the total upper energy associated with the model.
 * Recall that the lattice model fictitious hamiltonian is made of three different terms.
 * The fermionic component of the ficitious energy is computed by means of the field 
 *                                                      \Phi = ( D )_ij  \Chi_j
 * ===============================================================
 * H =  (P^2) / 2     +       S_G        +        S_PF
 * ===============================================================
 */
double ComputeFullEnergyWithPhiWithoutSigma(const struct Lattice const * lattice, const struct Info const * info, const struct Fermioncoor const * fcor, const struct Momenta const * momenta, const struct Fermionkit const * kit)
{
    // initialize the variable that stores the energy
    double total_energy = 0.0;

    // for all momenta elements
    // add      (P^2) / 2
    for(int r=0; r<(DIM*VOLUME); r++) total_energy += pow(momenta->momenta[r], 2.0) / 2.0;

    // add the pure gauge contribution
    #ifdef PURE_GAUGE
    total_energy -= (info->beta * TotalRealTracePlaquette(lattice));
    #endif

    #ifdef PSEUDO_FERMION
    // for all the flavors, add the energy contribution
    for(int f=0; f<HALFNFLAV; f++)
    for(int x=0; x<VOLUME; x++)
    total_energy += creal(conj(kit->phi.flav[f].fermion[x]) * kit->inverted_phi.flav[f].fermion[x]);
    #endif

    return total_energy;
}

/**
 * Compute the total fictitious energy required for the update of HMC of the sigma fields.
 */
double ComputeFullEnergyForSigma(const struct Lattice const * lattice, const struct Info const * info, const struct Fermioncoor const * fcor, const struct Momentasigma const * momenta, const struct Fermionkit const * kit)
{
    double total_energy = 0.0;

    // for all sigma momenta elements
    // add      (P^2) / 2
    for(int r=0; r<VOLUME; r++) total_energy += (pow(momenta->momenta[r], 2) / 2.0);

    #ifdef PSEUDO_FERMION
    // for all the flavors, add the energy contribution
    for(int f=0; f<HALFNFLAV; f++)
    for(int x=0; x<VOLUME; x++)
    total_energy += ((0.5) * creal(conj(kit->phi.flav[f].fermion[x]) * kit->inverted_phi.flav[f].fermion[x]));
    #endif

    #ifdef SIGMA_INTERACTION
    total_energy += ((NFLAV * info->beta_scalar / 4.0) * ComputeTotalSigmaSquared(lattice));
    #endif

    return total_energy;
}

/**
 * Compute the total fictitious energy required for the update of HMC of the sigma fields with the chi field.
 */
double ComputeFullEnergyForSigmaWithChi(const struct Lattice const * lattice, const struct Info const * info, const struct Fermioncoor const * fcor, const struct Momentasigma const * momenta, const struct Fermionkit const * kit)
{
    double total_energy = 0.0;

    // for all sigma momenta elements
    // add      (P^2) / 2
    for(int r=0; r<VOLUME; r++) total_energy += (pow(momenta->momenta[r], 2) / 2.0);

    #ifdef PSEUDO_FERMION
    // for all the flavors, add the energy contribution
    for(int f=0; f<HALFNFLAV; f++)
    for(int x=0; x<VOLUME; x++)
    total_energy += ((0.5) * creal(conj(kit->chi.flav[f].fermion[x]) * kit->chi.flav[f].fermion[x]));
    #endif

    #ifdef SIGMA_INTERACTION
    total_energy += ((NFLAV * info->beta_scalar / 4.0) * ComputeTotalSigmaSquared(lattice));
    #endif

    return total_energy;
}

/**
 * Compute the full energy associated with the molecular dynamics of the Pi field.
 */
double ComputeFullEnergyForPi(const struct Lattice const * lattice, const struct Info const * info, const struct Fermioncoor const * fcor, const struct Momentapi const * momenta, const struct Fermionkit const * kit)
{
    double total_energy = 0.0;

    // for all sigma momenta elements
    // add      (P^2) / 2
    for(int r=0; r<VOLUME; r++) total_energy += (pow(momenta->momenta[r], 2.0) / 2.0);

    #ifdef PSEUDO_FERMION
    // for all the flavors, add the energy contribution
    for(int f=0; f<HALFNFLAV; f++)
    for(int x=0; x<VOLUME; x++)
    total_energy += ((0.5) * creal(conj(kit->phi.flav[f].fermion[x]) * kit->inverted_phi.flav[f].fermion[x]));
    #endif

    #ifdef PI_INTERACTION
    total_energy += ((NFLAV * info->beta_scalar / 8.0) * ComputeTotalPiSquared(lattice));
    #endif

    return total_energy;
}

///////////////////////////////////////////////////////////////////////
//                evolution functions
////////////////////////////////////////////////////////////////////////
/**
 * Evolve all the coordinates for a fictitious amount of time equal to epsilon (this is passed to the function).
 * The hamilton equations related to the evolution of the coordinates follows.
 * ----------------------------------------------------------------------------
 * dQ / dt = dH / dp = P            =>            Q(\epsilon) = Q(0) + P \epsilon
 * ----------------------------------------------------------------------------
 */
void EvolveCoordinates(struct Lattice * lattice, const struct Fermioncoor const * fcor, const struct Momenta const * momenta, double epsilon)
{
    // for all sites in the lattice
    // for all the directions in the lattice
    // x is in fermionic coordinates!!!
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        lattice->dlink[COOR(x, d)].link = cexp(I_UNIT * (carg(lattice->dlink[COOR(x, d)].link) + Q_CHARGE * epsilon * momenta->momenta[COOR(lattice->fermionsite[x], d)])); 
    }
}

/**
 * Evolve all the momenta for a fictitious amount of time equal to epsilon (this is passed to the function).
 * The hamilton equations related to the evolution of the momenta follows.
 * -------------------------------------------------------------------------------------------------------
 * dP / dt = - dH / dq = F_G + F_PF            =>            P(\epsilon) = P(0) + (F_G + F_PF) \epsilon
 * -------------------------------------------------------------------------------------------------------
 */
void EvolveMomenta(const struct Lattice const * lattice, const struct Force const * force, struct Momenta * momenta, double epsilon)
{
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        momenta->momenta[COOR(x, d)] += Q_CHARGE * (force->force[COOR(x, d)] * epsilon);
    }
}

/**
 * 2^nd order LeapFrog (LF2).
 * The function initially draw the momenta (randomly from a gaussian distribution), then evolve globally the link variables.
 * The fundamental infinitesimal step in the LF2 algorithm is given by the following three steps:
 *                      - P(t)                 ->               P( t  +  \epsilon / 2  )
 *                      - Q(t)                 ->               Q( t  +  \epsilon  )
 *                      - P(t+\epsilon/2)      ->               P( t  +   \epsilon  )
 * These steps are repeated an integer amount of times equal to "steps", which is passed to the function.
 * The function return TRUE if the LeapFrog2Link has been carried out correctly.
 */
bool LeapFrog2Link(struct Lattice * lattice, struct Info * info, struct Force * force, struct Momenta * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file)
{
    #ifdef HYBRID_MC
    // file where we will write the Leapfrog acceptance parameters
    FILE * opf = fopen(acceptance_file, "a");
    // first draw initial momenta
    DrawMomentaFromGaussianDistribution(momenta);
    // draw also a random fermionic field to be associated with the chi field
    // for all flavors in the euclidean action
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    {
        DrawFullFermionFromGaussianDistribution(&(kit->chi.flav[f]));
        // compute phi at the beginning of the leapfrog
        DiracFermion(lattice, &(kit->chi.flav[f]), &(kit->phi.flav[f]), info, fermioncoor);
    }
    #endif
    
    // register a copy of all the links
    complex link_copy[DIM * VOLUME];
    CopyAllLinkFromLattice(lattice, link_copy);

    // we have all needed to compute the initial fictitious energy
    // invert all fermions
    // if the inversion is false, return FALSE
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    if(!ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f))
    {
        CopyAllLinkToLattice(link_copy, lattice);
        return false;
    }
    #endif
    
    // compute the starting energy of the molecular dynamics
    double starting_energy = ComputeFullEnergyWithPhiWithoutSigma(lattice, info, fermioncoor, momenta, kit);
    fprintf(opf, "# Slink\t%.16f\n", starting_energy);

    //           P(t)        ->      P( t  +  \epsilon / 2  )
    // first compute the force, then evolve in time
    if(!ComputeForce(lattice, info, force, fermioncoor, kit))
    {
        CopyAllLinkToLattice(link_copy, lattice);
        return false;
    }
    EvolveMomenta(lattice, force, momenta, epsilon / 2.0);

    // for the number of steps passed to the function
    for(int s=0; s<step-1; s++)
    {        
        //          Q(t)         ->      Q( t  +  \epsilon  )
        EvolveCoordinates(lattice, fermioncoor, momenta, epsilon);

        //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
        // first, compute again the force, then evolve in time
        if(!ComputeForce(lattice, info, force, fermioncoor, kit))
        {
            CopyAllLinkToLattice(link_copy, lattice);
            return false;
        }
        EvolveMomenta(lattice, force, momenta, epsilon);
    }

    //          Q(t)         ->      Q( t  +  \epsilon  )
    EvolveCoordinates(lattice, fermioncoor, momenta, epsilon);

    //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
    // first, compute again the force, then evolve in time
    if(!ComputeForce(lattice, info, force, fermioncoor, kit))
    {
        CopyAllLinkToLattice(link_copy, lattice);
        return false;
    }
    EvolveMomenta(lattice, force, momenta, epsilon / 2.0);

    // compute the final fictitious energy
    double final_energy = ComputeFullEnergyWithPhiWithoutSigma(lattice, info, fermioncoor, momenta, kit);
        
    // update delta_H_link in info
    info->delta_H_link = (final_energy - starting_energy);
    
    // accept - reject test
    if(Casuale() <= exp( - (final_energy - starting_energy) ) )
    {
        // the lattice has already been updated
        info->hybrid_accepted_link++;
        info->hybrid_done_link++;
        info->hybrid_done_acceptance_link = (double) (info->hybrid_accepted_link) / (double) (info->hybrid_done_link);

        // update story link
        // if the story is in a non-negative trend, add a positive event
        // else set the story to -1
        if(info->story_link >= 0) info->story_link++;
        else info->story_link = -1;
    }
    else
    {
        // if rejected copy back the all links to the lattice
        CopyAllLinkToLattice(link_copy, lattice);

        // update info
        info->hybrid_done_link++;
        info->hybrid_done_acceptance_link = (double) (info->hybrid_accepted_link) / (double) (info->hybrid_done_link);

        // update the story link
        // if the story is in a positive trend, set the trend to -1
        // else add a negative event
        if(info->story_link >= 0) info->story_link = -1;
        else info->story_link--;
    }

    // ifndef SIGMA_INTERACTION print 
    #ifndef SIGMA_INTERACTION
    fprintf(opf, "%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\n", info->beta, info->beta_scalar, info->mass, info->epsilon_link_HMC, info->steps_link_HMC, (info->epsilon_link_HMC * info->steps_link_HMC), info->delta_H_link, info->hybrid_done_acceptance_link, info->epsilon_sigma_HMC, info->steps_sigma_HMC, (info->epsilon_sigma_HMC * info->steps_sigma_HMC), final_energy-starting_energy, info->hybrid_done_acceptance_sigma);
    #endif
    // print the acceptance to acceptance file
    
    // close the acceptance file
    fclose(opf);

    // HMC correctly done, return TRUE
    return true;
    #endif
}

////////////////////////////////////////////////////////////////////
//                  force functions (sigma)
////////////////////////////////////////////////////////////////////

/**
 * Allocate the memory for the force of the sigma field.
 */
void AllocateForceSigma(struct Forcesigma * force)
{
    force->force = malloc(sizeof(double) * VOLUME);
}

/**
 * Set the Sigma force to zero.
 */
void SetForceSigmaToZero(struct Forcesigma * force);

/**
 * Destroy the memory of a Force sigma field.
 */
void DestroyForceSigma(struct Forcesigma * force)
{
    free(force->force);
    force->force = NULL;
}

/**
 * Add the force from the sigma sector.
 * The force is added to the pointer force.
 * =====================================================================================
 * F(\tilde{x}) = -(NFLAV * beta_scalar / 2.0) \sigma(\tilde{x})
 * ===================================================================================*/
void AddForceFromSigmaSector(const struct Lattice const * lat, const struct Info const * info, struct Forcesigma * force)
{
    // for all VOLUME sites
    for(int x=0; x<VOLUME; x++)
    force->force[x] -= (NFLAV * info->beta_scalar / 2.0) * lat->sigma[x];
}

/**
 * Add the force associated with the kinetic fermionic sector with respect to the sigma field.
 * The function computes the force for all the VOLUME force elements associated with the dual lattice sites where
 * the sigma field is defined.
 * =================================================================================
 * F(\tilde{x}) = ...
 * ===============================================================================*/
bool AddForceFromFermionSectorRespectToSigma(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Forcesigma * force, struct Fermionkit * kit, double inversion_precision)
{
    // first, for all HALFNFLAV fermion fields, invert the matrix
    for(int f=0; f<HALFNFLAV; f++)
    {
        if(!ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fcor, inversion_precision, kit, f)) return false;
    }

    // for all flavors
    // compute ( D^\dagger )_{ij}( inv_phi )_j
    for(int f=0; f<HALFNFLAV; f++)
    AdjDiracFermion(lat, &(kit->inverted_phi.flav[f]), &(kit->inverted_phi_obs.flav[f]), info, fcor);

    // the phi field has been inverted and stored into inverted_phi
    // EVEN + ODD sites contribution
    // (2.0 / 8.0)(\delta_{x,x} + \delta_{x+0,x+0} + \delta_{x+1,x+1} + \delta_{x+2,x+2} + \delta_{x+0+1,x+0+1} + \delta_{x+0+2,x+0+2} + \delta_{x+1+2, x+1+2} \delta_{x+0+1+2,x+0+1+2})(m+sigmamass[x of delta])
    for(int f=0; f<HALFNFLAV; f++)
    for(int x=0; x<VOLUME; x++)
    {
        
        // sum the 4 terms associated with the even sites
        force->force[x] += (0.125) * creal(conj(kit->inverted_phi.flav[f].fermion[x]) * kit->inverted_phi_obs.flav[f].fermion[x]);
        force->force[x] += (0.125) * creal(conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)]]);
        force->force[x] += (0.125) * creal(conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 2)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 2)]]);
        force->force[x] += (0.125) * creal(conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 1)], 2)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 1)], 2)]]);
        // sum the 4 terms associated with the odd sites
        force->force[x] += (0.125) * creal(conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(x, 0)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(x, 0)]]);
        force->force[x] += (0.125) * creal(conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(x, 1)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(x, 1)]]);
        force->force[x] += (0.125) * creal(conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(x, 2)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(x, 2)]]);
        force->force[x] += (0.125) * creal(conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)], 2)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)], 2)]]);
    }
    return true;
}

/**
 * Compute the total force of the sigma field, which is defined on the dual lattice.
 * This function computes the contribution stemming from the Dirac matrix.
 * This contains a \sigma(\tilde{x}) term and the \sigma^2 term.
 * ====================================================================================
 * F[\sigma(\tilde{x})] = F_{\sigma^2} + F_{Dirac}
 * ===================================================================================*/
bool ComputeSigmaForce(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Forcesigma * force, struct Fermionkit * kit, double inversion_precision)
{
    // initialize the force to zero
    SetForceSigmaToZero(force);

    // ifdef the fermion determinant
    #ifdef PSEUDO_FERMION
    if(!AddForceFromFermionSectorRespectToSigma(lat, info, fcor, force, kit, inversion_precision)) return false;
    #endif
    // ifdef SIGMA_INTERACTION
    #ifdef SIGMA_INTERACTION
    AddForceFromSigmaSector(lat, info, force);
    #endif

    // force computed correctly, return TRUE
    return true;
}

////////////////////////////////////////////////////////////////////
//                   force functions (pi)
////////////////////////////////////////////////////////////////////

/**
 * Allocate the memory for the force of the pi field.
 */
void AllocateForcePi(struct Forcepi * force)
{
    force->force = malloc(sizeof(double) * VOLUME);
}

/**
 * Set the Pi force to zero.
 */
void SetForcePiToZero(struct Forcepi * force);

/**
 * Destroy the memory of a pi force field.
 */
void DestroyForcePi(struct Forcepi * force)
{
    free(force->force);
    force->force = NULL;
}

/**
 * Add the force from the pi sector.
 * The force is added to the pointer force.
 * =====================================================================================
 * F(\tilde{x}) = -(NFLAV * beta_scalar / 4.0) \pi(\tilde{x})
 * ===================================================================================*/
void AddForceFromPiSector(const struct Lattice const * lat, const struct Info const * info, struct Forcepi * force)
{
    // for all VOLUME sites
    for(int x=0; x<VOLUME; x++)
    force->force[x] -= (NFLAV * info->beta_scalar / 4.0) * lat->pi[x];
}

/**
 * Add the force associated with the kinetic fermionic sector with respect to the sigma field.
 * The function computes the force for all the VOLUME force elements associated with the dual lattice sites where
 * the sigma field is defined.
 * =================================================================================
 * F(\tilde{x}) = ...
 * ===============================================================================*/
bool AddForceFromFermionSectorRespectToPi(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Forcepi * force, struct Fermionkit * kit)
{
    // first, for all HALFNFLAV fermion fields, invert the matrix
    for(int f=0; f<HALFNFLAV; f++)
    {
        if(!ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fcor, HIGH_PRECISION_INVERSION, kit, f)) return false;
    }

    // for all flavors
    // compute ( D^\dagger )_{ij}( inv_phi )_j
    for(int f=0; f<HALFNFLAV; f++)
    AdjDiracFermion(lat, &(kit->inverted_phi.flav[f]), &(kit->inverted_phi_obs.flav[f]), info, fcor);

    // the phi field has been inverted and stored into inverted_phi
    // EVEN + ODD sites contribution
    // (2.0 / 8.0)(\delta_{x,x} + \delta_{x+0,x+0} + \delta_{x+1,x+1} + \delta_{x+2,x+2} + \delta_{x+0+1,x+0+1} + \delta_{x+0+2,x+0+2} + \delta_{x+1+2, x+1+2} \delta_{x+0+1+2,x+0+1+2})(m+sigmamass[x of delta])
    for(int f=0; f<HALFNFLAV; f++)
    for(int x=0; x<VOLUME; x++)
    {
        
        // sum the 4 terms associated with the even sites
        force->force[x] += (0.125) * creal(I_UNIT * lat->epsilon_pi[x] * conj(kit->inverted_phi.flav[f].fermion[x]) * kit->inverted_phi_obs.flav[f].fermion[x]);
        force->force[x] += (0.125) * creal(I_UNIT * lat->epsilon_pi[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)]] * conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)]]);
        force->force[x] += (0.125) * creal(I_UNIT * lat->epsilon_pi[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 2)]] * conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 2)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 2)]]);
        force->force[x] += (0.125) * creal(I_UNIT * lat->epsilon_pi[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 1)], 2)]] * conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 1)], 2)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 1)], 2)]]);
        // sum the 4 terms associated with the odd sites
        force->force[x] += (0.125) * creal(I_UNIT * lat->epsilon_pi[fcor->fermionpp[COOR(x, 0)]] * conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(x, 0)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(x, 0)]]);
        force->force[x] += (0.125) * creal(I_UNIT * lat->epsilon_pi[fcor->fermionpp[COOR(x, 1)]] * conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(x, 1)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(x, 1)]]);
        force->force[x] += (0.125) * creal(I_UNIT * lat->epsilon_pi[fcor->fermionpp[COOR(x, 2)]] * conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(x, 2)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(x, 2)]]);
        force->force[x] += (0.125) * creal(I_UNIT * lat->epsilon_pi[fcor->fermionpp[COOR(fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)], 2)]] * conj(kit->inverted_phi.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)], 2)]]) * kit->inverted_phi_obs.flav[f].fermion[fcor->fermionpp[COOR(fcor->fermionpp[COOR(fcor->fermionpp[COOR(x, 0)], 1)], 2)]]);
    }
    return true;
}

/**
 * Compute the total force of the pi field, which is defined on the dual lattice.
 * This function computes the contribution stemming from the Dirac matrix and the one that is quadratic in pi.
 * The dirac matrix contains a term proportional to \pi(\tilde{x}).
 * ====================================================================================
 * F[\sigma(\tilde{x})] = F_{\sigma^2} + F_{Dirac}
 * ===================================================================================*/
bool ComputePiForce(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Forcepi * force, struct Fermionkit * kit)
{
    // initialize the force to zero
    SetForcePiToZero(force);

    // ifdef the fermion determinant
    #ifdef PSEUDO_FERMION
    if(!AddForceFromFermionSectorRespectToPi(lat, info, fcor, force, kit)) return false;
    #endif
    // ifdef PI_INTERACTION
    #ifdef PI_INTERACTION
    AddForceFromPiSector(lat, info, force);
    #endif

    // force computed correctly, return TRUE
    return true;
}



///////////////////////////////////////////////////////////////////
//                 momenta functions
///////////////////////////////////////////////////////////////////
//                   (sigma)
/**
 * Allocate the memory for the momenta of the sigma field
 */
void AllocateMomentaSigma(struct Momentasigma * momenta)
{
    momenta->momenta = malloc(sizeof(double) * VOLUME);
}

/**
 * Destroy the momenta of the sigma field.
 */
void DestroyMomentaSigma(struct Momentasigma * momenta)
{
    free(momenta->momenta);
    momenta->momenta = NULL;
}

/**
 * Draw all the VOLUME sigma momenta from a gaussian distribution with zero mean and unit variance.
 * The mean and the variance are passed to the function.
 */
void DrawSigmaMomentaFromGaussianDistribution(struct Momentasigma * momenta)
{
    // for all VOLUME sites
    for(int x=0; x<VOLUME; x++)
    momenta->momenta[x] = Gauss1();
}

//                      (pi) 
/**
 * Allocate the memory for the momenta of the pi field
 */
void AllocateMomentaPi(struct Momentapi * momenta)
{
    momenta->momenta = malloc(sizeof(double) * VOLUME);
}

/**
 * Destroy the momenta of the pi field.
 */
void DestroyMomentaPi(struct Momentapi * momenta)
{
    free(momenta->momenta);
    momenta->momenta = NULL;
}

/**
 * Draw all the VOLUME pi momenta from a gaussian distribution with zero mean and unit variance.
 * The mean and the variance are passed to the function.
 */
void DrawPiMomentaFromGaussianDistribution(struct Momentapi * momenta)
{
    // for all VOLUME sites
    for(int x=0; x<HALFVOLUME; x++)
    Gauss2MeanVariance(&(momenta->momenta[x]), &(momenta->momenta[x+HALFVOLUME]), 0.0, 1.0);
}

//////////////////////////////////////////////////////////////
//                evolve sigma functions
//////////////////////////////////////////////////////////////

/**
 * Evolve the sigma field in time of an infinitesimal step epsilon.
 */
void EvolveCoordinatesSigma(struct Lattice * lat, const struct Momentasigma const * momenta, const struct Fermioncoor const * fcor, double epsilon)
{
    // for all VOLUME fermionic sites
    // evolve the sigma field of an infinitesimal step epsilon * momenta[\tilde{x}]
    for(int x=0; x<VOLUME; x++)
    lat->sigma[x] += (epsilon * momenta->momenta[x]);

    // after the sigma field has been evolved, compute the updated effective sigma mass
    ComputeFullSigmaMass(lat, fcor, lat->sigmamass);

    // keep track of the distance effectively done (this must be equal to the one reporeted in the acceptance file)
    (lat->distance_sigma) += epsilon;
}

/**
 * Evolve all the momenta of the sigma fields of an infinitesimal time step epsilon.
 * ==================================================================================
 * P(t + epsilon) = P(t) + F(t) * epsilon
 * ================================================================================*/
void EvolveMomentaSigma(struct Momentasigma * momenta, const struct Forcesigma * force, double epsilon)
{
    // for all momenta in VOLUME
    for(int x=0; x<VOLUME; x++)
    momenta->momenta[x] += (force->force[x] * epsilon);
}

/**
 * 2^nd order LeapFrog (LF2) for the sigma field.
 * The function initially draw the momenta (randomly from a gaussian distribution), then evolve globally the sigma field.
 * The fundamental infinitesimal step in the LF2 algorithm is given by the following three steps:
 *                      - P(t)                 ->               P( t  +  \epsilon / 2  )
 *                      - Q(t)                 ->               Q( t  +  \epsilon  )
 *                      - P(t+\epsilon/2)      ->               P( t  +   \epsilon  )
 * These steps are repeated an integer amount of times equal to "steps", which is passed to the function.
 */
bool LeapFrog2Sigma(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file)
{
    #ifdef HYBRID_MC
    // file where we will write the Leapfrog acceptance parameters
    FILE * opf = fopen(acceptance_file, "a");
    // first draw initial momenta
    DrawSigmaMomentaFromGaussianDistribution(momenta);

    #ifdef DEBUG_LF2_SIGMA
    printf("S : Total sigma squared is %.16f\n", ComputeTotalSigmaSquared(lattice));
    #endif
    
    // draw also a random fermionic field to be associated with the chi field
    // for all flavors in the euclidean action
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    {
        DrawFullFermionFromGaussianDistribution(&(kit->chi.flav[f]));
        // compute phi at the beginning of the leapfrog
        DiracFermion(lattice, &(kit->chi.flav[f]), &(kit->phi.flav[f]), info, fermioncoor);
    }
    #endif


    // register a copy of the sigma variables
    double sigma_copy[VOLUME];
    CopyFromSigma(lattice->sigma, sigma_copy);


    // we have all needed to compute the initial fictitious energy
    // invert all fermions
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    if(!StartX0ImposedConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    #endif
    
    // compute the starting energy of the molecular dynamics
    double starting_energy = ComputeFullEnergyForSigma(lattice, info, fermioncoor, momenta, kit);
    #ifdef DEBUG_STARTING_ENERGY_HMC
    printf("Delta energy is %.16f\n", ComputeFullEnergyForSigmaWithChi(lattice, info, fermioncoor, momenta, kit) - starting_energy);
    #endif

    //           P(t)        ->      P( t  +  \epsilon / 2  )
    // first compute the force, then evolve in time
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, (epsilon / 2.0));

    // for the number of steps passed to the function
    for(int s=0; s<step-1; s++)
    {        
        //          Q(t)         ->      Q( t  +  \epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon);

        //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
        // first, compute again the force, then evolve in time
        if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
        {
            CopyToSigma(sigma_copy, lattice->sigma);
            ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
            return false;
        }
        EvolveMomentaSigma(momenta, force, epsilon);
    }

    //          Q(t)         ->      Q( t  +  \epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon);

    //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
    // first, compute again the force, then evolve in time
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, HIGH_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, epsilon / 2.0);

    #ifdef DEBUG_LF2_SIGMA
    for(int x=0; x<VOLUME; x++) momenta->momenta[x] = -momenta->momenta[x];
    // evolve in time with all momenta reversed
    EvolveSigmaInTimeLF2(lattice, info, force, momenta, fermioncoor, kit, epsilon, step);
    printf("F : Total sigma squared is %.16f\n\n\n", ComputeTotalSigmaSquared(lattice));
    //--------------------------------------------------------------------------------------------------------------------------
    #endif    

    // compute the final fictitious energy
    double final_energy = ComputeFullEnergyForSigma(lattice, info, fermioncoor, momenta, kit);
    
    // accept - reject test
    if(Casuale() <= exp( - (final_energy - starting_energy) ) )
    {
        // the lattice has already been updated
        info->hybrid_accepted_sigma++;
        info->hybrid_done_sigma++;
        info->hybrid_done_acceptance_sigma = (double) (info->hybrid_accepted_sigma) / (double) (info->hybrid_done_sigma);
        // local info
        info->local_accepted_sigma++;
        info->local_done_sigma++;
        info->local_acceptance_sigma = (double) (info->local_accepted_sigma) / (double) (info->local_done_sigma);
        // lattice done ++
        lattice->done_sigma++;
        // update story sigma
        // if the story is in a non-negative trend, add a positive event
        // else set the story to -1
        if(info->story_sigma >= 0) info->story_sigma++;
        else info->story_sigma = -1;
    }
    else
    {
        // if rejected copy back the all sigma fields
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);

        // update info
        info->hybrid_done_sigma++;
        info->hybrid_done_acceptance_sigma = (double) (info->hybrid_accepted_sigma) / (double) (info->hybrid_done_sigma);
        // local info
        info->local_done_sigma++;
        info->local_acceptance_sigma = (double) (info->local_accepted_sigma) / (double) (info->local_done_sigma);
        // lattice done ++
        lattice->done_sigma++;
        
        // update the story sigma
        // if the story is in a positive trend, set the trend to -1
        // else add a negative event
        if(info->story_sigma >= 0) info->story_sigma = -1;
        else info->story_sigma--;
    }

    fprintf(opf, "%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\n", info->beta, info->beta_scalar, info->mass, info->epsilon_link_HMC, info->steps_link_HMC, (info->epsilon_link_HMC * info->steps_link_HMC), info->delta_H_link, info->hybrid_done_acceptance_link, info->epsilon_sigma_HMC, info->steps_sigma_HMC, (info->epsilon_sigma_HMC * info->steps_sigma_HMC), final_energy-starting_energy, info->hybrid_done_acceptance_sigma);

    // close the acceptance file
    fclose(opf);

    return true;
    #endif
}

/**
 * 2^nd order LeapFrogMinimumNorm (LFMN2) for the sigma field follows.
 * The function initially draw the momenta (randomly from a gaussian distribution), then evolve globally the sigma field.
 * The fundamental infinitesimal step in the LFMN2 algorithm is given by the following three steps:
 *                      - Q(t)                 ->               Q( t  + lambda * epsilon  )
 *                      - P(t)                 ->               P( t  +  \epsilon / 2 )
 *                      - Q( t + lambda * epsilon)      ->             Q( t  +  (1-lambda) \epsilon  )
 *                      - P(t + epsilon/2)              ->             P( t   +  epsilon)
 *                      - Q( t + (1-lambda) * epsilon)   ->          Q(t + epsilon)
 * These steps are repeated an integer amount of times equal to "steps", which is passed to the function.
 */
bool LeapFrog2MinimumNormSigma(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file)
{
    #ifdef HYBRID_MC
    // file where we will write the Leapfrog acceptance parameters
    FILE * opf = fopen(acceptance_file, "a");
    // first draw initial momenta
    DrawSigmaMomentaFromGaussianDistribution(momenta);

    #ifdef DEBUG_LF2_SIGMA
    printf("S : Total sigma squared is %.16f\n", ComputeTotalSigmaSquared(lattice));
    #endif
    
    // draw also a random fermionic field to be associated with the chi field
    // for all flavors in the euclidean action
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    {
        DrawFullFermionFromGaussianDistribution(&(kit->chi.flav[f]));
        // compute phi at the beginning of the leapfrog
        DiracFermion(lattice, &(kit->chi.flav[f]), &(kit->phi.flav[f]), info, fermioncoor);
    }
    #endif


    // register a copy of the sigma variables
    double sigma_copy[VOLUME];
    CopyFromSigma(lattice->sigma, sigma_copy);


    // we have all we need to compute the initial fictitious energy
    // invert all fermions
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    if(!StartX0ImposedConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    #endif
    
    // compute the starting energy of the molecular dynamics
    double starting_energy = ComputeFullEnergyForSigma(lattice, info, fermioncoor, momenta, kit);
    #ifdef DEBUG_STARTING_ENERGY_HMC
    printf("Delta energy is %.16f\n", ComputeFullEnergyForSigmaWithChi(lattice, info, fermioncoor, momenta, kit) - starting_energy);
    #endif
    // -----------------  STEP ZERO --------------------------------------------------------
    //           Q(t)        ->      Q( t  +  lambda * epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (MAGIC_LAMBDA * epsilon));

    //          P(t)         ->       P(t + epsilon / 2)
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, (epsilon / 2.0));

    //           Q(t + lambda * epsilon)        ->      Q( t  + (1 - lambda) * epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, ((1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon));   

    //          P(t)         ->       P(t + epsilon / 2)
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, (epsilon / 2.0));

    // ---------------------      CYCLE     -------------------------------------
    for(int s=1; s<step-1; s++)
    {        
        //          Q(t)         ->      Q( t  +  2 * lambda * epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (2.0 * MAGIC_LAMBDA * epsilon));

        //      P(t)   ->     P( t  +   epsilon / 2  )
        // first, compute again the force, then evolve in time
        if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
        {
            CopyToSigma(sigma_copy, lattice->sigma);
            ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
            return false;
        }
        EvolveMomentaSigma(momenta, force, (epsilon/2.0));

        //          Q(t + 2 * lambda * epsilon)         ->      Q( t  +  epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, ((1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon));   

        //      P(t + espilon/2)   ->     P( t  +   epsilon )
        // first, compute again the force, then evolve in time
        if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
        {
            CopyToSigma(sigma_copy, lattice->sigma);
            ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
            return false;
        }
        EvolveMomentaSigma(momenta, force, (epsilon/2.0));             
    }

    //---------------------- Last step ------------------------------------------------
    //          Q(t - lambda * epsilon)         ->      Q( t  +  lambda * epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (2.0 * MAGIC_LAMBDA * epsilon));
    
    //         P(t)         ->       P(t + epsilon/2)
    //          P(t)         ->       P(t + epsilon / 2)
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, (epsilon / 2.0));

    //          Q(t + lambda * epsilon)         ->      Q( t  +  ( 1 - lambda) \epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, ((1 - (2.0 * MAGIC_LAMBDA)) * epsilon));

    //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
    // first, compute again the force, then evolve in time
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, (epsilon / 2.0));

    //          Q(t + (1 - lambda) epsilon)         ->      Q( t  +  \epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (MAGIC_LAMBDA * epsilon));

    // we have all needed to compute the initial fictitious energy
    // invert all fermions
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    if(!StartX0ImposedConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    #endif

    #ifdef DEBUG_LF2_SIGMA
    for(int x=0; x<VOLUME; x++) momenta->momenta[x] = -momenta->momenta[x];
    // Evolve with time epsilon with all momenta reversed
    EvolveSigmaInTimeMinimumNormLF2(lattice, info, force, momenta, fermioncoor, kit, epsilon, step);
    printf("F : Total sigma squared is %.16f\n\n\n", ComputeTotalSigmaSquared(lattice));
    //--------------------------------------------------------------------------------------------------------------------------
    #endif    

    // compute the final fictitious energy
    double final_energy = ComputeFullEnergyForSigma(lattice, info, fermioncoor, momenta, kit);
    
    // accept - reject test
    if(Casuale() <= exp( - (final_energy - starting_energy) ) )
    {
        // the lattice has already been updated
        info->hybrid_accepted_sigma++;
        info->hybrid_done_sigma++;
        info->hybrid_done_acceptance_sigma = (double) (info->hybrid_accepted_sigma) / (double) (info->hybrid_done_sigma);
        // local info
        info->local_accepted_sigma++;
        info->local_done_sigma++;
        info->local_acceptance_sigma = (double) (info->local_accepted_sigma) / (double) (info->local_done_sigma);
        // lattice done ++
        lattice->done_sigma++;

        // update story sigma
        // if the story is in a non-negative trend, add a positive event
        // else set the story to -1
        if(info->story_sigma >= 0) info->story_sigma++;
        else info->story_sigma = -1;
    }
    else
    {
        // if rejected copy back the all sigma fields
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);

        // update info
        info->hybrid_done_sigma++;
        info->hybrid_done_acceptance_sigma = (double) (info->hybrid_accepted_sigma) / (double) (info->hybrid_done_sigma);
        // local info
        info->local_done_sigma++;
        info->local_acceptance_sigma = (double) (info->local_accepted_sigma) / (double) (info->local_done_sigma);
        // lattice done ++
        lattice->done_sigma++;

        // update the story sigma
        // if the story is in a positive trend, set the trend to -1
        // else add a negative event
        if(info->story_sigma >= 0) info->story_sigma = -1;
        else info->story_sigma--;
    }

    fprintf(opf, "%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\n", info->beta, info->beta_scalar, info->mass, info->epsilon_link_HMC, info->steps_link_HMC, (info->epsilon_link_HMC * info->steps_link_HMC), info->delta_H_link, info->hybrid_done_acceptance_link, info->epsilon_sigma_HMC, info->steps_sigma_HMC, (info->epsilon_sigma_HMC * info->steps_sigma_HMC), final_energy-starting_energy, info->hybrid_done_acceptance_sigma);

    // close the acceptance file
    fclose(opf);

    return true;
    #endif
}

/**
 * 2^nd order LeapFrogMinimumNormVelocityFirst (LFMN2VF) for the sigma field follows.
 * The function initially draw the momenta (randomly from a gaussian distribution), then evolve globally the sigma field.
 * The fundamental infinitesimal step in the LFMN2 algorithm is given by the following three steps:
 *                      - P(t)                 ->               P( t  + lambda * epsilon  )
 *                      - Q(t)                 ->               Q( t  +  \epsilon / 2 )
 *                      - P( t + lambda * epsilon)      ->             P( t  +  (1-lambda) \epsilon  )
 *                      - Q(t + epsilon/2)              ->             Q( t   +  epsilon)
 *                      - P( t + (1-lambda) * epsilon)   ->          P(t + epsilon)
 * These steps are repeated an integer amount of times equal to "steps", which is passed to the function.
 */
bool LeapFrog2MinimumNormSigmaVelocityFirst(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file)
{
    #ifdef HYBRID_MC
    // file where we will write the Leapfrog acceptance parameters
    FILE * opf = fopen(acceptance_file, "a");
    // first draw initial momenta
    DrawSigmaMomentaFromGaussianDistribution(momenta);
    
    // draw also a random fermionic field to be associated with the chi field
    // for all flavors in the euclidean action
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    {
        DrawFullFermionFromGaussianDistribution(&(kit->chi.flav[f]));
        // compute phi at the beginning of the leapfrog
        DiracFermion(lattice, &(kit->chi.flav[f]), &(kit->phi.flav[f]), info, fermioncoor);
    }
    #endif


    // register a copy of the sigma variables
    double sigma_copy[VOLUME];
    CopyFromSigma(lattice->sigma, sigma_copy);


    // we have all we need to compute the initial fictitious energy
    // invert all fermions
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    if(!StartX0ImposedConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    #endif
    
    // compute the starting energy of the molecular dynamics
    double starting_energy = ComputeFullEnergyForSigma(lattice, info, fermioncoor, momenta, kit);
    #ifdef DEBUG_STARTING_ENERGY_HMC
    printf("Delta energy is %.16f\n", ComputeFullEnergyForSigmaWithChi(lattice, info, fermioncoor, momenta, kit) - starting_energy);
    #endif
    // -----------------  STEP ZERO --------------------------------------------------------

    //          P(t)         ->       P(t + lambda * epsilon)
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, (MAGIC_LAMBDA * epsilon));

    //           Q(t)        ->      Q( t  +  epsilon / 2 )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (epsilon/2.0));

    //          P(t + lambda * epsilon)         ->       P(t + (1 - lambda) * epsilon)
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force,  ((1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon));

    //           Q(t + epsilon / 2)        ->      Q( t  + epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (epsilon / 2.0));   

    // ---------------------      CYCLE     -------------------------------------
    for(int s=1; s<step-1; s++)
    {        

        //      P(t)   ->     P( t  + 2 * lambda * epsilon)
        // first, compute again the force, then evolve in time
        if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
        {
            CopyToSigma(sigma_copy, lattice->sigma);
            ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
            return false;
        }
        EvolveMomentaSigma(momenta, force, (2.0 * MAGIC_LAMBDA * epsilon));
        
        //          Q(t)         ->      Q( t  +  epsilon / 2)
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (epsilon / 2.0));

        //      P(t + 2 * lambda * espilon)   ->     P( t  +   epsilon )
        // first, compute again the force, then evolve in time
        if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
        {
            CopyToSigma(sigma_copy, lattice->sigma);
            ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
            return false;
        }
        EvolveMomentaSigma(momenta, force, ((1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon));    

        //          Q(t + epsilon/2)         ->      Q( t  +  epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (epsilon / 2.0));            
    }

    //---------------------- Last step ------------------------------------------------
    //          P(t - lambda)         ->       P(t + lambda * epsilon)
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, (2.0 * MAGIC_LAMBDA * epsilon));

    //          Q(t)         ->      Q( t  +  epsilon / 2 )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (epsilon / 2.0));

    //      P(t+ lambda * epsilon)   ->     P( t  + (1 - lambda) * epsilon  )
    // first, compute again the force, then evolve in time
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, ((1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon));

    //          Q(t + epsilon / 2)         ->      Q( t  + epsilon )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (epsilon / 2.0));


    //      P(t - lambda * epsilon)   ->     P( t  + epsilon  )
    // first, compute again the force, then evolve in time
    if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, HIGH_PRECISION_INVERSION))
    {
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        return false;
    }
    EvolveMomentaSigma(momenta, force, (MAGIC_LAMBDA * epsilon));
    

    #ifdef DEBUG_LF2_SIGMA
    // Evolve with time all momemnta reversed
    for(int x=0; x<VOLUME; x++) momenta->momenta[x] = - momenta->momenta[x];
    EvolveSigmaInTimeMinimumNormSigmaVelocityFirst(lattice, info, force, momenta, fermioncoor, kit, epsilon, step);
    //--------------------------------------------------------------------------------------------------------------------------
    #endif    

    // compute the final fictitious energy
    double final_energy = ComputeFullEnergyForSigma(lattice, info, fermioncoor, momenta, kit);
    #ifdef DEBUG_LF2_SIGMA
    if(fabs(final_energy - starting_energy) > 1e-14) printf("F: Delta energy is %.16f\n\n\n", fabs(final_energy - starting_energy));
    #endif   

    // accept - reject test
    if(Casuale() <= exp( - (final_energy - starting_energy) ) )
    {
        // the lattice has already been updated
        info->hybrid_accepted_sigma++;
        info->hybrid_done_sigma++;
        info->hybrid_done_acceptance_sigma = (double) (info->hybrid_accepted_sigma) / (double) (info->hybrid_done_sigma);
        // local info
        info->local_accepted_sigma++;
        info->local_done_sigma++;
        info->local_acceptance_sigma = (double) (info->local_accepted_sigma) / (double) (info->local_done_sigma);
        // lattice done ++
        lattice->done_sigma++;

        // update story sigma
        // if the story is in a non-negative trend, add a positive event
        // else set the story to -1
        if(info->story_sigma >= 0) info->story_sigma++;
        else info->story_sigma = -1;
    }
    else
    {
        // if rejected copy back the all sigma fields
        CopyToSigma(sigma_copy, lattice->sigma);
        ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);

        // update info
        info->hybrid_done_sigma++;
        info->hybrid_done_acceptance_sigma = (double) (info->hybrid_accepted_sigma) / (double) (info->hybrid_done_sigma);
        // local info
        info->local_done_sigma++;
        info->local_acceptance_sigma = (double) (info->local_accepted_sigma) / (double) (info->local_done_sigma);
        // lattice done ++
        lattice->done_sigma++;

        // update the story sigma
        // if the story is in a positive trend, set the trend to -1
        // else add a negative event
        if(info->story_sigma >= 0) info->story_sigma = -1;
        else info->story_sigma--;
    }

    fprintf(opf, "%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\n", info->beta, info->beta_scalar, info->mass, info->epsilon_link_HMC, info->steps_link_HMC, (info->epsilon_link_HMC * info->steps_link_HMC), info->delta_H_link, info->hybrid_done_acceptance_link, info->epsilon_sigma_HMC, info->steps_sigma_HMC, (info->epsilon_sigma_HMC * info->steps_sigma_HMC), final_energy-starting_energy, info->hybrid_done_acceptance_sigma);

    // close the acceptance file
    fclose(opf);

    return true;
    #endif

}

/**
 * 4^nd order LeapFrog4MinimumNormFirstPosition (LF4MN4FP) for the sigma field in terms of the constant rho, theta and lambda follows.
 * The function initially draw the momenta and fermion fields (randomly from a gaussian distribution), then evolve globally the sigma field.
 * The fundamental infinitesimal step in the LF4MN4FP algorithm is given by the following three steps:
 *                      - Q(t)                 ->               Q( t  + rho * epsilon  )
 *                      - P(t)                 ->               P( t  +  lambda * epsilon )
 *                      - Q( t + rho * epsilon)      ->             Q( t  +  (rho + theta) * epsilon  )
 *                      - P(t + lambda * epsilon)              ->             P( t   + 0.5 * epsilon)
 *                      - Q( t + (rho + theta) * epsilon)   ->          Q(t + (1 - (rho + theta)) * epsilon)
 *                      - P(t + 0.5 * epsilon)                 ->               P( t  +  (1 - lambda) * epsilon )
 *                      - Q( t + (1 - (rho + theta)) * epsilon)      ->             Q( t  +  (1 - rho) * epsilon  )
 *                      - P(t + (1 - lambda) * epsilon)              ->             P( t + epsilon)
 *                      - Q( t + (1 - rho) * epsilon)      ->             Q( t  + epsilon  )
 * These steps are repeated an integer amount of times equal to "steps", which is passed to the function.
 */
bool LeapFrog4MinimumNorm4SigmaFirstPosition(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file)
{
    ;
}

/**
 * The function is useful to check the reveresibility of the HMC.
 * If the momenta sigma have been flipped, the evolution in time of sigma should take us back to the starting sigma.
 */
void EvolveSigmaInTimeLF2(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step)
{
    //           P(t)        ->      P( t  +  \epsilon / 2  )
    // first compute the force, then evolve in time
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, HIGH_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, epsilon / 2.0);

    // for the number of steps passed to the function
    for(int s=0; s<step-1; s++)
    {        
        //          Q(t)         ->      Q( t  +  \epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon);

        //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
        // first, compute again the force, then evolve in time
        ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
        EvolveMomentaSigma(momenta, force, epsilon);
    }

    //          Q(t)         ->      Q( t  +  \epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon);

    //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
    // first, compute again the force, then evolve in time
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, epsilon / 2.0);

}

void EvolveSigmaInTimeMinimumNormLF2(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step)
{
    // -----------------  STEP ZERO --------------------------------------------------------
    //           Q(t)        ->      Q( t  +  lambda * epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, MAGIC_LAMBDA * epsilon);

    //          P(t)         ->       P(t + epsilon / 2)
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, HIGH_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, epsilon/2.0);

    //           Q(t + lambda * epsilon)        ->      Q( t  + (1 - lambda) * epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon);   

    //          P(t + epsilon/2)         ->       P(t + epsilon)
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, epsilon/2.0);

    // ---------------------      CYCLE     -------------------------------------
    for(int s=1; s<step-1; s++)
    {        
        //          Q(t)         ->      Q( t  +  2 * lambda * epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, 2.0 * MAGIC_LAMBDA * epsilon);

        //      P(t)   ->     P( t  +   epsilon / 2  )
        // first, compute again the force, then evolve in time
        ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
        EvolveMomentaSigma(momenta, force, epsilon/2.0);

        //          Q(t + 2 * lambda * epsilon)         ->      Q( t  +  epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon);   

        //      P(t + espilon/2)   ->     P( t  +   epsilon )
        // first, compute again the force, then evolve in time
        ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
        EvolveMomentaSigma(momenta, force, epsilon/2.0);  
           
    }

    //---------------------- Last step ------------------------------------------------
    //          Q(t - lambda * epsilon)         ->      Q( t  +  lambda * epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, 2.0 * MAGIC_LAMBDA * epsilon);
    
    //          P(t)         ->       P(t + epsilon / 2)
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, epsilon / 2.0);

    //          Q(t + lambda * epsilon)         ->      Q( t  +  ( 1 - lambda) \epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, (1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon);

    //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
    // first, compute again the force, then evolve in time
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, epsilon / 2.0);

    //          Q(t + (1 - lambda) epsilon)         ->      Q( t  +  \epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, MAGIC_LAMBDA * epsilon);
}

/**
 * Evolve in time the sigma field according to the integrator LF2MN2VF.
 * No energy is computed as the function should only check the reversibility of the HMC algorithm.
 */
void EvolveSigmaInTimeMinimumNormSigmaVelocityFirst(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step)
{
    // we have all we need to compute the initial fictitious energy
    // invert all fermions
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    StartX0ImposedConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f);
    #endif
    
    // -----------------  STEP ZERO --------------------------------------------------------

    //          P(t)         ->       P(t + lambda * epsilon)
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, MAGIC_LAMBDA * epsilon);

    //           Q(t)        ->      Q( t  +  epsilon / 2 )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon/2.0);

    //          P(t + lambda * epsilon)         ->       P(t + (1 - lambda) * epsilon)
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force,  (1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon);

    //           Q(t + epsilon / 2)        ->      Q( t  + epsilon  )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon / 2.0);   

    // ---------------------      CYCLE     -------------------------------------
    for(int s=1; s<step-1; s++)
    {        

        //      P(t)   ->     P( t  + 2 * lambda * epsilon)
        // first, compute again the force, then evolve in time
        ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
        EvolveMomentaSigma(momenta, force, 2.0 * MAGIC_LAMBDA * epsilon);
        
        //          Q(t)         ->      Q( t  +  epsilon / 2)
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon/2.0);

        //      P(t + 2 * lambda * espilon)   ->     P( t  +   epsilon )
        // first, compute again the force, then evolve in time
        ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
        EvolveMomentaSigma(momenta, force, (1.0 - (2.0 * MAGIC_LAMBDA)) * epsilon);    

        //          Q(t + epsilon/2)         ->      Q( t  +  epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon/2.0);            
    }

    //---------------------- Last step ------------------------------------------------
    //          P(t - lambda)         ->       P(t + lambda * epsilon)
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, 2.0 * MAGIC_LAMBDA * epsilon);

    //          Q(t)         ->      Q( t  +  epsilon / 2 )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon / 2.0);

    //      P(t+ lambda * epsilon)   ->     P( t  + (1 - lambda) * epsilon  )
    // first, compute again the force, then evolve in time
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, (1 - (2.0 * MAGIC_LAMBDA)) * epsilon);

    //          Q(t + epsilon / 2)         ->      Q( t  + epsilon )
    EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon / 2.0);


    //      P(t - lambda * epsilon)   ->     P( t  + epsilon  )
    // first, compute again the force, then evolve in time
    ComputeSigmaForce(lattice, info, fermioncoor, force, kit, HIGH_PRECISION_INVERSION);
    EvolveMomentaSigma(momenta, force, MAGIC_LAMBDA * epsilon);

    // we have all needed to compute the initial fictitious energy
    // invert all fermions
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f);
    #endif
}

void FreeLeapFrog2SigmaWithoutAcceptance(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step)
{
    for(int num=0; num<10; num++)
    {
        #ifdef HYBRID_MC
        // first draw initial momenta
        DrawSigmaMomentaFromGaussianDistribution(momenta);
        
        // draw also a random fermionic field to be associated with the chi field
        // for all flavors in the euclidean action
        #ifdef PSEUDO_FERMION
        for(int f=0; f<HALFNFLAV; f++)
        {
            DrawFullFermionFromGaussianDistribution(&(kit->chi.flav[f]));
            // compute phi at the beginning of the leapfrog
            DiracFermion(lattice, &(kit->chi.flav[f]), &(kit->phi.flav[f]), info, fermioncoor);
        }
        #endif


        // register a copy of the sigma variables
        double sigma_copy[VOLUME];
        CopyFromSigma(lattice->sigma, sigma_copy);


        // we have all needed to compute the initial fictitious energy
        // invert all fermions
        #ifdef PSEUDO_FERMION
        for(int f=0; f<HALFNFLAV; f++)
        if(!StartX0ImposedConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f))
        {
            CopyToSigma(sigma_copy, lattice->sigma);
            ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        }
        #endif
        
        //           P(t)        ->      P( t  +  \epsilon / 2  )
        // first compute the force, then evolve in time
        if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
        {
            CopyToSigma(sigma_copy, lattice->sigma);
            ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        }
        EvolveMomentaSigma(momenta, force, epsilon / 2.0);

        // for the number of steps passed to the function
        for(int s=0; s<step-1; s++)
        {        
            //          Q(t)         ->      Q( t  +  \epsilon  )
            EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon);

            //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
            // first, compute again the force, then evolve in time
            if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, MEDIUM_PRECISION_INVERSION))
            {
                CopyToSigma(sigma_copy, lattice->sigma);
                ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
            }
            EvolveMomentaSigma(momenta, force, epsilon);
        }

        //          Q(t)         ->      Q( t  +  \epsilon  )
        EvolveCoordinatesSigma(lattice, momenta, fermioncoor, epsilon);

        //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
        // first, compute again the force, then evolve in time
        if(!ComputeSigmaForce(lattice, info, fermioncoor, force, kit, HIGH_PRECISION_INVERSION))
        {
            CopyToSigma(sigma_copy, lattice->sigma);
            ComputeFullSigmaMass(lattice, fermioncoor, lattice->sigmamass);
        }
        EvolveMomentaSigma(momenta, force, epsilon / 2.0);

        #endif
    }
}

//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////

/**
 * Evolve the pi field in time of an infinitesimal step epsilon.
 */
void EvolveCoordinatesPi(struct Lattice * lat, const struct Momentapi const * momenta, const struct Fermioncoor const * fcor, double epsilon)
{
    // for all VOLUME fermionic sites
    // evolve the pi field of an infinitesimal step epsilon * momenta[\tilde{x}]
    for(int x=0; x<VOLUME; x++)
    lat->pi[x] += (epsilon * momenta->momenta[x]);

    // after the pi field has been evolved, compute the updated effective pi mass
    ComputeFullPiMass(lat, fcor, lat->pimass);
}

/**
 * Evolve all the momenta of the pi fields of an infinitesimal time step epsilon.
 * ==================================================================================
 * P(t + epsilon) = P(t) + F(t) * epsilon
 * ================================================================================*/
void EvolveMomentaPi(struct Momentapi * momenta, const struct Forcepi * force, double epsilon)
{
    // for all momenta in VOLUME
    for(int x=0; x<VOLUME; x++)
    momenta->momenta[x] += (force->force[x] * epsilon);
}

/**
 * 2^nd order LeapFrog (LF2) for the pi field.
 * The function initially draw the momenta (randomly from a gaussian distribution), then evolve globally the pi field.
 * The fundamental infinitesimal step in the LF2 algorithm is given by the following three steps:
 *                      - P(t)                 ->               P( t  +  \epsilon / 2  )
 *                      - Q(t)                 ->               Q( t  +  \epsilon  )
 *                      - P(t+\epsilon/2)      ->               P( t  +   \epsilon  )
 * These steps are repeated an integer amount of times equal to "steps", which is passed to the function.
 */
bool LeapFrog2Pi(struct Lattice * lattice, struct Info * info, struct Forcepi * force, struct Momentapi * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file)
{
    #ifdef HYBRID_MC
    // file where we will write the Leapfrog acceptance parameters
    FILE * opf = fopen(acceptance_file, "a");
    // first draw initial momenta
    DrawPiMomentaFromGaussianDistribution(momenta);

    #ifdef DEBUG_LF2_PI
    printf("S : Energy is %.16f\n", ComputeTotalPiSquared(lattice));
    #endif
    
    // draw also a random fermionic field to be associated with the chi field
    // for all flavors in the euclidean action
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    {
        DrawFullFermionFromGaussianDistribution(&(kit->chi.flav[f]));
        // compute phi at the beginning of the leapfrog
        DiracFermion(lattice, &(kit->chi.flav[f]), &(kit->phi.flav[f]), info, fermioncoor);
    }
    #endif


    // register a copy of the sigma variables
    double pi_copy[VOLUME];
    CopyFromPi(lattice->pi, pi_copy);


    // we have all needed to compute the initial fictitious energy
    // invert all fermions
    #ifdef PSEUDO_FERMION
    for(int f=0; f<HALFNFLAV; f++)
    if(!ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lattice, &(kit->phi.flav[f]), &(kit->inverted_phi.flav[f]), info, fermioncoor, HIGH_PRECISION_INVERSION, kit, f))
    {
        CopyToPi(pi_copy, lattice->pi);
        ComputeFullPiMass(lattice, fermioncoor, lattice->pimass);
        return false;
    }
    #endif
    
    // compute the starting energy of the molecular dynamics
    double starting_energy = ComputeFullEnergyForPi(lattice, info, fermioncoor, momenta, kit);
    
    //           P(t)        ->      P( t  +  \epsilon / 2  )
    // first compute the force, then evolve in time
    if(!ComputePiForce(lattice, info, fermioncoor, force, kit))
    {
        CopyToPi(pi_copy, lattice->pi);
        ComputeFullPiMass(lattice, fermioncoor, lattice->pimass);
        return false;
    }
    EvolveMomentaPi(momenta, force, epsilon / 2.0);

    // for the number of steps passed to the function
    for(int s=0; s<step-1; s++)
    {        
        //          Q(t)         ->      Q( t  +  \epsilon  )
        EvolveCoordinatesPi(lattice, momenta, fermioncoor, epsilon);

        //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
        // first, compute again the force, then evolve in time
        if(!ComputePiForce(lattice, info, fermioncoor, force, kit))
        {
            CopyToPi(pi_copy, lattice->pi);
            ComputeFullPiMass(lattice, fermioncoor, lattice->pimass);
            return false;
        }
        EvolveMomentaPi(momenta, force, epsilon);
    }

    //          Q(t)         ->      Q( t  +  \epsilon  )
    EvolveCoordinatesPi(lattice, momenta, fermioncoor, epsilon);

    //      P(t+\epsilon/2)   ->     P( t  +   \epsilon  )
    // first, compute again the force, then evolve in time
    if(!ComputePiForce(lattice, info, fermioncoor, force, kit))
    {
        CopyToPi(pi_copy, lattice->pi);
        ComputeFullPiMass(lattice, fermioncoor, lattice->pimass);
        return false;
    }
    EvolveMomentaPi(momenta, force, epsilon / 2.0);

    #ifdef DEBUG_LF2_PI
    // invert all momenta
    for(int x=0; x<VOLUME; x++) momenta->momenta[x] = -momenta->momenta[x];

    // Evolve in time with momenta inverted
    EvolvePiInTime(lattice, info, force, momenta, fermioncoor, kit, epsilon, step);
    printf("F : Energy is %.16f\n\n\n", ComputeTotalPiSquared(lattice));
    //--------------------------------------------------------------------------------------------------------------------------
    #endif    

    // compute the final fictitious energy
    double final_energy = ComputeFullEnergyForPi(lattice, info, fermioncoor, momenta, kit);
    printf("Delta energy is %.16f\n", final_energy - starting_energy);
    
    // accept - reject test
    if(Casuale() <= exp( - (final_energy - starting_energy) ) )
    {
        // the lattice has already been updated
        info->hybrid_accepted_pi++;
        info->hybrid_done_pi++;
        info->hybrid_done_acceptance_pi = (double) (info->hybrid_accepted_pi) / (double) (info->hybrid_done_pi);
        
        // update story pi
        // if the story is in a non-negative trend, add a positive event
        // else set the story to -1
        if(info->story_pi >= 0) info->story_pi++;
        else info->story_pi = -1;
    }
    else
    {
        // if rejected copy back the all sigma fields
        CopyToPi(pi_copy, lattice->pi);
        ComputeFullPiMass(lattice, fermioncoor, lattice->pimass);

        // update info
        info->hybrid_done_pi++;
        info->hybrid_done_acceptance_pi = (double) (info->hybrid_accepted_pi) / (double) (info->hybrid_done_pi);

        // update the story pi
        // if the story is in a positive trend, set the trend to -1
        // else add a negative event
        if(info->story_pi >= 0) info->story_pi = -1;
        else info->story_pi--;
    }
    printf("%d/%d\n", info->hybrid_accepted_pi, info->hybrid_done_pi);

    fprintf(opf, "%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\t%.16f\t%d\t%.16f\t%.16f\t%.16f\n", info->beta, info->beta_scalar, info->mass, info->epsilon_link_HMC, info->steps_link_HMC, (info->epsilon_link_HMC * info->steps_link_HMC), info->delta_H_link, info->hybrid_done_acceptance_link, info->epsilon_sigma_HMC, info->steps_sigma_HMC, (info->epsilon_sigma_HMC * info->steps_sigma_HMC), final_energy-starting_energy, info->hybrid_done_acceptance_sigma);

    // close the acceptance file
    fclose(opf);

    return true;
    #endif
}
