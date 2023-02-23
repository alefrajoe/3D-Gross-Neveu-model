#include "dirac.h"

///////////////////////////////////////////////////////////////////
//                  Fast Dirac functions
///////////////////////////////////////////////////////////////////

/**
 * Print a fast half dmatrix to stdout.
 * A fast half matrix is composed by the first HALFVOLUME columns of a complete Dirac matrix
 * The sites are ordered according to fermioncoordinates site.
 */
void PrintDirac(const struct Lattice const * lat, const struct Fermioncoor const * fermioncoor, const struct Info const * info)
{
    printf("=========== FastHalfDmatrix ===========\n");
    // print to stdout all the columns contained in the matrix
    for(int x=0; x<HALFVOLUME; x++)
    {
        printf("========== fermionic site %d ============\n", x);
        // mass term
        printf("FastDirac [%d][%d] is %.16f\n", x, x, info->mass);

        int site = fermioncoor->site[x]; // stereo is now an even site

        // U_mu(x)
        for(int d=0; d<(DIM); d++)
        {
            complex temp = (0.5) * lat->eta_pp[COOR(lat->mm[COOR(site, d)], d)] * lat->dlink[COOR(lat->mm[COOR(site, d)], d)].link;
            printf("FastDirac [%d][%d] is %.16f + i %.16f\n", fermioncoor->fermionmm[COOR(x, d)], x, creal(temp), cimag(temp));
        }

        // U^*_mu(x)
        for(int d=0; d<DIM; d++)
        {
            complex temp = - (0.5) * lat->eta_pp[COOR(lat->pp[COOR(site, d)], d)] * conj(lat->dlink[COOR(site, d)].link);
            printf("FastDirac [%d][%d] is %.16f + i %.16f\n", fermioncoor->fermionpp[COOR(x, d)], x, creal(temp), cimag(temp));
        }
        printf("============================\n");
    }
}

/**
 * Compute the product of a Dirac matrix and a full fermion vector.
 * The resulting fermion is saved into target.
 * The function exploits the knowledge of only M (the mass matrix) and ( D_oe )
 */
void DiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor)
{
    // mass term for both the upper part of the fermion vector (U) and the bottom part (D)
    // the mass term initialize the target fermion
    // if SIGMA_INTERACTION is defined, consider also the presence of an effective mass term due to sigma
    for(int x=0; x<VOLUME; x++)
    {
        #ifndef SIGMA_INTERACTION
        #ifndef PI_INTERACTION
        target->fermion[x] = 0.0;
        #endif
        #endif

        #ifndef PI_INTERACTION
        #ifdef SIGMA_INTERACTION
        target->fermion[x] = (lat->sigmamass[x]) * fermion->fermion[x];
        #endif
        #endif

        #ifdef SIGMA_INTERACTION
        #ifdef PI_INTERACTION
        target->fermion[x] = (lat->sigmamass[x] + (I_UNIT * lat->epsilon_pi[x] * lat->pimass[x])) * fermion->fermion[x];
        #endif
        #endif
    }
    

    // Out-of-diagonal D terms
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        // + (0.5) \eta_mu(x) * U_mu(x) \delta_{x, x+mu}
        target->fermion[x] += (0.5) * fermioncoor->eta_pp[COOR(x, d)] * fermion->fermion[fermioncoor->fermionpp[COOR(x, d)]];
        // - (0.5) \eta_mu(x) * U_mu(x-mu) * delta_{x, x-mu}
        target->fermion[x] -= (0.5) * fermioncoor->eta_mm[COOR(x, d)] * fermion->fermion[fermioncoor->fermionmm[COOR(x, d)]]; 
    }
}

/**
 * Compute the product of ( D^\dagger )_{ij} ( F )_{j}.
 * The result is saved into the target fermion that does not require to be initialized to zero when passed to the function.
 */
void AdjDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor)
{
    // mass term for both the upper part of the fermion vector (U) and the bottom part (D)
    // if SIGMA_INTERACTION is defined, consider also the presence of an effective mass term due to sigma
    for(int x=0; x<VOLUME; x++)
    {
        #ifndef SIGMA_INTERACTION
        #ifndef PI_INTERACTION
        target->fermion[x] = 0.0;
        #endif
        #endif

        #ifndef PI_INTERACTION
        #ifdef SIGMA_INTERACTION
        target->fermion[x] = (lat->sigmamass[x]) * fermion->fermion[x];
        #endif
        #endif

        #ifdef SIGMA_INTERACTION
        #ifdef PI_INTERACTION
        target->fermion[x] = (lat->sigmamass[x] - (I_UNIT * lat->epsilon_pi[x] * lat->pimass[x])) * fermion->fermion[x];
        #endif
        #endif
    }
    
    // Out-of-diagonal D terms
    // ( D )^\dagger = - ( D )
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        // + (0.5) \eta_mu(x) * U_mu(x) \delta_{x, x+mu}
        target->fermion[fermioncoor->fermionpp[COOR(x, d)]] += (0.5) * fermioncoor->eta_pp[COOR(x, d)] * fermion->fermion[x];
        // - (0.5) \eta_mu(x) * U_mu(x-mu) * delta_{x, x-mu}
        target->fermion[fermioncoor->fermionmm[COOR(x, d)]] -= (0.5) * fermioncoor->eta_mm[COOR(x, d)] * fermion->fermion[x]; 
    }
}

/**
 * Compute the product of a Dirac matrix and a full fermion vector.
 * The resulting fermion is saved into target.
 * The function exploits the knowledge of only M (the mass matrix) and ( D_oe )
 */
void PBCDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor)
{
    // mass term for both the upper part of the fermion vector (U) and the bottom part (D)
    // the mass term initialize the target fermion
    // if SIGMA_INTERACTION is defined, consider also the presence of an effective mass term due to sigma
    for(int x=0; x<VOLUME; x++)
    {
        #ifndef SIGMA_INTERACTION
        #ifndef PI_INTERACTION
        target->fermion[x] = 0.0;
        #endif
        #endif

        #ifndef PI_INTERACTION
        #ifdef SIGMA_INTERACTION
        target->fermion[x] = (lat->sigmamass[x]) * fermion->fermion[x];
        #endif
        #endif

        #ifdef SIGMA_INTERACTION
        #ifdef PI_INTERACTION
        target->fermion[x] = (lat->sigmamass[x] + (I_UNIT * lat->epsilon_pi[x] * lat->pimass[x])) * fermion->fermion[x];
        #endif
        #endif
    }
    

    // Out-of-diagonal D terms
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        // + (0.5) \eta_mu(x) * U_mu(x) \delta_{x, x+mu}
        target->fermion[x] += (0.5) * fermioncoor->eta_pbc[COOR(x, d)] * fermion->fermion[fermioncoor->fermionpp[COOR(x, d)]];
        // - (0.5) \eta_mu(x) * U_mu(x-mu) * delta_{x, x-mu}
        target->fermion[x] -= (0.5) * fermioncoor->eta_pbc[COOR(x, d)] * fermion->fermion[fermioncoor->fermionmm[COOR(x, d)]]; 
    }
}

/**
 * Compute the product of ( D^\dagger )_{ij} ( F )_{j}.
 * The result is saved into the target fermion that does not require to be initialized to zero when passed to the function.
 */
void PBCAdjDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor)
{
    // mass term for both the upper part of the fermion vector (U) and the bottom part (D)
    // if SIGMA_INTERACTION is defined, consider also the presence of an effective mass term due to sigma
    for(int x=0; x<VOLUME; x++)
    {
        #ifndef SIGMA_INTERACTION
        #ifndef PI_INTERACTION
        target->fermion[x] = 0.0;
        #endif
        #endif

        #ifndef PI_INTERACTION
        #ifdef SIGMA_INTERACTION
        target->fermion[x] = (lat->sigmamass[x]) * fermion->fermion[x];
        #endif
        #endif

        #ifdef SIGMA_INTERACTION
        #ifdef PI_INTERACTION
        target->fermion[x] = (lat->sigmamass[x] - (I_UNIT * lat->epsilon_pi[x] * lat->pimass[x])) * fermion->fermion[x];
        #endif
        #endif
    }
    
    // Out-of-diagonal D terms
    // ( D )^\dagger = - ( D )
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        // + (0.5) \eta_mu(x) * U_mu(x) \delta_{x, x+mu}
        target->fermion[fermioncoor->fermionpp[COOR(x, d)]] += (0.5) * fermioncoor->eta_pbc[COOR(x, d)] * fermion->fermion[x];
        // - (0.5) \eta_mu(x) * U_mu(x-mu) * delta_{x, x-mu}
        target->fermion[fermioncoor->fermionmm[COOR(x, d)]] -= (0.5) * fermioncoor->eta_pbc[COOR(x, d)] * fermion->fermion[x]; 
    }
}

/**
 * Compute the product of ( D )_{ij} ( F )_{j}.
 * The result is saved into the target fermion that does not require to be initialized to zero when passed to the function.
 */
void SubtractDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor)
{
    // mass term for both the upper part of the fermion vector (U) and the bottom part (D)
    // if SIGMA_INTERACTION is defined, consider also the presence of an effective mass term due to sigma
    for(int x=0; x<VOLUME; x++)
    {
        #ifndef SIGMA_INTERACTION
        #ifndef PI_INTERACTION
        target->fermion[x] -= 0.0;
        #endif
        #endif

        #ifndef PI_INTERACTION
        #ifdef SIGMA_INTERACTION
        target->fermion[x] -= (lat->sigmamass[x]) * fermion->fermion[x];
        #endif
        #endif

        #ifdef SIGMA_INTERACTION
        #ifdef PI_INTERACTION
        target->fermion[x] -= (lat->sigmamass[x] + (I_UNIT * lat->epsilon_pi[x] * lat->pimass[x])) * fermion->fermion[x];
        #endif
        #endif
    }
    
    // Out-of-diagonal D terms
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        // + (0.5) \eta_mu(x) * U_mu(x) \delta_{x, x+mu}
        target->fermion[x] -= (0.5) * fermioncoor->eta_pp[COOR(x, d)] * fermion->fermion[fermioncoor->fermionpp[COOR(x, d)]];
        // - (0.5) \eta_mu(x) * U_mu(x-mu) * delta_{x, x-mu}
        target->fermion[x] += (0.5) * fermioncoor->eta_mm[COOR(x, d)] * fermion->fermion[fermioncoor->fermionmm[COOR(x, d)]]; 
    }
}

/**
 * Compute the product of ( D^\dagger )_{ij} ( F )_{j}.
 * The result is saved into the target fermion that does not require to be initialized to zero when passed to the function.
 */
void SubtractAdjDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor)
{
    // mass term for both the upper part of the fermion vector (U) and the bottom part (D)
    // if SIGMA_INTERACTION is defined, consider also the presence of an effective mass term due to sigma
    for(int x=0; x<VOLUME; x++)
    {
        #ifndef SIGMA_INTERACTION
        #ifndef PI_INTERACTION
        target->fermion[x] -= 0.0;
        #endif
        #endif

        #ifndef PI_INTERACTION
        #ifdef SIGMA_INTERACTION
        target->fermion[x] -= (lat->sigmamass[x]) * fermion->fermion[x];
        #endif
        #endif

        #ifdef SIGMA_INTERACTION
        #ifdef PI_INTERACTION
        target->fermion[x] -= (lat->sigmamass[x] - (I_UNIT * lat->epsilon_pi[x] * lat->pimass[x])) * fermion->fermion[x];
        #endif
        #endif
    }
    
    // Out-of-diagonal D terms
    // ( D )^\dagger = - ( D )
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        // + (0.5) \eta_mu(x) * U_mu(x) \delta_{x, x+mu}
        target->fermion[fermioncoor->fermionpp[COOR(x, d)]] -= (0.5) * fermioncoor->eta_pp[COOR(x, d)] * fermion->fermion[x];
        // - (0.5) \eta_mu(x) * U_mu(x-mu) * delta_{x, x-mu}
        target->fermion[fermioncoor->fermionmm[COOR(x, d)]] += (0.5) * fermioncoor->eta_mm[COOR(x, d)] * fermion->fermion[x]; 
    }
}