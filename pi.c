#include "pi.h"

//////////////////////////////////////////////////////////////////////////////////////
//                             pi functions
//////////////////////////////////////////////////////////////////////////////////////

/**
 * Print all the pi fields in ascending order of the fermionic field.
 */
void PrintPi(const double const * pi)
{
    for(int x=0; x<VOLUME; x++)
    printf("Pi[%d] is %.16f\n", x, pi[x]);
}

/**
 * Start the sigma variables according to a specific macro.
 */
void StartPi(struct Lattice * lat, const struct Fermioncoor const * fcor)
{
    #ifdef PI_INTERACTION

    #ifdef START_PI_ZERO
    // start all pi fields to zero
    for(int x=0; x<VOLUME; x++)
    lat->pi[x] = 0.0;
    #endif
    #ifdef START_PI_GAUSSIAN
    for(int x=0; x<HALFVOLUME; x++)
    Gauss2MeanVariance(&(lat->pi[x]), &(lat->pi[x+HALFVOLUME]), 0.0, 1.0);
    // the pi field has been updated -> update pimass
    ComputeFullPiMass(lat, fcor, lat->pimass);
    #endif
    #ifdef START_WEAK_PI
    for(int x=0; x<VOLUME; x++)
    lat->pi[x] = START_WEAK_PI;
    // the pi field has been updated to START_WEAK_PI
    ComputeFullPiMass(lat, fcor, lat->pimass);
    #endif

    #endif
}

/**
 * Compute the effective mass that will appear in the Dirac matrix.
 * =====================================================================
 * Return : (1.0/8.0) \delta_{x,y} \sum_<x,\tilde{x}> \pi_{\tilde{x}}
 * ===================================================================*/
void ComputeFullPiMass(struct Lattice * lat, const struct Fermioncoor const * fcor, double * pimass)
{
    #if DIM == 3
    // for all VOLUME diagonal matrix elements
    // x is in fermionic coordinates
    for(int x=0; x<VOLUME; x++)
    {
        // compute the sum of all the 8 terms
        pimass[x] = lat->pi[x]; //1
        pimass[x] += lat->pi[fcor->fermionmm[COOR(x, 0)]]; //2
        pimass[x] += lat->pi[fcor->fermionmm[COOR(x, 1)]]; //3
        pimass[x] += lat->pi[fcor->fermionmm[COOR(x, 2)]]; //4
        pimass[x] += lat->pi[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 1)]]; //5
        pimass[x] += lat->pi[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 2)]]; //6
        pimass[x] += lat->pi[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 1)], 2)]]; //7
        pimass[x] += lat->pi[fcor->fermionmm[COOR(fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 1)], 2)]]; //8
        pimass[x] = pimass[x] / 8.0;
    }
    #endif
}

/**
 * Return the sum of the squared mass of the pi field.
 * ===============================================================
 * Return : \sum_<\tilde{x}> \pi^2(\tilde{x})
 * =============================================================*/
double ComputeTotalPiSquared(const struct Lattice * lat)
{
    double total = 0.0;
    // for all VOLUME sites
    for(int x=0; x<VOLUME; x++)
    total += pow(lat->pi[x], 2);

    return total;
}

/**
 * Return the expectation value of the condensate evaluated by means of the expectation value of the pi field.
 * =======================================================================================================
 * Return : (1.0 / VOLUME) \sum_{\tilde{x}=1}^{VOLUME} < \sigma(\tilde{x}) >
 * =====================================================================================================*/
double CondensatePi(const struct Lattice * lat, const struct Info const * info)
{
    double total = 0.0;

    #ifdef PI_INTERACTION
    for(int x=0; x<VOLUME; x++)
    total += lat->pi[x];
    #endif

    return (total / VOLUME);
}

/**
 * Copy all the VOLUME elements from pi vector to a copy
 */
void CopyFromPi(const double * from, double * target)
{
    // copy all the VOLUME elements into copy
    for(int x=0; x<VOLUME; x++)
    target[x] = from[x];
}

/**
 * Copy all the VOLUME elements from the copy into pi
 */
void CopyToPi(const double * from, double * target)
{
    // copy all the VOLUME elements from the copy into pi
    for(int x=0; x<VOLUME; x++)
    target[x] = from[x];
}