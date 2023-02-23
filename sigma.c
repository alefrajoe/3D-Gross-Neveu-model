#include "sigma.h"

//////////////////////////////////////////////////////////////////////////////////////
//                             sigma functions
//////////////////////////////////////////////////////////////////////////////////////

/**
 * Print all the sigma fields in ascending order of the fermionic field.
 */
void PrintSigma(const double const * sigma)
{
    for(int x=0; x<VOLUME; x++)
    printf("Sigma[%d] is %.16f\n", x, sigma[x]);
}

/**
 * Start the sigma variables according to a specific macro.
 */
void StartSigma(struct Lattice * lat, const struct Fermioncoor const * fcor)
{
    #ifdef SIGMA_INTERACTION

    #ifdef START_SIGMA_ZERO
    // start all sigma fields to zero
    for(int x=0; x<VOLUME; x++)
    lat->sigma[x] = 0.0;
    #endif
    #ifdef START_SIGMA_GAUSSIAN
    for(int x=0; x<HALFVOLUME; x++)
    Gauss2MeanVariance(&(lat->sigma[x]), &(lat->sigma[x+HALFVOLUME]), 0.0, 1.0);
    // the sigma field has been updated -> update sigmamass
    ComputeFullSigmaMass(lat, fcor, lat->sigmamass);
    #endif
    #ifdef START_WEAK_SIGMA
    for(int x=0; x<VOLUME; x++)
    lat->sigma[x] = START_WEAK_SIGMA;
    // the sigma field has been updated to START_WEAK_SIGMA
    ComputeFullSigmaMass(lat, fcor, lat->sigmamass);
    #endif

    #endif
}

/**
 * Compute half of matrix with respect to the effective mass that will appear in the Dirac matrix.
 * =====================================================================
 * Return : (1.0/8.0) \delta_{x,y} \sum_<x,\tilde{x}> \sigma_{\tilde{x}}
 * ===================================================================*/
void ComputeHalfSigmaMass(struct Lattice * lat, const struct Fermioncoor const * fcor, double * sigmamass, int start)
{
    #if DIM == 3
    // for all VOLUME diagonal matrix elements
    // x is in fermionic coordinates
    for(int x=start; x<start+HALFVOLUME; x++)
    {
        // compute the sum of all the 8 terms
        sigmamass[x] = lat->sigma[x]; //1
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(x, 0)]]; //2
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(x, 1)]]; //3
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(x, 2)]]; //4
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 1)]]; //5
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 2)]]; //6
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 1)], 2)]]; //7
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 1)], 2)]]; //8
        sigmamass[x] = sigmamass[x] / 8.0;
    }
    #endif
}

/**
 * Compute the effective mass that will appear in the Dirac matrix.
 * =====================================================================
 * Return : (1.0/8.0) \delta_{x,y} \sum_<x,\tilde{x}> \sigma_{\tilde{x}}
 * ===================================================================*/
void ComputeFullSigmaMass(struct Lattice * lat, const struct Fermioncoor const * fcor, double * sigmamass)
{
    #if DIM == 3
    // for all VOLUME diagonal matrix elements
    // x is in fermionic coordinates
    for(int x=0; x<VOLUME; x++)
    {
        // compute the sum of all the 8 terms
        sigmamass[x] = lat->sigma[x]; //1
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(x, 0)]]; //2
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(x, 1)]]; //3
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(x, 2)]]; //4
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 1)]]; //5
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 2)]]; //6
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 1)], 2)]]; //7
        sigmamass[x] += lat->sigma[fcor->fermionmm[COOR(fcor->fermionmm[COOR(fcor->fermionmm[COOR(x, 0)], 1)], 2)]]; //8
        sigmamass[x] = sigmamass[x] / 8.0;
    }
    #endif
}

/**
 * Return the sum of the squared mass of the sigma field.
 * ===============================================================
 * Return : \sum_<\tilde{x}> \sigma^2(\tilde{x})
 * =============================================================*/
double ComputeTotalSigmaSquared(const struct Lattice * lat)
{
    double total = 0.0;
    // for all VOLUME sites
    for(int x=0; x<VOLUME; x++)
    total += pow(lat->sigma[x], 2);

    return total;
}

/**
 * Copy all the VOLUME elements from sigma vector to a copy
 */
void CopyFromSigma(const double * from, double * target)
{
    // copy all the VOLUME elements into copy
    for(int x=0; x<VOLUME; x++)
    target[x] = from[x];
}

/**
 * Copy all the VOLUME elements from the copy into sigma
 */
void CopyToSigma(const double * from, double * target)
{
    // copy all the VOLUME elements from the copy into sigma
    for(int x=0; x<VOLUME; x++)
    target[x] = from[x];
}