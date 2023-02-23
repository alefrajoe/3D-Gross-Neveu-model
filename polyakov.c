#include "polyakov.h"

//////////////////////////////////////////////////////////////////
//                 polyakov functions
//////////////////////////////////////////////////////////////////

/**
 * Compute the value of the Polyakov loop associated with a site of the lattice.
 */
double PolyakovSite(const struct Lattice const * lat, int x)
{
    double poly = 0.0;
    // start temp from 1.0 because at step zero it is multiplied by the link U_mu(x)
    complex link = 1.0;
    int temp_pos = x;
    // for all sites in the temporal direction
    for(int t=0; t<LTIME; t++)
    {
        // multiply the "next" link
        poly *= lat->dlink[COOR(temp_pos, DIM-1)].link;
        // update the lattice position
        temp_pos = lat->pp[COOR(temp_pos, DIM-1)];
    }

    return creal(poly);
}

/**
 * The function returns the expectation value of the Polyakov loop, evaluated for all LSPACE^{DIM-1} lattice points
 * This is the order parameter of center symmetry.
 * Fix the temporal slice to zero and then compute all the Polyakov
 * ===========================================================================
 * < Poly > = 1/LSPACE^{DIM-1}   \sum_x^{LSPACE^{DIM-1}}  RE \Tr { \Prod_{i=1}^LTIME U_t(x+i t) }
 * ===========================================================================
 */
double PolyakovValue(const struct Lattice const * lat)
{
    double poly = 0.0;

    // for all sites in VOLUME
    for(int x=0; x<VOLUME; x++)
    {
        // if the temporal coordinate of the lattice is 0
        // the zero temporal slice is always present in all simulations
        if((x % LTIME) == 0)
        poly += (PolyakovSite(lat, x) / (VOLUME)) * LTIME;
    }
    return poly;
}
