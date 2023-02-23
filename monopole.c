#include "monopole.h"

////////////////////////////////////////////////////////////
//                monopole functions
////////////////////////////////////////////////////////////
/**
 * Given a double number x, the function returns
 * -----------------------------------------------------
 *    m(x) = x - |_ x + 0.5 _|
 * ----------------------------------------------------*/
double MFunction(double x)
{
    // return the result of the M function, exploiting the definition of the floor function
    return x - (floor((x + 0.5)));
}

/**
 * The function computes the m(x) function applied to a surface.
 * Note that if mu-nu are passed in cyclic order, the surface at site x points always to
 * the center of the cube of the dual site associated with x.
 * Here x is defined as follows
 * ---------------------------------------------------------------------
 *   x = theta_mu(x) + theta_nu(x+mu) - theta_mu(x+nu) - theta_nu(x)
 * --------------------------------------------------------------------
 */
double MFunctionSurface(const struct Lattice const * lat, int x, int mu, int nu)
{
    double angle = 0.0;

    // compute the sum of the angles
    // remember that the carg function always return an angle in [-pi, pi)
    angle = carg(lat->dlink[COOR(x, mu)].link) + carg(lat->dlink[COOR(lat->pp[COOR(x, mu)], nu)].link);
    angle -= (carg(lat->dlink[COOR(lat->pp[COOR(x, nu)], mu)].link) + carg(lat->dlink[COOR(x, nu)].link));
    
    return MFunction(angle / (2.0 * PI));
}

/**
 * The function returns the number of monopoles-antimonopoles associated with the upper-right corner of the site 
 * that is dual to x.
 * The function is defined only in DIM == 3, otherwise it returns back 0.
 * It computes the monopoles pointing outward with respect of the cube of the dual lattice.
 */
double AbsMonopoleSite(const struct Lattice const * lat, int x)
{
    #if DIM == 3
    double monopole = 0.0;

    // surface 0-1
    monopole += MFunctionSurface(lat, x, 0, 1);
    monopole -= MFunctionSurface(lat, lat->pp[COOR(x, 2)], 0, 1);
    // surface 1-2
    monopole += MFunctionSurface(lat, x, 1, 2);
    monopole -= MFunctionSurface(lat, lat->pp[COOR(x, 0)], 1, 2);
    //surface 2-0
    monopole += MFunctionSurface(lat, x, 2, 0);
    monopole -= MFunctionSurface(lat, lat->pp[COOR(x, 1)], 2, 0);

    return fabs(monopole);
    #endif

    // return 0 if DIM != 3
    return 0;
}

/**
 * Return the monopole density associated with the lattice.
 * Note that a 2.0 factor is present as we only count the number of M+ (M-) monopoles.
 * --------------------------------------------------------------
 * Return : sum_x MonopoleSite(x) / (2.0 * VOLUME)
 * ------------------------------------------------------------*/
double MonopoleDensity(const struct Lattice const * lat)
{
    double density = 0.0;
    
    // compute the total number of monopoles into the lattice
    for(int x=0; x<VOLUME; x++) density += AbsMonopoleSite(lat, x);

    return density / (2.0 * VOLUME);
}
