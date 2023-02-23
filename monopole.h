#ifndef MONOPOLE_H
#define MONOPOLE_H

#include "lattice.h"

// monopole functions
double MFunction(double x);
double MFunctionSurface(const struct Lattice const * lat, int x, int mu, int nu);
double AbsMonopoleSite(const struct Lattice const * lat, int x);
double MonopoleDensity(const struct Lattice const * lat);

#endif