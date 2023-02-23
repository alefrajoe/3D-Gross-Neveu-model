#ifndef POLYAKOV_H
#define POLYAKOV_H

#include "linklattice.h"

// polyakov functions
double PolyakovSite(const struct Lattice const * lat, int x);
double PolyakovValue(const struct Lattice const * lat);

#endif