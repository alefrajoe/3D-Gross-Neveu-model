#ifndef PI_H
#define PI_H

#include "lattice.h"
#include "dirac.h"

// pi functions
void PrintPi(const double const * pi);
void StartPi(struct Lattice * lat, const struct Fermioncoor const * fcor);
void ComputeFullPiMass(struct Lattice * lat, const struct Fermioncoor const * fcor, double * pimass);
double ComputeTotalPiSquared(const struct Lattice * lat);
double CondensatePi(const struct Lattice * lat, const struct Info const * info);
void CopyFromPi(const double * from, double * target);
void CopyToPi(const double * from, double * target);
#endif