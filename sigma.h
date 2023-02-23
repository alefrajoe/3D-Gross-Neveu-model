#ifndef SIGMA_H
#define SIGMA_H

#include "lattice.h"
#include "dirac.h"

// sigma functions
void PrintSigma(const double const * sigma);
void StartSigma(struct Lattice * lat, const struct Fermioncoor const * fcor);
void ComputeHalfSigmaMass(struct Lattice * lat, const struct Fermioncoor const * fcor, double * sigmamass, int start);
void ComputeFullSigmaMass(struct Lattice * lat, const struct Fermioncoor const * fcor, double * sigmamass);
double ComputeTotalSigmaSquared(const struct Lattice * lat);
void CopyFromSigma(const double * from, double * target);
void CopyToSigma(const double * from, double * target);
#endif