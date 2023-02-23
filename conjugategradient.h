#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "lattice.h"
#include "dirac.h"

// conjugate gradient

bool ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target_inverted, const struct Info const * info, const struct Fermioncoor const * fermioncoor, double inversion_precision, struct Fermionkit * kit, int f);
bool StartX0ImposedConjugateGradientFastDiracAdjFastDiracToTheMinusOne(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target_inverted, const struct Info const * info, const struct Fermioncoor const * fermioncoor, double inversion_precision, struct Fermionkit * kit, int f);
bool PBCConjugateGradientFastDiracAdjFastDiracToTheMinusOne(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target_inverted, const struct Info const * info, const struct Fermioncoor const * fermioncoor, double inversion_precision, struct Fermionkit * kit, int f);
double BMinusAX0(const struct Lattice const * lat, const struct Fermion const * B, const struct Fermion const * x0, const struct Info const * info, const struct Fermioncoor const * fcor, struct Fermionkit * kit);
#endif