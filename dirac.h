#ifndef DIRAC_H
#define DIRAC_H

#include "fermion.h"
#include "info.h"
#include "lattice.h"

// Fast Dirac functions
void PrintDirac(const struct Lattice const * lat, const struct Fermioncoor const * fermioncoor, const struct Info const * info);
void DiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor);
void AdjDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor);
void PBCDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor);
void PBCAdjDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor);
void SubtractDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor);
void SubtractAdjDiracFermion(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target, const struct Info const * info, const struct Fermioncoor const * fermioncoor);
#endif