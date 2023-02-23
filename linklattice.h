#ifndef LINKLATTICE_H
#define LINKLATTICE_H

#include "lattice.h"
#include "info.h"


// link energy
double RealTracePlaquette(const struct Lattice const * lat, int site, int mu, int nu);
double ImagTracePlaquette(const struct Lattice const * lat, int site, int mu, int nu);
double TotalRealTracePlaquette(const struct Lattice const * lat);
double TotalImagTracePlaquette(const struct Lattice const * lat);
void ComputeStaple(const struct Lattice const * lat, int site, int mu, struct Link * staple);
double ComputeImagAllPlaquettesIncludingLink(const struct Lattice const * lat, int site, int mu);
inline double RealTraceLinkStaple(const struct Link const * link, const struct Link const * staple){return creal(link->link * staple->link);};
inline double ImagTraceLinkStaple(const struct Link const * link, const struct Link const * staple){return cimag(link->link * staple->link);};
double RealTraceLinkSiteStapleSameSite(const struct Lattice const * lat, int site, int mu);
void ObsPlaquette(const struct Lattice const * lat, double * plaq1, double * plaq2, double * plaq3);
void CopyAllLinkFromLattice(const struct Lattice * lat, complex * copy);
void CopyAllLinkToLattice(const complex * copy, struct Lattice * lat);

// metropolis link
void MetropolisLink(struct Lattice * lat, struct Info * info, int site);
void CorrectLinkLattice(struct Lattice * lat);
#endif