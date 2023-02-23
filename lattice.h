#ifndef LATTICE_H
#define LATTICE_H

#include "link.h"
#include "coordinate.h"

/** =====================================================
 * Lattice structure.
 * A lattice is composed of VOLUME sites.
 * A site would only made of DIM links.
 =======================================================*/
typedef struct Lattice
{
    // One link for each direction
    struct Link * dlink;
    // Define also the sigma field, which lives in the dual lattice
    // this variable is always defined, but it is used just whether the macro SIGMA_INTERACTION is defined
    // The sigma field is the unique field saved in FERMIONIC coordinates
    double * sigma;
    double * sigmamass;
    // define the pi field
    double * pi;
    double * epsilon_pi;
    double * pimass;

    // Site in Lessicographic coordinates
    int * site;
    // Site in lattice coordinates [1]..[DIM]. This is a VOLUME * DIM dimensional array, requiring the use of COOR(x, d)
    int * lattice_site;
    // Nearest neighbors in the forward direction (Lessicographic coordinates)
    int * pp;
    // Nearest neighbors in the backward direction (Lessicographic coordinates)
    int * mm;
    // Site in fermionic coordinate (this is X = parity * V/2 + x_stereo / 2)
    int * fermionsite;
    // Nearest neighbors in the forward direction (Lessicographic fermion coordinates)
    int * fermionpp;
    // Nearest neighbors in the backward direction (Lessicographic fermion coordinates)
    int * fermionmm;
    // \eta_mu(x) is the usual phase related to staggered fermions: \eta_\mu(x) = (-)^{n_1 + .. + n_{\mu-1}}
    double * eta_pp;
    // \eta_mu(x) is the usual phase related to staggered fermions: \eta_\mu(x) = (-)^{n_1 + .. + n_{\mu-1}}
    double * eta_mm;
    // \eta_mu(x) is the usual phase related to staggered fermions: \eta_\mu(x) = (-)^{n_1 + .. + n_{\mu-1}}
    double * eta_pbc;

    // distances  and done HMCs
    double distance_link;
    double distance_sigma;
    double distance_pi;

    int done_link;
    int done_sigma;
    int done_pi;

}Lattice;

// sites functions
void StartSite(struct Lattice * lattice, int site_number);
void PrintSite(const struct Lattice const * lattice, int site);

// lattice functions
void AllocateLattice(struct Lattice * lattice);
void StartLattice(struct Lattice * lattice);
void DestroyLattice(struct Lattice * lattice);
void PrintLattice(const struct Lattice const * lattice);
#endif