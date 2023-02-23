#ifndef FERMION_H
#define FERMION_H

#include "macro.h"
#include "coordinate.h"
#include "random.h"
#include "lattice.h"

/**==================================================================
 * Structure containing the fermion sites in ascending order (ordered in fermion Lessicographic coordinates).
 * ================================================================*/
typedef struct Fermioncoor
{
    // Site in Lessicographic coordinates
    int * site;
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

}Fermioncoor;

// Fermioncoor functions
void AllocateFermioncoor(struct Fermioncoor * fermioncoor);
void StartFermioncoor(struct Fermioncoor * fermioncoor);
void StartFermioncoorFromLattice(struct Fermioncoor * fcor, const struct Lattice const * lat);
void PrintFermioncoorSite(const struct Fermioncoor const * fermioncoor, int site);
void PrintFermioncoor(const struct Fermioncoor const * fermioncoor);
void DestroyFermioncoor(struct Fermioncoor * fermioncoor);


/**===================================
 * Fermionic vector.
 * In the case of staggered fermions, only the fermion fields associated with even (or odd) 
 * sites are really necessary. By the way, all the VOLUME components of the fields are stored.
 *====================================*/
typedef struct Fermion
{
    complex * fermion;
}Fermion;

// fermion function 
void AllocateFermion(struct Fermion * fermion);
void DestroyFermion(struct Fermion * fermion);
inline void SetHalfFermionToZero(struct Fermion * fermion, int start){for(int i=start; i<start+HALFVOLUME; i++) fermion->fermion[i] = 0.0;};
inline void SetFermionToZero(struct Fermion * fermion){for(int r=0; r<VOLUME; r++) fermion->fermion[r] = 0.0;};
inline void CopyHalfFermion(const struct Fermion const * from, struct Fermion * target, int start){for(int i=start; i<start+HALFVOLUME; i++) target->fermion[i] = from->fermion[i];};
inline void CopyFermion(const struct Fermion const * from, struct Fermion * target){for(int x=0; x<VOLUME; x++) target->fermion[x] = from->fermion[x];};
inline void AddHalfFermionWithCoefficientToTarget(const struct Fermion const * fermion, complex coeff, struct Fermion * target, int start){for(int x=start; x<start+HALFVOLUME; x++) target->fermion[x] += (coeff * fermion->fermion[x]);};
inline void AddFermionWithCoefficientToTarget(const struct Fermion const * fermion, complex coeff, struct Fermion * target){for(int x=0; x<VOLUME; x++) target->fermion[x] += (coeff * fermion->fermion[x]);};
void PrintFermion(const struct Fermion const * fermion);
double NormHalfFermion(const struct Fermion const * fermion, int start);
double NormFullFermion(const struct Fermion const * fermion);
double NormUpperHalfFermion(const struct Fermion const * fermion);
inline double NormFermion(const struct Fermion const * fermion){double total = 0.0; for(int x=0; x<VOLUME; x++) total += creal(conj(fermion->fermion[x]) * fermion->fermion[x]); return pow(total, 0.5);};

// fermion for hybrid
void DrawHalfFermionsFromGaussianDistribution(struct Fermion * fermion);
void DrawFullFermionFromGaussianDistribution(struct Fermion * fermion);

/**===========================================
 * ManyFermions.
 * This is a collection of HALFNFLAV fermions. 
 * ==========================================*/
typedef struct ManyFermions
{
    struct Fermion * flav;
}ManyFermions;

// many fermions function
void AllocateManyFermions(struct ManyFermions * fermions);
inline void SetToZeroManyFermions(struct ManyFermions * fermions){for(int f=0; f<HALFNFLAV; f++)SetFermionToZero(&(fermions->flav[f]));};
void PrintFermions(const struct ManyFermions * fermions);
void DestroyManyFermions(struct ManyFermions * fermions);

/**====================================
 * Fermion kit.
 * This is a structure containing several fermion vectors
 * useful for the conjugate gradient method.
 * ==================================*/
typedef struct Fermionkit
{
    // fermion variables
    struct ManyFermions r_k;
    struct ManyFermions p_k;
    struct ManyFermions chi;
    struct ManyFermions phi;
    struct ManyFermions inverted_phi;
    struct ManyFermions inverted_phi_obs;
    struct Fermion Ddag_inverted_phi_obs;
    struct Fermion temp1;
    struct Fermion temp2;
    struct Fermion tempbminusx0;
    struct Fermion bminusax0;

    // CG info
    double iterations;
    int CG_call;
    
}Fermionkit;

// fermionkit functions

void AllocateFermionkit(struct Fermionkit * fermionkit);
void DestroyFermionkit(struct Fermionkit * fermionkit);
#endif