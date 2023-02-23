#ifndef CONDENSATE_H
#define CONDENSATE_H

#include "conjugategradient.h"
#include "inputoutput.h"

/**==============================================================
 * Struct RandomVec.
 * These are random vector numbers (complex), whose length is given by VOLUME.
 *=============================================================*/
typedef struct RandomVec
{
    struct Fermion * vec;

}RandomVec;

// random vector functions
void AllocateRandomVec(struct RandomVec * vec);
void DestroyRandomVec(struct RandomVec * vec);
void SetFullRandomVecToZero(struct RandomVec * vec);
void Z2NoiseHalfRandomVec(struct RandomVec * vec, int start);
void Z2NoiseFullRandomVec(struct RandomVec * vec);
void FullVecOneOverVolume(struct RandomVec * vec);
void FullVecAllPhaseOverVolume(struct RandomVec * vec, const struct Lattice const * lattice);
void FullVecPhasesForChipOverVolume(struct RandomVec * vec, int d, const struct Lattice const * lat);
void PhaseFullNoise(struct RandomVec * vec);
void ComputeCondensateDirac(const struct Lattice const * lat, double * condensate_tr1, double * condensate_tr1_abs, double * condensate_tr2, double * condensate_tr4, double * condensate_mu2, double * condensate_mu4, double * condensate_chi, double * condensate_chi_abs, double * condensate_chip, double * condensate12_nosigma, double * condensate12_withsigma, struct RandomVec * vec, struct Fermion * target_inverted, const struct Info const * info, const struct Fermioncoor const * fermioncoor, double inversion_precision, struct Fermionkit * kit, int f);
#endif