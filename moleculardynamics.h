#ifndef MOLECULARDYNAMICS_H
#define MOLECULARDYNAMICS_H

#include "macro.h"
#include "random.h"
#include "lattice.h"
#include "linklattice.h"
#include "conjugategradient.h"
#include "sigma.h"
#include "pi.h"

/**==========================================================
 * Struct containing the momenta related to the gauge fields.
 * The momenta are stored into fermionic coordinates.
 * ==========================================================
 */
typedef struct Momenta
{
    double * momenta;
}Momenta;

/**==================================================
 * Struct containing the forces for the gauge fields.
 * The forces are stored into fermionic coordiantes.
 * ==================================================
 */
typedef struct Force
{
    double * force;
}Force;

// momenta functions
void AllocateMomenta(struct Momenta * momenta);
void PrintMomenta(const struct Momenta const * momenta);
void DrawMomentaFromGaussianDistribution(struct Momenta * momenta);
void DestroyMomenta(struct Momenta * momenta);

// force functions
void AllocateForce(struct Force * force);
inline void SetForceToZero(struct Force * force){for(int r=0; r<DIM*VOLUME; r++) force->force[r] = 0.0;};
void PrintForce(const struct Force const * force);
void DestroyForce(struct Force * force);
double NormForce(const struct Force const * force);

// force computation
void AddForceFromPureGaugeSector(const struct Lattice const * lattice, const struct Info const * info, struct Force * force);
void PartialAdjDiracFermion(const struct Lattice const * lat, const struct Fermion const * start, struct Fermion * target, int fermion_x, int mu, const struct Fermioncoor const * fermioncoor);
bool AddForceFromFermionSectorRespectToLink(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Force * force, struct Fermionkit * kit);
bool ComputeForce(const struct Lattice const * lattice, const struct Info const * info, struct Force * force, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit);

/**========================================================
 * typedef struct Momentasigma
 * the structure contains the momenta associated with the sigma field
 * ======================================================*/
typedef struct Momentasigma
{
    double * momenta;
}Momentasigma;

/**========================================================
 * typedef struct Momentapi
 * the structure contains the momenta associated with the pi field
 * ======================================================*/
typedef struct Momentapi
{
    double * momenta;
}Momentapi;

// energy functions
double ComputeTotalFictitiousEnergyWithChi(const struct Lattice const * lattice, const struct Info const * info, const struct Momenta const * momenta, const struct Fermionkit const * kit);
double ComputeFullEnergyWithPhiWithoutSigma(const struct Lattice const * lattice, const struct Info const * info, const struct Fermioncoor const * fcor, const struct Momenta const * momenta, const struct Fermionkit const * kit);
double ComputeFullEnergyForSigma(const struct Lattice const * lattice, const struct Info const * info, const struct Fermioncoor const * fcor, const struct Momentasigma const * momenta, const struct Fermionkit const * kit);
double ComputeFullEnergyForSigmaWithChi(const struct Lattice const * lattice, const struct Info const * info, const struct Fermioncoor const * fcor, const struct Momentasigma const * momenta, const struct Fermionkit const * kit);
double ComputeFullEnergyForPi(const struct Lattice const * lattice, const struct Info const * info, const struct Fermioncoor const * fcor, const struct Momentapi const * momenta, const struct Fermionkit const * kit);

// evolution functions
void EvolveCoordinates(struct Lattice * lattice, const struct Fermioncoor const * fcor, const struct Momenta const * momenta, double epsilon);
void EvolveMomenta(const struct Lattice const * lattice, const struct Force const * force, struct Momenta * momenta, double epsilon);
bool LeapFrog2Link(struct Lattice * lattice, struct Info * info, struct Force * force, struct Momenta * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file);

/**=======================================================
 * typedef struct Forcesigma.
 * This struct contains the forces of the sigma field
 * =====================================================*/
typedef struct Forcesigma
{
    double * force;
}Forcesigma;

/**=======================================================
 * typedef struct Forcepi.
 * This struct contains the forces of the pi field
 * =====================================================*/
typedef struct Forcepi
{
    double * force;
}Forcepi;

// force functions (sigma)
void AllocateForceSigma(struct Forcesigma * force);
inline void SetForceSigmaToZero(struct Forcesigma * force){for(int x=0; x<VOLUME; x++) force->force[x] = 0.0;};
void DestroyForceSigma(struct Forcesigma * force);
void AddForceFromSigmaSector(const struct Lattice const * lat, const struct Info const * info, struct Forcesigma * force);
bool AddForceFromFermionSectorRespectToSigma(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Forcesigma * force, struct Fermionkit * kit, double inversion_precision);
bool ComputeSigmaForce(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Forcesigma * force, struct Fermionkit * kit, double inversion_precision);

// force functions (pi)
void AllocateForcePi(struct Forcepi * force);
inline void SetForcePiToZero(struct Forcepi * force){for(int x=0; x<VOLUME; x++) force->force[x] = 0.0;};
void DestroyForcePi(struct Forcepi * force);
void AddForceFromPiSector(const struct Lattice const * lat, const struct Info const * info, struct Forcepi * force);
bool AddForceFromFermionSectorRespectToPi(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Forcepi * force, struct Fermionkit * kit);
bool ComputePiForce(const struct Lattice const * lat, const struct Info const * info, const struct Fermioncoor const * fcor, struct Forcepi * force, struct Fermionkit * kit);


// momenta functions (sigma)
void AllocateMomentaSigma(struct Momentasigma * momenta);
void DestroyMomentaSigma(struct Momentasigma * momenta);
void DrawSigmaMomentaFromGaussianDistribution(struct Momentasigma * momenta);
// momenta functions (pi)
void AllocateMomentaPi(struct Momentapi * momenta);
void DestroyMomentaPi(struct Momentapi * momenta);
void DrawPiMomentaFromGaussianDistribution(struct Momentapi * momenta);

// evolution sigma functions
void EvolveCoordinatesSigma(struct Lattice * lat, const struct Momentasigma const * momenta, const struct Fermioncoor const * fcor, double epsilon);
void EvolveMomentaSigma(struct Momentasigma * momenta, const struct Forcesigma * force, double epsilon);
bool LeapFrog2Sigma(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file);
bool LeapFrog2MinimumNormSigma(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file);
bool LeapFrog2MinimumNormSigmaVelocityFirst(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file);
bool LeapFrog4MinimumNorm4SigmaFirstPosition(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file);
void EvolveSigmaInTimeLF2(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step);
void EvolveSigmaInTimeMinimumNormLF2(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step);
void EvolveSigmaInTimeMinimumNormSigmaVelocityFirst(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step);
void FreeLeapFrog2SigmaWithoutAcceptance(struct Lattice * lattice, struct Info * info, struct Forcesigma * force, struct Momentasigma * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step);

// evolution pi functions
void EvolveCoordinatesPi(struct Lattice * lat, const struct Momentapi const * momenta, const struct Fermioncoor const * fcor, double epsilon);
void EvolveMomentaPi(struct Momentapi * momenta, const struct Forcepi * force, double epsilon);
bool LeapFrog2Pi(struct Lattice * lattice, struct Info * info, struct Forcepi * force, struct Momentapi * momenta, const struct Fermioncoor const * fermioncoor, struct Fermionkit * kit, double epsilon, const int step, const char const * acceptance_file);

#endif