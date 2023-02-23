#ifndef MACRO_H
#define MACRO_H


#include<limits.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<complex.h>
#include<time.h>
#include<stdbool.h>
#include<time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/////////////////////////////////////////////////////////////
//           Lattice parameters                            
/////////////////////////////////////////////////////////////
#define LSPACE (8)
#define LTIME  (LSPACE)
#define DIM (3)
#define VOLUME (LSPACE * LSPACE * LTIME)
#define NFLAV (4)
#define Q_CHARGE (1)


#define HALFNFLAV (NFLAV/2)
#define HALFVOLUME  (VOLUME/2)


#define LATTICE_UPDATE (1)
#define IMEAS (100)
#define MEASURE (10000)
/////////////////////////////////////////////////////////////
//                  Coordinates
/////////////////////////////////////////////////////////////
#define COOR(x, d)  ((x) * DIM + (d))

/////////////////////////////////////////////////////////////////////
//                    Hamiltonian
/////////////////////////////////////////////////////////////////////
/**
 * PURE_GAUGE : if defined it means that the pure gauge action is present.
 */
//#define PURE_GAUGE
#define PSEUDO_FERMION
#define SIGMA_INTERACTION
//#define PI_INTERACTION
////////////////////////////////////////////////////////////////////
//               Initialization parameters
////////////////////////////////////////////////////////////////////
////////               Link
/**
 * START_LINK_IDENTITY : Links (the group elements) are initialized to the identity.
 * START_LINK_RANDOM : Links are initialized randomly in (-PI, PI).
 */
#define START_LINK_IDENTITY
//#define START_LINK_RANDOM

//////////////          Sigma           ////////////////////
//#define START_SIGMA_ZERO
//#define START_SIGMA_GAUSSIAN
#define START_WEAK_SIGMA (0.3)

// ---------- Integrator ----------------------------

//#define LF2SIGMA
//#define LFMN2SIGMA
#define LFMN2SIGMAVF
//#define LF4MN4FPSIGMA

//////////////          Pi        ///////////////////////
//#define START_PI_ZERO
#define START_PI_GAUSSIAN
//#define START_WEAK_PI (0.3)

////////////////////////////////////////////////////////////////////
//                    Observables
////////////////////////////////////////////////////////////////////
//#define OBS_MONOPOLE
#define OBS_CONDENSATE
#define OBS_CONDENSATE_SIGMA
//#define OBS_PLAQUETTE

#define CONDENSATE_PRECISION (1)
///////////////////////////////////////////////////////////////////
//                   Boundary conditions
///////////////////////////////////////////////////////////////////
/**
 * If any macro associated with the boundary conditions is defined,
 * boundary conditions are defined PERIODIC for all fields in all directions.
 * -------------------------------------------------------------------
 * #define APBC_FERMION_TIME : the fermion field is defined with APBC along time direction.
 *                             The BC is encapsulated into the eta_pp arrays.      
 * #define APBC_ALL : the fermion field has APBC along all lattice directions.
 *                      The BCs are encapsulated in the eta_pp and eta_mm arrays.                       
 */
#define APBC_FERMION_TIME
//#define APBC_ALL

///////////////////////////////////////////////////////////////////
//                    Info parameters
//////////////////////////////////////////////////////////////////
/**
 * METROPOLIS_LINK : It is defined if a metropolis update for the link variables is considered.
 * METROPOLIS_LINK_THRESHOLD : If the number of metropolis_link_done is larger than this threshold
 *                             the Metropolis link info are updated.
 * AUTOREGULATION : if autoregulation is defined, the HMC parameters epsilon and step autoregulate theireslves 
 *                  consistently in order to avoid 2 rejected HMC consecutively and \order 10 accepted HMC consecutively
 */
//#define METROPOLIS_LINK
#define HYBRID_MC
//#define AUTOREGULATION
#define TERMALIZATION_STEP (10000)

#define METROPOLIS_LINK_THRESHOLD (100000)
#define METROPOLIS_LINK_ACCEPTANCE (0.3)
#define PRINT_LINK_ACCEPTANCE
////////////////////////////////////////////////////////////////
//                 Conjugate gradient
/////////////////////////////////////////////////////////////////
/**
 * START_X0_NULL : set the starting vector for the inversion to zero.
 * START_X0_PREVIOUS : leave the starting vector as it is
 * MAX_CONJUGATE_GRADIENT_ITERATION : define the maximum number of iterations in the conjugate gradient method
 */
#define START_X0_NULL
//#define START_X0_PREVIOUS

// this is the maximum number of iterations allowed to invert a pure EVEN fermion field
// for the inversion of the full fermion field this number is doubled
#define MAX_CONJUGATE_GRADIENT_ITERATION (VOLUME + 1)
#define MEDIUM_PRECISION_INVERSION (1e-3)
#define HIGH_PRECISION_INVERSION (1e-9)

//////////////////////////////////////////////////////////////
//                Condensate
//////////////////////////////////////////////////////////////
/**
 * Z2_NOISE : the RandomVec is initialized with Z2 noise
 * PHASE_NOISE : the RandomVec is initialized with complex phases exp(i theta)
 */
#define Z2_NOISE
//#define PHASE_NOISE

//////////////////////////////////////////////////////////////
//                  Debug
//////////////////////////////////////////////////////////////
/**
 * DEBUG_LINK : The debugger is active only for functions associated with the link variables.
 * DEBUG_COORDINATE : The debgger is active only for functions associated with different coordinates 
 *                      and the lattice geometry.
 */
//#define DEBUG_LINK
//#define DEBUG_COORDINATE
//#define DEBUG_CONJUGATE_GRADIENT
//#define DEBUG_ADJ_FERMION_FAST_DIRAC
//#define DEBUG_LF2_SIGMA
//#define DEBUG_STARTING_X0
//#define DEBUG_PRINT_ITERATIONS
//#define DEBUG_STARTING_ENERGY_HMC
///////////////////////////////////////////////////////////////////
//             Constants
///////////////////////////////////////////////////////////////////
#define PI (3.141592653589793238462643383279502884197169399375105820974944)
// LF2MN constant
#define MAGIC_LAMBDA (0.1931833275037836)
// LF4MN4FPSIGMA constants
#define LF4MN4_RHO  (0.1786178958448091)
#define LF4MN4_THETA (-0.06626458266981843)
#define LF4MN4_LAMBDA (0.7123418310626056)

#define I_UNIT (I)

// define the minimum momentum always equal to (2 pi / L)
#define P_MIN ((2 * PI) / LSPACE)
#define P_MIN_APBC ((PI) / LSPACE)

#endif
