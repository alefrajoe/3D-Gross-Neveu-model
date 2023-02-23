#include "condensate.h"

//////////////////////////////////////////////////////////
//             random vector functions
//////////////////////////////////////////////////////////
/**
 * Allocate the memory of a Random vector.
 */
void AllocateRandomVec(struct RandomVec * vec)
{
    vec->vec = malloc(sizeof(Fermion));
    AllocateFermion(vec->vec);
}

/**
 * Destroy the Random vector (and in particular its fermion field).
 */
void DestroyRandomVec(struct RandomVec * vec)
{
    DestroyFermion(vec->vec);
    free(vec->vec);
    vec->vec = NULL;
}

/**
 * Set the full random vector to zero.
 */
void SetFullRandomVecToZero(struct RandomVec * vec)
{
    for(int x=0; x<VOLUME; x++) vec->vec->fermion[x] = 0.0;
}

/**
 * Initialize half Random vector with Z2 noise.
 * This means that each vector element of Random Vec is 1.0 or -1.0 with probability 0.5.
 * The HALFVOLUME complementary elements are not initialized.
 */
void Z2NoiseHalfRandomVec(struct RandomVec * vec, int start)
{
    for(int x=start; x<start+HALFVOLUME; x++)
    {
        // with probability 0.5 set 1.0
        if(Casuale() <= 0.5) vec->vec->fermion[x] = 1.0;
        else vec->vec->fermion[x] = -1.0;
    }
}

/**
 * Initialize the full Random vector with Z2 noise.
 * This means that each vector element of Random Vec is 1.0 or -1.0 with probability 0.5.
 */
void Z2NoiseFullRandomVec(struct RandomVec * vec)
{
    for(int x=0; x<VOLUME; x++)
    {
        // with probability 0.5 set 1.0
        if(Casuale() <= 0.5) vec->vec->fermion[x] = 1.0;
        else vec->vec->fermion[x] = -1.0;
    }
}

/**
 * Set all components of a Vec (it is no more random) to 1.0 / VOLUME.
 */
void FullVecOneOverVolume(struct RandomVec * vec)
{
    for(int x=0; x<VOLUME; x++) vec->vec->fermion[x] = 1.0 / VOLUME;
}

/**
 * Set all components of a vec (it is no more random) to exp(i (\pi / L) x_3).
 */
void FullVecAllPhaseOverVolume(struct RandomVec * vec, const struct Lattice const * lattice)
{
    for(int x=0; x<VOLUME; x++)
    vec->vec->fermion[lattice->fermionsite[x]] = cexp(I_UNIT * (P_MIN_APBC) * ( lattice->lattice_site[COOR(x, 0)] + lattice->lattice_site[COOR(x, 1)] + lattice->lattice_site[COOR(x, 2)])) / VOLUME;
}

/**
 * Set all components of a vec (it is no more random) to the proper phase exp(i p x) where
 * p is the minimum momentum on the lattice and x is the position along the direction d passed
 * to the function.
 */
void FullVecPhasesForChipOverVolume(struct RandomVec * vec, int d, const struct Lattice const * lat)
{
    // for all lattice sites in lessicographic coordinates
    for(int x=0; x<VOLUME; x++)
    vec->vec->fermion[lat->fermionsite[x]] = cexp(I_UNIT * (P_MIN_APBC) * lat->lattice_site[COOR(x, d)]) / VOLUME;
}

/**
 * Initialize the full random vector with complex U(1) phases.
 * This means that each element in the vector is of the form R_x = exp(i * theta_x).
 */
void PhaseFullNoise(struct RandomVec * vec)
{
    for(int x=0; x<VOLUME; x++)
    {
        // extract a random phase in [-PI, PI) and save into the fermion as exp(I * phase)
        double phase = 0.0;
        SymmetricCasuale(PI, &phase);
        vec->vec->fermion[x] = cexp(I_UNIT * phase);
    }
}

/**
 * Compute the trace of the inverse Dirac matrix and its cumulants.
 * To evaluate Tr( D^{-1} ), first note that Tr( D^{-1} ) = Tr( D^\dagger * ((D^\dagger)(D))^{-1} ).
 * The function uses the conjugate gradient method used to invert positive definite matrix.
 * --------------------------------------------------------------------------------------
 * Return : Tr( (D^\dagger) * (DD^\dagger)^{-1} ).
 * --------------------------------------------------------------------------------------
 */
void ComputeCondensateDirac(const struct Lattice const * lat, double * condensate_tr1, double * condensate_tr1_abs, double * condensate_tr2, double * condensate_tr4, double * condensate_mu2, double * condensate_mu4, double * condensate_chi, double * condensate_chi_abs, double * condensate_chip, double * condensate12_nosigma, double * condensate12_withsigma, struct RandomVec * vec, struct Fermion * target_inverted, const struct Info const * info, const struct Fermioncoor const * fermioncoor, double inversion_precision, struct Fermionkit * kit, int f)
{
    // initialize condensate variables
    (*condensate_tr1) = 0.0; (*condensate_tr1_abs) = 0.0; (*condensate_tr2) = 0.0; (*condensate_tr4) = 0.0; (*condensate_mu2) = 0.0; (*condensate_mu4) = 0.0; (*condensate_chi) = 0.0; (*condensate_chi_abs) = 0.0; (*condensate_chip) = 0.0; (*condensate12_nosigma) = 0.0; (*condensate12_withsigma) = 0.0;

    // as we have to compute the fourth cumulant of the distribution p->4, create a 4-length temp variable
    // the temp variables are initialized in the next cycle
    // we want to compute until the 4^th power of the trace of the condensate
    double temp[4*CONDENSATE_PRECISION];
    // initialize temp array
    for(int p=0; p<4*CONDENSATE_PRECISION; p++) temp[p] = 0.0;

    // calculate 4 * CONDENSATE_PRECISION arrays 
    for(int p=0; p<4*CONDENSATE_PRECISION; p++)
    {
        // Initialize the random vector to compute the trace of the inverse of the Dirac matrix
        #ifdef Z2_NOISE
        Z2NoiseFullRandomVec(vec);
        #endif
        #ifdef PHASE_NOISE
        PhaseFullNoise(vec);
        #endif

        // invert all the elements of "vec"
        ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, vec->vec, target_inverted, info, fermioncoor, inversion_precision, kit, f);
        
        // multiply the inverted vec (which is stored into target_inverted) for the adjoint Dirac matrix
        AdjDiracFermion(lat, target_inverted, &(kit->Ddag_inverted_phi_obs), info, fermioncoor);

        // compute one among the 4 * CONDENSATE_PRECISION elements of temp[p] that is given by conj(v_x) D^\dagger(DD^\dagger)^{-1}_{xy} v_y
        for(int x=0; x<VOLUME; x++) temp[p] += creal(conj(vec->vec->fermion[x]) * kit->Ddag_inverted_phi_obs.fermion[x]) / VOLUME;
    }
    
    // we can extract 4 different values for the expectation value of Tr M^{-1}
    double tr0 = Mean(&(temp[0*CONDENSATE_PRECISION]), CONDENSATE_PRECISION);
    double tr1 = Mean(&(temp[1*CONDENSATE_PRECISION]), CONDENSATE_PRECISION);
    double tr2 = Mean(&(temp[2*CONDENSATE_PRECISION]), CONDENSATE_PRECISION);
    double tr3 = Mean(&(temp[3*CONDENSATE_PRECISION]), CONDENSATE_PRECISION);
    // save the value of the power of the trace of the condensate in the proper observables
    // the trace of the inverse Dirac matrix is given by the average over CONDENSATE_PRECISION measures
    (*condensate_tr1) = (tr0 + tr1 + tr2 + tr3)/4.0;
    (*condensate_tr1_abs) = fabs((*condensate_tr1));
    (*condensate_tr2) = ((tr0 * tr1) + (tr0 * tr2) + (tr0 * tr3) + (tr1 * tr2) + (tr1 * tr3) + (tr2 * tr3))/6.0;
    (*condensate_tr4) = (tr0 * tr1 * tr2 * tr3);


    // compute momentum of D^{-1}_{xy}
    // in particular, use the fact that mu_2 = V^{-2}\sum_{xy}D^{-1}_{xy}
    //                                  mu_22 = (V^{-2}\sum_{xy})^2
    //                                  chi = V^{-1}\sum_{xy}D^{-1}_{xy}
    //                                  chip = V^{-1}\sum_{xy}exp(i p_m (x-y))D^{-1}_{xy}
    
    // ----------------------- COMPUTE CHI -------------------------------------
    FullVecOneOverVolume(vec);
    // invert "vec"
    ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, vec->vec, target_inverted, info, fermioncoor, inversion_precision, kit, f);
        
    // multiply the inverted vec (which is stored into target_inverted) for the adjoint Dirac matrix
    AdjDiracFermion(lat, target_inverted, &(kit->Ddag_inverted_phi_obs), info, fermioncoor);

    // compute mu_2
    for(int x=0; x<VOLUME; x++) {(*condensate_mu2) += creal(conj(vec->vec->fermion[x]) * kit->Ddag_inverted_phi_obs.fermion[x]);}
    //(*condensate_mu2) += creal(kit->Ddag_inverted_phi_obs.fermion[0]);

    // mu_2 is computed, calculate chi and mu2^2 ~ \mu_4
    (*condensate_mu4) = pow((*condensate_mu2), 2);
    (*condensate_chi) = VOLUME * (*condensate_mu2);
    (*condensate_chi_abs) = fabs((*condensate_chi));


    // ----------------------- COMPUTE CHI -------------------------------------
    FullVecPhasesForChipOverVolume(vec, 2, lat);
    // invert "vec"
    ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, vec->vec, target_inverted, info, fermioncoor, inversion_precision, kit, f);
        
    // multiply the inverted vec (which is stored into target_inverted) for the adjoint Dirac matrix
    AdjDiracFermion(lat, target_inverted, &(kit->Ddag_inverted_phi_obs), info, fermioncoor);

    // compute mu_2
    double temp_chip = 0.0;
    for(int x=0; x<VOLUME; x++) {temp_chip += cimag(conj(vec->vec->fermion[x]) * kit->Ddag_inverted_phi_obs.fermion[x]);}
    //(*condensate_mu2) += creal(kit->Ddag_inverted_phi_obs.fermion[0]);    
    // save the value of mixed chi-sigma correlator
    (*condensate_chip) = (temp_chip * VOLUME);

    //---------------------------------- compute susceptibility of \overline{\chi}^1\chi^2 ----------------------------------
    // create the z_2 random vector
    Z2NoiseFullRandomVec(vec);
    // copy of the z2 random vectors
    double copy_z2[VOLUME];
    for(int x=0; x<VOLUME; x++) copy_z2[x] = creal(vec->vec->fermion[x]);

    // invert the matrix so that sigma_x \to M^{-2}_{yx} sigma_x
    ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, vec->vec, target_inverted, info, fermioncoor, inversion_precision, kit, f);
    AdjDiracFermion(lat, target_inverted, &(kit->Ddag_inverted_phi_obs), info, fermioncoor);
    ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, &(kit->Ddag_inverted_phi_obs), target_inverted, info, fermioncoor, inversion_precision, kit, f);
    AdjDiracFermion(lat, target_inverted, &(kit->Ddag_inverted_phi_obs), info, fermioncoor);

    for(int x=0; x<VOLUME; x++) (*condensate12_nosigma) += kit->Ddag_inverted_phi_obs.fermion[x] * copy_z2[x] / VOLUME;

    //---------------------------------- compute susceptibility of \overline{\chi}^1\chi^2\sigma ----------------------------------
    // create the z_2 random vector
    Z2NoiseFullRandomVec(vec);
    // copy of the z2 random vectors
    for(int x=0; x<VOLUME; x++) copy_z2[x] = creal(vec->vec->fermion[x]);

    // multiply by sigmamass
    for(int x=0; x<VOLUME; x++) vec->vec->fermion[x] = vec->vec->fermion[x] * lat->sigmamass[x]; 
    // invert the matrix so that sigma_x \to M^{-1}_{yx} sigma_x
    ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, vec->vec, target_inverted, info, fermioncoor, inversion_precision, kit, f);
    AdjDiracFermion(lat, target_inverted, &(kit->Ddag_inverted_phi_obs), info, fermioncoor);
    for(int y=0; y<VOLUME; y++) kit->Ddag_inverted_phi_obs.fermion[y] = kit->Ddag_inverted_phi_obs.fermion[y] * lat->sigmamass[y];

    ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, &(kit->Ddag_inverted_phi_obs), target_inverted, info, fermioncoor, inversion_precision, kit, f);
    AdjDiracFermion(lat, target_inverted, &(kit->Ddag_inverted_phi_obs), info, fermioncoor);

    for(int x=0; x<VOLUME; x++) (*condensate12_withsigma) += kit->Ddag_inverted_phi_obs.fermion[x] * copy_z2[x] / VOLUME;

    /*
    // ------------------------ COMPUTE CHIP ----------------------------------
    // initialize the vec as a vector of phases (along the temporal direction)
    FullVecPhasesForChipOverVolume(vec, (DIM-1), lat);
    // invert "vec"
    ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(lat, vec->vec, target_inverted, info, fermioncoor, inversion_precision, kit, f);
        
    // multiply the inverted vec (which is stored into target_inverted) for the adjoint Dirac matrix
    AdjDiracFermion(lat, target_inverted, &(kit->Ddag_inverted_phi_obs), info, fermioncoor);
    // compute mu_2
    for(int x=0; x<VOLUME; x++) {(*condensate_chi_abs) += creal(conj(vec->vec->fermion[x]) * kit->Ddag_inverted_phi_obs.fermion[x]) / VOLUME;}
    */
}
