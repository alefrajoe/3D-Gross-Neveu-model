#include "conjugategradient.h"

///////////////////////////////////////////////////////////////////
//                   conjugate gradient
///////////////////////////////////////////////////////////////////

/**
 * Compute in a single step the inverted Phi fermion vector.
 * Note that differently from the half analogous function, this conjugate gradient function does not require
 * the start integer number.
 * If the inversion has been carried out correctly return TRUE, else FALSE.
 */
bool ConjugateGradientFastDiracAdjFastDiracToTheMinusOne(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target_inverted, const struct Info const * info, const struct Fermioncoor const * fermioncoor, double inversion_precision, struct Fermionkit * kit, int f)
{
    // set the correct flavor for inversion
    #ifdef START_X0_NULL
    // ************************** check if passed fermion is already inverted ********************
    CopyFermion(fermion, &(kit->r_k.flav[f])); 
    // compute temp Ax_0
    AdjDiracFermion(lat, target_inverted, &(kit->temp1), info, fermioncoor);
    SubtractDiracFermion(lat, &(kit->temp1), &(kit->r_k.flav[f]), info, fermioncoor);
    if(NormFermion(&(kit->r_k.flav[f])) < inversion_precision) { 
                                                                 #ifdef DEBUG_STARTING_X0
                                                                 printf("Iteration : NULL \n");
                                                                 #endif
                                                                 // if zero moves are done, do not consider them for CG iteration count!!!
                                                                 return true;}
    //********************************************************************************************
    // start the inverted vector to X_0 = 0
    // As X_0 is equal to 0, the starting residues is given by  r_0 = F
    SetFermionToZero(target_inverted);
    CopyFermion(fermion, &(kit->r_k.flav[f]));
    CopyFermion(fermion, &(kit->p_k.flav[f]));
    #endif 

    #ifdef START_X0_PREVIOUS
    // starting r_k is given by r_k = b - A x_0
    CopyFermion(fermion, &(kit->r_k.flav[f])); 
    // compute temp Ax_0
    AdjDiracFermion(lat, target_inverted, &(kit->temp1), info, fermioncoor);
    SubtractDiracFermion(lat, &(kit->temp1), &(kit->r_k.flav[f]), info, fermioncoor);
    if(NormFermion(&(kit->r_k.flav[f])) < inversion_precision) { 
                                                                 #ifdef DEBUG_STARTING_X0
                                                                 printf("Iteration : -\n");
                                                                 #endif
                                                                 // if zero moves are done, do not consider them for CG iteration count!!!
                                                                 return true;}
    CopyFermion(&(kit->r_k.flav[f]), &(kit->p_k.flav[f]));
    #endif

    #ifdef DEBUG_STARTING_X0
    printf("Starting distance is %.16f\n", BMinusAX0(lat, fermion, target_inverted, info, fermioncoor, kit));
    #endif

    // for all families of fermion flavors
    // start the iteration to look for the inverted vector
    // the maximum number of iterations that we will perform is equal to VOLUME
    for(int iteration=0; iteration<MAX_CONJUGATE_GRADIENT_ITERATION; iteration++)
    {        
        // compute alpha_k
        double alpha_k = 0.0, norm_r_k = 0.0, norm_p_k = 0.0;

        // Numerator : norm of half vector r_k -> r_k^\dagger r_k
        for(int x=0; x<VOLUME; x++)
        norm_r_k += pow(cabs(kit->r_k.flav[f].fermion[x]), 2.0);

        // Denominator : norm of half ( D ) r_k -> r_k^\dagger ( D )( D^\dagger ) r_k
        AdjDiracFermion(lat, &(kit->p_k.flav[f]), &(kit->temp1), info, fermioncoor);
        for(int x=0; x<VOLUME; x++) norm_p_k += pow(cabs(kit->temp1.fermion[x]), 2.0);
        //norm_p_k = creal(AdjFermionAdjDiracDiracHalfFermion(&(kit->p_k), lat, &(kit->p_k), info, fermioncoor, start));
        
        // new \alpha_k is given by the ratio between the numerator and the denominator
        alpha_k = norm_r_k / norm_p_k;

        // the created trial vector X_{k+1} is given by X_{k+1} = X_{k} + \alpha_k p_k
        AddFermionWithCoefficientToTarget(&(kit->p_k.flav[f]), alpha_k, target_inverted);

        
        // continue the iteration for the residues.
        // at this point the residues is r_{k+1} = r_k - \alpha_k A p_k
        // ...
        // first compute A p_K
        DiracFermion(lat, &(kit->temp1), &(kit->temp2), info, fermioncoor);
        AddFermionWithCoefficientToTarget(&(kit->temp2), - alpha_k, &(kit->r_k.flav[f]));

        // if the residue is sufficiently small, exit the function
        double intermediate_norm = NormFermion(&(kit->r_k.flav[f]));
        // if the residue is sufficiently small, exit the loop
        if(intermediate_norm < inversion_precision){

                                                    // update CG info
                                                    kit->iterations += (iteration+1);
                                                    kit->CG_call++;
                                                    return true;}

        // compute \beta_{\kappa} 
        double beta_k = 0.0, norm_r_kplus1 = 0.0;
        for(int x=0; x<VOLUME; x++) norm_r_kplus1 += pow(cabs(kit->r_k.flav[f].fermion[x]), 2.0);
        beta_k = norm_r_kplus1 / norm_r_k;

        // compute p_{k+1} = r_{k+1} + \beta_{k} p_{k}
        for(int x=0; x<VOLUME; x++) kit->p_k.flav[f].fermion[x] = kit->r_k.flav[f].fermion[x] + beta_k * kit->p_k.flav[f].fermion[x];
    }

    // compute the sum of the norms of all vectors
    return false;
}

/**
 * Compute in a single step the inverted Phi fermion vector.
 * Note that differently from the half analogous function, this conjugate gradient function does not require
 * the start integer number.
 * If the inversion has been carried out correctly return TRUE, else FALSE.
 */
bool StartX0ImposedConjugateGradientFastDiracAdjFastDiracToTheMinusOne(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target_inverted, const struct Info const * info, const struct Fermioncoor const * fermioncoor, double inversion_precision, struct Fermionkit * kit, int f)
{
    // set the correct flavor for inversion
    // start the inverted vector to X_0 = 0
    // As X_0 is equal to 0, the starting residues is given by  r_0 = F
    SetFermionToZero(target_inverted);
    CopyFermion(fermion, &(kit->r_k.flav[f]));
    // if the residues is sufficiently small return
    if(NormFermion(&(kit->r_k.flav[f])) < inversion_precision) {
                                                                 #ifdef DEBUG_STARTING_X0
                                                                 printf("Iteration : -\n");
                                                                 #endif
                                                                 // if zero moves are done, do not consider them for CG iteration count!!!
                                                                 return true;}
    CopyFermion(fermion, &(kit->p_k.flav[f]));

    // for all families of fermion flavors
    // start the iteration to look for the inverted vector
    // the maximum number of iterations that we will perform is equal to VOLUME
    for(int iteration=0; iteration<MAX_CONJUGATE_GRADIENT_ITERATION; iteration++)
    {        
        // compute alpha_k
        double alpha_k = 0.0, norm_r_k = 0.0, norm_p_k = 0.0;

        // Numerator : norm of half vector r_k -> r_k^\dagger r_k
        for(int x=0; x<VOLUME; x++)
        norm_r_k += pow(cabs(kit->r_k.flav[f].fermion[x]), 2.0);

        // Denominator : norm of half ( D ) r_k -> r_k^\dagger ( D )( D^\dagger ) r_k
        AdjDiracFermion(lat, &(kit->p_k.flav[f]), &(kit->temp1), info, fermioncoor);
        for(int x=0; x<VOLUME; x++) norm_p_k += pow(cabs(kit->temp1.fermion[x]), 2.0);
        //norm_p_k = creal(AdjFermionAdjDiracDiracHalfFermion(&(kit->p_k), lat, &(kit->p_k), info, fermioncoor, start));
        
        // new \alpha_k is given by the ratio between the numerator and the denominator
        alpha_k = norm_r_k / norm_p_k;

        // the created trial vector X_{k+1} is given by X_{k+1} = X_{k} + \alpha_k p_k
        AddFermionWithCoefficientToTarget(&(kit->p_k.flav[f]), alpha_k, target_inverted);

        
        // continue the iteration for the residues.
        // at this point the residues is r_{k+1} = r_k - \alpha_k A p_k
        // ...
        // first compute A p_K
        DiracFermion(lat, &(kit->temp1), &(kit->temp2), info, fermioncoor);
        AddFermionWithCoefficientToTarget(&(kit->temp2), - alpha_k, &(kit->r_k.flav[f]));

        // if the residue is sufficiently small, exit the function
        double intermediate_norm = NormFermion(&(kit->r_k.flav[f]));
        // if the residue is sufficiently small, exit the loop
        if(intermediate_norm < inversion_precision){
                                                    // update CG info
                                                    kit->iterations += (iteration+1);
                                                    kit->CG_call++;
                                                    return true;}

        // compute \beta_{\kappa} 
        double beta_k = 0.0, norm_r_kplus1 = 0.0;
        for(int x=0; x<VOLUME; x++) norm_r_kplus1 += pow(cabs(kit->r_k.flav[f].fermion[x]), 2.0);
        beta_k = norm_r_kplus1 / norm_r_k;

        // compute p_{k+1} = r_{k+1} + \beta_{k} p_{k}
        for(int x=0; x<VOLUME; x++) kit->p_k.flav[f].fermion[x] = kit->r_k.flav[f].fermion[x] + beta_k * kit->p_k.flav[f].fermion[x];
    }

    // compute the sum of the norms of all vectors
    return false;
}

/**
 * Compute in a single step the inverted Phi fermion vector.
 * This function always uses the dirac and the adjdirac matrix with PBC intended (no additional phases stored into eta).
 * Note that differently from the half analogous function, this conjugate gradient function does not require
 * the start integer number.
 * If the inversion has been carried out correctly return TRUE, else FALSE.
 */
bool PBCConjugateGradientFastDiracAdjFastDiracToTheMinusOne(const struct Lattice const * lat, const struct Fermion const * fermion, struct Fermion * target_inverted, const struct Info const * info, const struct Fermioncoor const * fermioncoor, double inversion_precision, struct Fermionkit * kit, int f)
{
    // set the correct flavor for inversion
    // start the inverted vector to X_0 = 0
    // As X_0 is equal to 0, the starting residues is given by  r_0 = F
    SetFermionToZero(target_inverted);
    CopyFermion(fermion, &(kit->r_k.flav[f]));
    CopyFermion(fermion, &(kit->p_k.flav[f]));

    #ifdef DEBUG_STARTING_X0
    printf("Starting distance is %.16f\n", BMinusAX0(lat, fermion, target_inverted, info, fermioncoor, kit));
    #endif

    // for all families of fermion flavors
    // start the iteration to look for the inverted vector
    // the maximum number of iterations that we will perform is equal to VOLUME
    for(int iteration=0; iteration<MAX_CONJUGATE_GRADIENT_ITERATION; iteration++)
    {        
        // compute alpha_k
        double alpha_k = 0.0, norm_r_k = 0.0, norm_p_k = 0.0;

        // Numerator : norm of half vector r_k -> r_k^\dagger r_k
        for(int x=0; x<VOLUME; x++)
        norm_r_k += pow(cabs(kit->r_k.flav[f].fermion[x]), 2.0);

        // Denominator : norm of half ( D ) r_k -> r_k^\dagger ( D )( D^\dagger ) r_k
        PBCAdjDiracFermion(lat, &(kit->p_k.flav[f]), &(kit->temp1), info, fermioncoor);
        for(int x=0; x<VOLUME; x++) norm_p_k += pow(cabs(kit->temp1.fermion[x]), 2.0);
        //norm_p_k = creal(AdjFermionAdjDiracDiracHalfFermion(&(kit->p_k), lat, &(kit->p_k), info, fermioncoor, start));
        
        // new \alpha_k is given by the ratio between the numerator and the denominator
        alpha_k = norm_r_k / norm_p_k;

        // the created trial vector X_{k+1} is given by X_{k+1} = X_{k} + \alpha_k p_k
        AddFermionWithCoefficientToTarget(&(kit->p_k.flav[f]), alpha_k, target_inverted);

        
        // continue the iteration for the residues.
        // at this point the residues is r_{k+1} = r_k - \alpha_k A p_k
        // ...
        // first compute A p_K
        PBCDiracFermion(lat, &(kit->temp1), &(kit->temp2), info, fermioncoor);
        AddFermionWithCoefficientToTarget(&(kit->temp2), - alpha_k, &(kit->r_k.flav[f]));

        // if the residue is sufficiently small, exit the function
        double intermediate_norm = NormFermion(&(kit->r_k.flav[f]));
        // if the residue is sufficiently small, exit the loop
        if(intermediate_norm < inversion_precision){

                                                    // update CG info
                                                    kit->iterations += (iteration+1);
                                                    kit->CG_call++;
                                                    return true;}

        // compute \beta_{\kappa} 
        double beta_k = 0.0, norm_r_kplus1 = 0.0;
        for(int x=0; x<VOLUME; x++) norm_r_kplus1 += pow(cabs(kit->r_k.flav[f].fermion[x]), 2.0);
        beta_k = norm_r_kplus1 / norm_r_k;

        // compute p_{k+1} = r_{k+1} + \beta_{k} p_{k}
        for(int x=0; x<VOLUME; x++) kit->p_k.flav[f].fermion[x] = kit->r_k.flav[f].fermion[x] + beta_k * kit->p_k.flav[f].fermion[x];
    }

    // compute the sum of the norms of all vectors
    return false;
}

/**
 * The function returns back the norm of the || b - A x_0 ||, where b and x_0 are two fermion vectors
 * and A is (D D^\dagger).
 * Note that when the sigma or pi field are defined, it does not hold anymore that (D D^\dagger) = (D^\dagger D)
 */
double BMinusAX0(const struct Lattice const * lat, const struct Fermion const * B, const struct Fermion const * x0, const struct Info const * info, const struct Fermioncoor * fcor, struct Fermionkit * kit)
{
    double distance = 0.0;

    // copy B into bminusAx0
    CopyFermion(B, &(kit->bminusax0));

    // calculate A x_0
    AdjDiracFermion(lat, x0, &(kit->tempbminusx0), info, fcor);
    SubtractDiracFermion(lat, &(kit->tempbminusx0), &(kit->bminusax0), info, fcor);

    // return the norm of the fermion
    return NormFermion(&(kit->bminusax0));
}