#include "fermion.h"

//////////////////////////////////////////////////////////
//            Fermioncoor functions
//////////////////////////////////////////////////////////
/**
 * Allocate the memory for a Fermioncoor struct.
 */
void AllocateFermioncoor(struct Fermioncoor * fermioncoor)
{
    fermioncoor->site = malloc(sizeof(int) * VOLUME);
    fermioncoor->pp = malloc(sizeof(int) * VOLUME * DIM);
    fermioncoor->mm = malloc(sizeof(int) * VOLUME * DIM);
    fermioncoor->fermionsite = malloc(sizeof(int) * VOLUME);
    fermioncoor->fermionpp = malloc(sizeof(int) * VOLUME * DIM);
    fermioncoor->fermionmm = malloc(sizeof(int) * VOLUME * DIM);
    fermioncoor->eta_pp = malloc(sizeof(double) * VOLUME * DIM);
    fermioncoor->eta_mm = malloc(sizeof(double) * VOLUME * DIM);
    fermioncoor->eta_pbc = malloc(sizeof(double) * VOLUME * DIM);

    // initialize trivially to zero all elements
    for(int x=0; x<VOLUME; x++)
    {
        fermioncoor->site[x] = 0;
        fermioncoor->fermionsite[x] = 0;

        for(int d=0; d<DIM; d++)
        {
            fermioncoor->pp[COOR(x, d)] = 0;
            fermioncoor->mm[COOR(x, d)] = 0;
            fermioncoor->fermionpp[COOR(x, d)] = 0;
            fermioncoor->fermionpp[COOR(x, d)] = 0;
            fermioncoor->eta_pp[COOR(x, d)] = 0.0;
            fermioncoor->eta_mm[COOR(x, d)] = 0.0;
            fermioncoor->eta_pbc[COOR(x, d)] = 0.0;
        }
    }
}

/**
 * Start all the sites in ascending order of fermion coordinates.
 */ 
void StartFermioncoor(struct Fermioncoor * fermioncoor)
{
    for(int x=0; x<VOLUME; x++)
    {
        // temporary variable needed for the construction of the coordinates
        int stereo;
        ConvertFermionLessicographicIntoLessicographic(x, &stereo);

        // start fermion site 
        fermioncoor->fermionsite[x] = x;
        // start fermion pp sites
        InitializeFermionPPArray(x, &(fermioncoor->fermionpp[COOR(x, 0)]));
        // start fermion mm sites
        InitializeFermionMMArray(x, &(fermioncoor->fermionmm[COOR(x, 0)]));
        // start site in Lessicographic
        fermioncoor->site[x] = stereo;
        // start pp sites in Lessicographic
        InitializePPArray(stereo, &(fermioncoor->pp[COOR(x, 0)]));
        // start mm sites in Lessicographic
        InitializeMMArray(stereo, &(fermioncoor->mm[COOR(x, 0)]));
        // begin etapp
        StartEtapp(fermioncoor->site[x], &(fermioncoor->eta_pp[COOR(x, 0)]));
        // begin etamm
        StartEtamm(fermioncoor->site[x], &(fermioncoor->eta_mm[COOR(x, 0)]));
        // begin etapbc
        StartEtapbc(fermioncoor->site[x], &(fermioncoor->eta_pbc[COOR(x, 0)]));
    }
}

/**
 * Start a fermioncoordinate from a lattice that has been already constructed.
 */
void StartFermioncoorFromLattice(struct Fermioncoor * fcor, const struct Lattice const * lat)
{
    // for all sites, find the correspondent fermionsite
    for(int x=0; x<VOLUME; x++)
    {
        // fermion site
        int ferm_x = lat->fermionsite[x];

        // ------- Build fcor -----------
        fcor->site[ferm_x] = lat->site[x];
        fcor->fermionsite[ferm_x] = lat->fermionsite[x];

        for(int d=0; d<DIM; d++)
        {
            fcor->pp[COOR(ferm_x, d)] = lat->pp[COOR(x, d)];
            fcor->mm[COOR(ferm_x, d)] = lat->mm[COOR(x, d)];
            fcor->fermionpp[COOR(ferm_x, d)] = lat->fermionpp[COOR(x, d)];
            fcor->fermionmm[COOR(ferm_x, d)] = lat->fermionmm[COOR(x, d)];
            fcor->eta_pp[COOR(ferm_x, d)] = lat->eta_pp[COOR(x, d)];
            fcor->eta_mm[COOR(ferm_x, d)] = lat->eta_mm[COOR(x, d)];
            fcor->eta_pbc[COOR(ferm_x, d)] = lat->eta_pbc[COOR(x, d)];
        }
    }
}

/**
 * Print a fermioncoor site to stdoutput.
 * The site is passed in fermionic coordinates.
 */
void PrintFermioncoorSite(const struct Fermioncoor const * fermioncoor, int site)
{
    printf("============= start ================\n");
    // site number
    printf("Site[%d]\n", fermioncoor->site[site]);

    // site in lattice coordinates
    int latcor[DIM];
    ConvertLessicographicIntoCoordinate(fermioncoor->site[site], latcor);
    printf("Lattice pos is ");
    for(int d=0; d<DIM; d++) printf("[%d]", latcor[d]);
    printf("\n");
    
    // pp array
    printf("array pp");
    for(int d=0; d<DIM; d++) printf("[%d]", fermioncoor->pp[COOR(site, d)]);
    printf("\n");

    // mm array
    printf("array mm");
    for(int d=0; d<DIM; d++) printf("[%d]", fermioncoor->mm[COOR(site, d)]);
    printf("\n");

    // fermion site
    printf("fermion site[%d]\n", fermioncoor->fermionsite[site]);

    // fermion pp array
    printf("array fermion pp");
    for(int d=0; d<DIM; d++) printf("[%d]", fermioncoor->fermionpp[COOR(site, d)]);
    printf("\n");

    // mm array
    printf("array fermion mm");
    for(int d=0; d<DIM; d++) printf("[%d]", fermioncoor->fermionmm[COOR(site, d)]);
    printf("\n");

    // eta_mu (each element is + or -)
    printf("eta_pp ");
    for(int d=0; d<DIM; d++) printf("[%.6f]", fermioncoor->eta_pp[COOR(site, d)]);
    printf("\n");
    // eta_mu (each element is + or -)
    printf("eta_mm ");
    for(int d=0; d<DIM; d++) printf("[%.6f]", fermioncoor->eta_mm[COOR(site, d)]);
    printf("\n");
    // eta_mu pbc (each element is + or -)
    printf("eta_pbc ");
    for(int d=0; d<DIM; d++) printf("[%.6f]", fermioncoor->eta_pbc[COOR(site, d)]);
    printf("\n");
}

/**
 * Print to stdout a fermioncoor structure.
 * Note that the sites are ordered in ascending order of fermionsite.
 */
void PrintFermioncoor(const struct Fermioncoor const * fermioncoor)
{
    for(int x=0; x<VOLUME; x++)
    {
        printf("=================================\n");
        PrintFermioncoorSite(fermioncoor, x);
        printf("=================================\n");
    }
}

/**
 * Destroy a Fermioncoor coordinate container.
 */
void DestroyFermioncoor(struct Fermioncoor * fermioncoor)
{
    free(fermioncoor->site);
    free(fermioncoor->pp);
    free(fermioncoor->mm);
    free(fermioncoor->fermionsite);
    free(fermioncoor->fermionpp);
    free(fermioncoor->fermionmm);
    free(fermioncoor->eta_pp);
    free(fermioncoor->eta_mm);
    free(fermioncoor->eta_pbc);

    fermioncoor->site = NULL;
    fermioncoor->pp = NULL;
    fermioncoor->mm = NULL;
    fermioncoor->fermionsite = NULL;
    fermioncoor->fermionpp = NULL;
    fermioncoor->fermionmm = NULL;
    fermioncoor->eta_pp = NULL;
    fermioncoor->eta_mm = NULL;
    fermioncoor->eta_pbc = NULL;
}

/////////////////////////////////////////////////////////
//              fermion functions
/////////////////////////////////////////////////////////

/**
 * Allocate the memory of a fermion field (sizeof(complex) * VOLUME) numbers.
 */
void AllocateFermion(struct Fermion * fermion)
{
    fermion->fermion = malloc(sizeof(complex) * (VOLUME));

    // initialize all elements to zero
    for(int x=0; x<(VOLUME); x++) fermion->fermion[x] = 0.0;
}

/**
 * Deallocate the memory of a fermion field.
 */
void DestroyFermion(struct Fermion * fermion)
{
    free(fermion->fermion);
    fermion->fermion = NULL;
}

/**
 * Set the upper or the lower part of a fermion to zero.
 */
void SetHalfFermionToZero(struct Fermion * fermion, int start);

/**
 * Set the fermion array to zero.
 */
void SetFermionToZero(struct Fermion * fermion);

/**
 * Copy the upper or the lower part of a fermion.
 */
void CopyHalfFermion(const struct Fermion const * from, struct Fermion * target, int start);

/**
 * Copy the "from" fermion into the target fermion.
 */
void CopyFermion(const struct Fermion const * from, struct Fermion * target);

/**
 * Add to target the upper part of a fermion multiplied by a coefficient
 */
void AddFermionWithCoefficientToTarget(const struct Fermion const * fermion, complex coeff, struct Fermion * target);

/**
 * Print to stdout all the matrix elements of the fermion vector.
 */
void PrintFermion(const struct Fermion const * fermion)
{
    printf("========== start ============\n");
    // for all VOLUME elements in the vector (both upper and lower halves)
    for(int x=0; x<VOLUME; x++)
    {
        // temporary vector wherein the lattice coordinates will be stored
        int temp[DIM];
        ConvertFermionLessicographicIntoCoordinate(x, temp);

        printf("Fermion vector [%d] at ", x);
        for(int d=0; d<DIM; d++) printf("[%d]", temp[d]);
        printf(" is %.16f + i %.16f\n", creal(fermion->fermion[x]), cimag(fermion->fermion[x]));
    }
    printf("=========== end =============\n");
}

/**
 * Return the norm of the upper or the lower part of the passed fermion.
 */
double NormHalfFermion(const struct Fermion const * fermion, int start)
{
    double total = 0.0;

    for(int x=start; x<start+HALFVOLUME; x++)
    total += pow(cabs(fermion->fermion[x]), 2.0);

    return pow(total, 0.5);
}

/**
 * Print the norm of the full fermion vector.
 * It is the result of NORM = \sum_{i = 0} ^ {VOLUME - 1} psi^*_i psi_i.
 */
double NormFullFermion(const struct Fermion const * fermion)
{
    double total = 0.0;

    // sum over the full fermion vector
    for(int x=0; x<VOLUME; x++)
    total += pow(cabs(fermion->fermion[x]), 2.0);

    return pow(total, 0.5);
}

/**
 * Print the norm of the upper half fermion vector.
 * It is the result of NORM = \sum_{i = 0} ^ {VOLUME / 2 - 1} psi^*_i psi_i.
 */
double NormUpperHalfFullFermion(const struct Fermion const * fermion)
{
    double total = 0.0;

    // sum over the full fermion vector
    for(int x=0; x<(VOLUME / 2); x++)
    total += pow(cabs(fermion->fermion[x]), 2.0);

    return pow(total, 0.5);
}

/**
 * Compute the full norm of a fermion vector.
 */
double NormFermion(const struct Fermion const * fermion);

///////////////////////////////////////////////////////////////
//                   fermion for hybrid
///////////////////////////////////////////////////////////////
/**
 * Draw a number VOLUME / 2 of fermions (only the ones associated with the upper half of a fermion vector).
 * Both the real and imaginary part of a fermion field component is drawn from a Gaussian distribution
 * with 0 mean and 1.0 / 2 variance (this is because the \chi field appears in the distribution as \sim \exp(- \chi^\dagger \chi)).
 * The lower half of the vector is not initialized to zero.
 */
void DrawHalfFermionsFromGaussianDistribution(struct Fermion * fermion)
{
    // for the upper half of the fermion vector
    for(int r=0; r<HALFVOLUME; r++)
    {
        // pick two random doubles from a Gaussian distribution
        double num1, num2;
        Gauss2MeanVariance(&num1, &num2, 0.0, 1.0);

        // save into the fermion field
        fermion->fermion[r] = num1 + I_UNIT * num2;
    }
}

/**
 * Draw all the VOLUME elements of a fermion field from a Gaussian distribution.
 */
void DrawFullFermionFromGaussianDistribution(struct Fermion * fermion)
{
    // for the upper half of the fermion vector
    for(int r=0; r<VOLUME; r++)
    {
        // if Z_2 GN model is simulated, save real numbers into the fermion field
        fermion->fermion[r] = Gauss1();
    }
}

//////////////////////////////////////////////////////////////
//               many fermions function
//////////////////////////////////////////////////////////////

/**
 * Allocate the memory for HALFNFLAV fermions.
 */
void AllocateManyFermions(struct ManyFermions * fermions)
{
    fermions->flav = malloc(sizeof(Fermion) * HALFNFLAV);

    for(int f=0; f<HALFNFLAV; f++)
    AllocateFermion(&(fermions->flav[f]));
}

/**
 * Set all the fermion flavors to zero.
 */
void SetToZeroManyFermions(struct ManyFermions * fermions);

/**
 * Print to stdoutput all the fermions.
 */
void PrintFermions(const struct ManyFermions * fermions)
{
    // for all flavors in many flavors
    for(int f=0; f<HALFNFLAV; f++)
    {
        printf("========= start flav %d =============\n", f);
        PrintFermion(&(fermions->flav[f]));
        printf("========== end flav %d ==============\n", f);
    }
}

/**
 * Destruct a many fermion structure.
 */
void DestroyManyFermions(struct ManyFermions * fermions)
{
    for(int f=0; f<HALFNFLAV; f++)
    DestroyFermion(&(fermions->flav[f]));

    free(fermions->flav);
    fermions->flav = NULL;
}

////////////////////////////////////////////////////////////////////////
//                     fermionkit functions
////////////////////////////////////////////////////////////////////////

/**
 * Allocate the memory for all the instances of the fermionkit structure.
 */
void AllocateFermionkit(struct Fermionkit * fermionkit)
{
    AllocateManyFermions(&(fermionkit->p_k));
    AllocateManyFermions(&(fermionkit->r_k));
    AllocateManyFermions(&(fermionkit->chi));
    AllocateManyFermions(&(fermionkit->phi));
    AllocateManyFermions(&(fermionkit->inverted_phi));
    AllocateManyFermions(&(fermionkit->inverted_phi_obs));
    AllocateFermion(&(fermionkit->Ddag_inverted_phi_obs));
    AllocateFermion(&(fermionkit->temp1));
    AllocateFermion(&(fermionkit->temp2));
    AllocateFermion(&(fermionkit->tempbminusx0));
    AllocateFermion(&(fermionkit->bminusax0));

    // CG info
    fermionkit->iterations = 0.0;
    fermionkit->CG_call = 0;
}

/**
 * Destroy a fermiokit structure.
 */
void DestroyFermionkit(struct Fermionkit * fermionkit)
{
    DestroyManyFermions(&(fermionkit->p_k));
    DestroyManyFermions(&(fermionkit->r_k));
    DestroyManyFermions(&(fermionkit->chi));
    DestroyManyFermions(&(fermionkit->phi));
    DestroyManyFermions(&(fermionkit->inverted_phi));
    DestroyManyFermions(&(fermionkit->inverted_phi_obs));
    DestroyFermion(&(fermionkit->Ddag_inverted_phi_obs));
    DestroyFermion(&(fermionkit->temp1));
    DestroyFermion(&(fermionkit->temp2));
    DestroyFermion(&(fermionkit->tempbminusx0));
    DestroyFermion(&(fermionkit->bminusax0));
}