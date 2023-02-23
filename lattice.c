#include "lattice.h"

/////////////////////////////////////////////////////////////
//               Sites functions
/////////////////////////////////////////////////////////////
/**
 * Construct the site passed to the lattice passed to the function. 
 * The site number, passed in Lessicographic COORDINATES, is saved into the site structure. 
 * Moreover, the links associated with the site are initialized according to the "START_LINK.." macro,
 * as well as the pp[DIM] and mm[DIM] arrays, containing the Lessicographic coordinates of the next and previous sites along 
 * the 0 <= \mu < DIM directions.
 * The fermion coordinates are related to the standard Lessicographic ones: the first HALFVOLUME sites are the even sites,
 *                                                                         the remnant ones are the odd ones.
 */
void StartSite(struct Lattice * lattice, int site_number)
{
    // site number (Lessicographic coordinates) and site in lattice coordinates
    lattice->site[site_number] = site_number;
    ConvertLessicographicIntoCoordinate(site_number, &(lattice->lattice_site[COOR(site_number, 0)]));
    // start all links
    for(int d=0; d<DIM; d++) StartLink(&(lattice->dlink[COOR(site_number, d)]));
    // save pp and mm arrays
    InitializePPArray(lattice->site[site_number], &(lattice->pp[COOR(site_number, 0)]));
    InitializeMMArray(lattice->site[site_number], &(lattice->mm[COOR(site_number, 0)]));

    ////////////// fermion sector //////////////////
    // temp array containing the site in lattice coordinate
    int temp[DIM];
    // Convert Lessicographic into fermion Lessicographic
    ConvertLessicographicIntoFermionLessicographic(site_number, &(lattice->fermionsite[site_number]));
    // save fermionic pp and mm arrays
    InitializeFermionPPArray(lattice->fermionsite[site_number], &(lattice->fermionpp[COOR(site_number, 0)]));
    InitializeFermionMMArray(lattice->fermionsite[site_number], &(lattice->fermionmm[COOR(site_number, 0)]));

    //InitializeFastFermionPPArray(site_number, &(lattice->fermionpp[COOR(site_number, 0)]), lattice->fermionsite);
    // start the eta_pp phase
    StartEtapp(lattice->site[site_number], &(lattice->eta_pp[COOR(site_number, 0)]));
    // start the eta_pp phase
    StartEtamm(lattice->site[site_number], &(lattice->eta_mm[COOR(site_number, 0)]));
    // start the eta_pbc phase
    StartEtapbc(lattice->site[site_number], &(lattice->eta_pbc[COOR(site_number, 0)]));
    // start the epsilon phases at the proper FERMION_SITE
    StartEpsilonPi(lattice->site[site_number], lattice->fermionsite[site_number], lattice->epsilon_pi);


    #ifdef DEBUG_COORDINATE
    // for all VOLUME Lessicographic coordinates
    for(int x=0; x<VOLUME; x++)
    {
        // consider x as an initial Lessicographic coordinate
        int stereo, fermion;
        ConvertLessicographicIntoFermionLessicographic(x, &fermion);
        ConvertFermionLessicographicIntoLessicographic(fermion, &stereo);
        if(stereo != x) printf("Error in bijection in stereo -> fermion -> stereo !!!\n");

        // consider x as an initial fermionic coordinate
        ConvertFermionLessicographicIntoLessicographic(x, &stereo);
        ConvertLessicographicIntoFermionLessicographic(stereo, &fermion);
        if(fermion != x) printf("Error in bijection in fermion -> stereo -> fermion !!! %d -> %d -> %d\n", x, stereo, fermion);
    }
    #endif

}

/**
 * Print a site to standard output.
 * The site is passed in Lessicographic coordinates.
 * In order, the following info are printed:
 * - site number
 * - pp[DIM] array
 * - mm[DIM] array
 * - fermionsite
 * - fermionpp[DIM] array
 * - fermionmm[DIM] array
 * - eta_\mu[DIM] array
 * - eta_mm_\mu[DIM] array
 * - link[DIM] array
 */
void PrintSite(const struct Lattice const * lattice, int site)
{
    printf("============= start ================\n");
    // site number
    printf("Site[%d]\n", lattice->site[site]);

    // site in lattice coordinates
    printf("Lattice pos is ");
    for(int d=0; d<DIM; d++) printf("[%d]", lattice->lattice_site[COOR(site, d)]);
    printf("\n");
    
    // pp array
    printf("array pp");
    for(int d=0; d<DIM; d++) printf("[%d]", lattice->pp[COOR(site, d)]);
    printf("\n");

    // mm array
    printf("array mm");
    for(int d=0; d<DIM; d++) printf("[%d]", lattice->mm[COOR(site, d)]);
    printf("\n");

    // fermion site
    printf("fermion site[%d]\n", lattice->fermionsite[site]);

    // fermion pp array
    printf("array fermion pp");
    for(int d=0; d<DIM; d++) printf("[%d]", lattice->fermionpp[COOR(site, d)]);
    printf("\n");

    // mm array
    printf("array fermion mm");
    for(int d=0; d<DIM; d++) printf("[%d]", lattice->fermionmm[COOR(site, d)]);
    printf("\n");

    // eta_mu (each element is + or -)
    printf("eta_pp ");
    for(int d=0; d<DIM; d++) printf("[%.6f]", lattice->eta_pp[COOR(site, d)]);
    printf("\n");
    // eta_mu (each element is + or -)
    printf("eta_mm ");
    for(int d=0; d<DIM; d++) printf("[%.6f]", lattice->eta_mm[COOR(site, d)]);
    printf("\n");
    // eta_mu (each element is + or -)
    printf("eta_pbc ");
    for(int d=0; d<DIM; d++) printf("[%.6f]", lattice->eta_pbc[COOR(site, d)]);
    printf("\n");


    // ifdef SIGMA_INTERACTION, print the sigma value and the sigma mass
    #ifdef SIGMA_INTERACTION
    printf("sigma %.16f\n", lattice->sigma[lattice->fermionsite[site]]);
    printf("sigmamass %.16f\n", lattice->sigmamass[lattice->fermionsite[site]]);
    #endif

    #ifdef PI_INTERACTION
    printf("pi %.16f\n", lattice->pi[lattice->fermionsite[site]]);
    printf("pimass %.16f\n", lattice->pi[lattice->fermionsite[site]]);
    printf("epsilon_pi[%.6f]\n", lattice->epsilon_pi[lattice->fermionsite[site]]);
    #endif

    #ifdef DEBUG_COORDINATE
    // fermion pp array in standard coordinate
    printf("array fermion pp in standard coordinate");
    for(int d=0; d<DIM; d++)
    {
        int temp;
        ConvertFermionLessicographicIntoLessicographic(site->fermionpp[d], &temp);
        printf("[%d]", temp);
    }
    printf("\n");
    // fermion mm array in standard coordinate
    printf("array fermion mm in standard coordinate");
    for(int d=0; d<DIM; d++)
    {
        int temp;
        ConvertFermionLessicographicIntoLessicographic(site->fermionmm[d], &temp);
        printf("[%d]", temp);
    }
    printf("\n");
    #endif

    // links
    for(int d=0; d<DIM; d++)
    {
        printf("Link %d: ", d);
        PrintLink(&(lattice->dlink[COOR(site, d)]));
    }
    printf("============== end =================\n");
}

///////////////////////////////////////////////////////////////////////////////
//                    Lattice functions
///////////////////////////////////////////////////////////////////////////////
/**
 * The function allocates the memory to the passed lattice.
 */
void AllocateLattice(struct Lattice * lattice)
{
    lattice->dlink = malloc(sizeof(Link) * VOLUME * DIM);
    lattice->site = malloc(sizeof(int) * VOLUME);
    lattice->lattice_site = malloc(sizeof(int) * VOLUME * DIM);
    lattice->pp = malloc(sizeof(int) * VOLUME * DIM);
    lattice->mm = malloc(sizeof(int) * VOLUME * DIM);
    lattice->fermionsite = malloc(sizeof(int) * VOLUME);
    lattice->fermionpp = malloc(sizeof(int) * VOLUME * DIM);
    lattice->fermionmm = malloc(sizeof(int) * VOLUME * DIM);
    lattice->eta_pp = malloc(sizeof(double) * VOLUME * DIM);
    lattice->eta_mm = malloc(sizeof(double) * VOLUME * DIM);
    lattice->eta_pbc = malloc(sizeof(double) * VOLUME * DIM);
    lattice->sigma = malloc(sizeof(double) * VOLUME);
    lattice->sigmamass = malloc(sizeof(double) * VOLUME);
    lattice->pi = malloc(sizeof(double) * VOLUME);
    lattice->epsilon_pi = malloc(sizeof(double) * VOLUME);
    lattice->pimass = malloc(sizeof(double) * VOLUME);

    // initialize trivially all to zero
    for(int x=0; x<VOLUME; x++)
    {
        lattice->site[x] = 0;
        lattice->fermionsite[x] = 0;
        lattice->sigma[x] = 0.0;
        lattice->sigmamass[x] = 0.0;
        lattice->pi[x] = 0.0;
        lattice->epsilon_pi[x] = 0.0;
        lattice->pimass[x] = 0.0;


        for(int d=0; d<DIM; d++)
        {
            lattice->lattice_site[COOR(x, d)] = 0;
            lattice->dlink[COOR(x, d)].link = 0.0;
            lattice->pp[COOR(x, d)] = 0;
            lattice->mm[COOR(x, d)] = 0;
            lattice->fermionpp[COOR(x, d)] = 0;
            lattice->fermionmm[COOR(x, d)] = 0;
            lattice->eta_pp[COOR(x, d)] = 0.0;
            lattice->eta_mm[COOR(x, d)] = 0.0;
            lattice->eta_pbc[COOR(x, d)] = 0.0;
        }
    }
}

/**
 * The function initializes the lattice. 
 * Each lattice sites is initialized with a site number (in Lessicographic coordinates) that ranges in [0, VOL).
 */
void StartLattice(struct Lattice * lattice)
{
    // first create the lattice sites and fermionsites
    for(int x=0; x<VOLUME; x++)
    {
        // site
        lattice->site[x] = x;
        // fermionsite
        int temp = 0;
        ConvertLessicographicIntoFermionLessicographic(x, &temp);
        lattice->fermionsite[x] = temp;
    }

    for(int x=0; x<VOLUME; x++)
    {
        // start all sites in lessicographic coordinates
        StartSite(lattice, x);
    }

    // initialize all distances to zero
    lattice->distance_link = 0.0;
    lattice->distance_sigma = 0.0;
    lattice->distance_pi = 0.0;

    lattice->done_link = 0;
    lattice->done_sigma = 0;
    lattice->done_pi = 0;
}

/**
 * The function frees the memory of a previously allocated lattice.
 */
void DestroyLattice(struct Lattice * lattice)
{
    free(lattice->dlink);
    free(lattice->site);
    free(lattice->lattice_site);
    free(lattice->pp);
    free(lattice->mm);
    free(lattice->fermionsite);
    free(lattice->fermionpp);
    free(lattice->fermionmm);
    free(lattice->eta_pp);
    free(lattice->eta_mm);
    free(lattice->eta_pbc);
    free(lattice->sigma);
    free(lattice->sigmamass);
    free(lattice->pi);
    free(lattice->epsilon_pi);
    free(lattice->pimass);

    lattice->dlink = NULL;
    lattice->site = NULL;
    lattice->lattice_site = NULL;
    lattice->pp = NULL;
    lattice->mm = NULL;
    lattice->fermionsite = NULL;
    lattice->fermionpp = NULL;
    lattice->fermionmm = NULL;
    lattice->eta_pp = NULL;
    lattice->eta_mm = NULL;
    lattice->eta_pbc = NULL;
    lattice->sigma = NULL;
    lattice->sigmamass = NULL;
    lattice->pi = NULL;
    lattice->epsilon_pi = NULL;
    lattice->pimass = NULL;
}

/**
 * Print each site of the lattice to stdout.
 */
void PrintLattice(const struct Lattice const * lattice)
{
    for(int x=0; x<VOLUME; x++)
    {
        printf("=================================\n");
        PrintSite(lattice, x);
        printf("=================================\n");
    }
}