#include "simulation.h"

//////////////////////////////////////////////////////////////////////////////
//                      simulation function
//////////////////////////////////////////////////////////////////////////////

/**
 * Allocate all the elements required for a lattice simulation.
 */
void AllocateSimulation(struct Simulation * simulation)
{
    // create the memory
    simulation->lattice = malloc(sizeof(Lattice));
    simulation->fcor = malloc(sizeof(Fermioncoor));
    simulation->kit = malloc(sizeof(Fermionkit));
    simulation->force = malloc(sizeof(Force));
    simulation->momenta = malloc(sizeof(Momenta));
    simulation->momentasigma = malloc((sizeof(Momentasigma)));
    simulation->forcesigma = malloc((sizeof(Forcesigma)));
    simulation->momentapi = malloc(sizeof(Momentapi));
    simulation->forcepi = malloc(sizeof(Forcepi));

    // allocate stuff
    AllocateLattice(simulation->lattice);
    AllocateFermionkit(simulation->kit);
    AllocateFermioncoor(simulation->fcor);
    AllocateMomenta(simulation->momenta);
    AllocateForce(simulation->force);
    AllocateMomentaSigma(simulation->momentasigma);
    AllocateForceSigma(simulation->forcesigma);
    AllocateMomentaPi(simulation->momentapi);
    AllocateForcePi(simulation->forcepi);
}

/**
 * Start the several istances of the simulation structure.
 * The function also start the string of the output files and write the "title"
 * of the different output files.
 */
void StartSimulation(struct Simulation * simulation, double beta, double beta_scalar, double mass, double epsilon_link_HMC, double epsilon_sigma_HMC, double epsilon_pi_HMC, int step_link_HMC, int step_sigma_HMC, int step_pi_HMC, int sosia)
{
    // start all items required
    StartLattice(simulation->lattice);
    StartFermioncoorFromLattice(simulation->fcor, simulation->lattice);
    StartInfo(&(simulation->info), beta, beta_scalar, mass, epsilon_link_HMC, epsilon_sigma_HMC, epsilon_pi_HMC, step_link_HMC, step_sigma_HMC, step_pi_HMC, sosia);
    StartSigma(simulation->lattice, simulation->fcor);
    StartPi(simulation->lattice, simulation->fcor);

    // observables string and title
    FILE * opf;
    StringObservables(simulation->observables_file, &(simulation->info));
    StringAcceptance(simulation->acceptance_file, &(simulation->info));
    StringTerma(simulation->terma_file, simulation->terma_file_copy, &(simulation->info));
    OpenAndWriteObservablesTitle(opf, simulation->observables_file);
    OpenAndWriteAcceptanceTitle(opf, simulation->acceptance_file, &(simulation->info));
}

/**
 * Destroy all the istances of the simulation structure.
 */
void DestroySimulation(struct Simulation * simulation)
{
    // destroy all the variables
    DestroyLattice(simulation->lattice);
    DestroyFermionkit(simulation->kit);
    DestroyFermioncoor(simulation->fcor);
    DestroyForce(simulation->force);
    DestroyMomenta(simulation->momenta);
    DestroyMomentaSigma(simulation->momentasigma);
    DestroyForceSigma(simulation->forcesigma);
    DestroyMomentaPi(simulation->momentapi);
    DestroyForcePi(simulation->forcepi);

    // free pointers
    free(simulation->lattice);
    simulation->lattice = NULL;
    free(simulation->kit);
    simulation->kit = NULL;
    free(simulation->fcor);
    simulation->fcor = NULL;
    free(simulation->force);
    simulation->force = NULL;
    free(simulation->momenta);
    simulation->momenta = NULL;
    free(simulation->momentasigma);
    simulation->momentasigma = NULL;
    free(simulation->forcesigma);
    simulation->forcesigma = NULL;
    free(simulation->momentapi);
    simulation->momentapi = NULL;
    free(simulation->forcepi);
    simulation->forcepi = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////
//                       read and write lattice
//////////////////////////////////////////////////////////////////////////////////////

/**
 * The function is defined to test whether the lattice has been well read when READ_TERMA is defined.
 * =================================================================================================
 * Return : \sum_{DIM}\sum_{VOLUME} ( carg{U_mu(x)} + sigmamass[x] )
 * ===============================================================================================*/
double TestWriteReadLattice(const struct Simulation const * simulation)
{
    double total = 0.0;

    // link
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    {
        total += carg(simulation->lattice->dlink[COOR(x, d)].link);
        #ifdef SIGMA_INTERACTION
        total += simulation->lattice->sigmamass[x] / DIM;
        #endif
        #ifdef PI_INTERACTION
        total += simulation->lattice->pimass[x] / DIM;
        #endif
    }
    return total;
}

/**
 * The function write the lattice in a binary file.
 * In the file the following stuff are saved in the following order:
 * -
 * -
 */
void WriteBinaryLattice(const struct Simulation const * simulation)
{
    // if a terma copy file exist, remove the file
    struct stat buffer;
    if(stat(simulation->terma_file_copy, &buffer) == 0) remove(simulation->terma_file_copy);
    // if a written terma file exist, first rename this file into a copy
    // try to rename the file
    rename(simulation->terma_file, simulation->terma_file_copy);

    // open and write the terma file in binary
    FILE * opf=fopen(simulation->terma_file, "w");

    // return if there is an error
    if(opf == NULL)
    {
        fprintf(stderr, "It is impossible to open the file file in writing mode\n");
        exit(1);
    }

    // write the lattice in binary mode to a file
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    fwrite(&(simulation->lattice->dlink[COOR(x, d)].link), sizeof(complex), 1, opf);
    // ifdef SIGMA_INTERACTION, write down the sigma field
    #ifdef SIGMA_INTERACTION
    for(int x=0; x<VOLUME; x++)
    fwrite(&(simulation->lattice->sigma[x]), sizeof(double), 1, opf);
    #endif
    // ifdef PI_INTERACTION, write down the pi field
    #ifdef PI_INTERACTION
    for(int x=0; x<VOLUME; x++)
    fwrite(&(simulation->lattice->pi[x]), sizeof(double), 1, opf);
    #endif
    // test result
    double test = TestWriteReadLattice(simulation);
    fwrite(&test, sizeof(double), 1, opf);
    // write all the HMC parameters to be read i future simulations
    // the length of the HMC algorithm is not read
    fwrite(&(simulation->info.epsilon_link_HMC), sizeof(double), 1, opf);
    fwrite(&(simulation->info.story_link), sizeof(double), 1, opf);
    fwrite(&(simulation->info.hybrid_accepted_link), sizeof(int), 1, opf);
    fwrite(&(simulation->info.hybrid_done_link), sizeof(int), 1, opf);

    fwrite(&(simulation->info.story_sigma), sizeof(double), 1, opf);
    fwrite(&(simulation->info.hybrid_accepted_sigma), sizeof(int), 1, opf);
    fwrite(&(simulation->info.hybrid_done_sigma), sizeof(int), 1, opf);

    fwrite(&(simulation->info.epsilon_pi_HMC), sizeof(double), 1, opf);
    fwrite(&(simulation->info.story_pi), sizeof(double), 1, opf);
    fwrite(&(simulation->info.hybrid_accepted_pi), sizeof(int), 1, opf);
    fwrite(&(simulation->info.hybrid_done_pi), sizeof(int), 1, opf);
    // write also the measure number    
    fwrite(&(simulation->info.measure), sizeof(int), 1, opf);   

    // close the file
    fclose(opf);
}

/**
 * Read a simulation structure from a binary file.
 * The link of the lattice, the sigma and the pi field, as well as the HMC parameters and the measure number are read.
 */
void ReadBinaryLattice(struct Simulation * simulation, const char * terma_file, bool try_copy)
{
    FILE * opf=fopen(terma_file, "r");

    // if the file does not exist, simply return: this means this is this is the first run with the passed parameters
    // in this case, the simulation is started from the beginning, so that it requires a thermalization process
    if(opf == NULL)
    {
        // try to open the copy of the file
        opf=fopen(simulation->terma_file_copy, "r");
        // if there is an error also in the copy_terma_file return and start the simulation from 0
        if(opf == NULL)
        return;
    }

    int reading_errors = 0;

    // read the binary file and store it into simulation
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    if(fread(&(simulation->lattice->dlink[COOR(x, d)].link), sizeof(complex), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    // ifdef SIGMA_INTERACTION, read also the sigma fields
    #ifdef SIGMA_INTERACTION
    for(int x=0; x<VOLUME; x++)
    if(fread(&(simulation->lattice->sigma[x]), sizeof(double), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    // after we have read all the sigma fields, compute sigmamass
    ComputeFullSigmaMass(simulation->lattice, simulation->fcor, simulation->lattice->sigmamass);
    #endif

    // ifdef PI_INTERACTION, read also the pi field
    #ifdef PI_INTERACTION
    for(int x=0; x<VOLUME; x++)
    if(fread(&(simulation->lattice->pi[x]), sizeof(double), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    // after we have read all the sigma fields, compute sigmamass
    ComputeFullPiMass(simulation->lattice, simulation->fcor, simulation->lattice->sigmamass);
    #endif

    double test_old = 0.0, test_new = TestWriteReadLattice(simulation);
    if(fread(&(test_old), sizeof(double), 1, opf) != 1 || test_old != test_new)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }

    // READ HMC PARAMETERS --------------------------------------------------------------------------------------------------------
    if(fread(&(simulation->info.epsilon_link_HMC), sizeof(double), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.story_link), sizeof(double), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.hybrid_accepted_link), sizeof(int), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.hybrid_done_link), sizeof(int), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.story_sigma), sizeof(double), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.hybrid_accepted_sigma), sizeof(int), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.hybrid_done_sigma), sizeof(int), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.epsilon_pi_HMC), sizeof(double), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.story_pi), sizeof(double), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.hybrid_accepted_pi), sizeof(int), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }
    if(fread(&(simulation->info.hybrid_done_pi), sizeof(int), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }

    // compute current acceptance
    simulation->info.hybrid_done_acceptance_link = ((double) simulation->info.hybrid_accepted_link) / ((double) simulation->info.hybrid_done_link);
    simulation->info.hybrid_done_acceptance_sigma = ((double) simulation->info.hybrid_accepted_sigma) / ((double) simulation->info.hybrid_done_sigma);
    simulation->info.hybrid_done_acceptance_pi = ((double) simulation->info.hybrid_accepted_pi) / ((double) simulation->info.hybrid_done_pi);

    // read the measure number
    if(fread(&(simulation->info.measure), sizeof(int), 1, opf) != 1)
    {
        reading_errors++;
        // if there are troubles, try the copy file, else return
        if(try_copy == false)
        {
            fprintf(stderr, "Error when reading the file!\n");
            exit(1);
        }
    }

    // close the file
    fclose(opf);

    // if there was at least one reading error, read the copy file from the beginnning (now try_copy is enforced to false)
    if(reading_errors > 0) ReadBinaryLattice(simulation, simulation->terma_file_copy, false);
}