#include "observables.h"

int main(int argc, char **argv)
{
    // check if the input is inserted correctly, otherwise exit(1)
    ReadCorrectInput(argv);
    CreateAllRequiredDirectories();

    // read all different simulation parameters
    double beta = atof(argv[1]);
    double beta_scalar = atof(argv[2]);
    double mass = atof(argv[3]);
    double epsilon_link_HMC = atof(argv[4]);
    double epsilon_sigma_HMC = atof(argv[5]);
    double epsilon_pi_HMC = atof(argv[6]);
    int step_link_HMC = atoi(argv[7]);
    int step_sigma_HMC = atoi(argv[8]);
    int step_pi_HMC = atoi(argv[9]);
    int sosia = atoi(argv[10])+10;

    // allocate variables -------------------------------------------------------
    struct Simulation simulation;
    AllocateSimulation(&simulation);
    struct Observables obs;
    AllocateObservables(&obs);

    // start simulation ---------------------------------------------------------
    StartSimulation(&simulation, beta, beta_scalar, mass, epsilon_link_HMC, epsilon_sigma_HMC, epsilon_pi_HMC, step_link_HMC, step_sigma_HMC, step_pi_HMC, sosia);
    // if the proper file exist in TERMA directory, read the link of the lattice and the measure number from the appropriate terma file
    ReadBinaryLattice(&simulation, simulation.terma_file, true);

    // termalization if no measure has been taken (first run) -------------------------------
    //if(simulation.info.measure == 0){UpdateLatticeTermalization(&simulation); exit(0);}

    // take measures ------------------------------------------------------------
    TakeMeasures(&simulation, &obs);

    // destroy simulation -------------------------------------------------------
    DestroySimulation(&simulation);
    DestroyObservables(&obs);
    return 0;
}
