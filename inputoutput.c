#include "inputoutput.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                           input output funtions
//////////////////////////////////////////////////////////////////////////////////////////

/**
 * Check whether the input parameters are inserted correctly.
 * Print to stderr the errors if the number of parameters is different from expected.
 */
bool ReadCorrectInput(char **argv)
{
    int i=1;
    // read until NULL is find
    while(argv[i] != NULL) i++;

    // expected 8 parameters
    if(i == 11) return true;
    else
    {
        fprintf(stderr, "Expected 10 parameters, %d inserted:\n"
        "1) double: beta\n"
        "2) double: beta_scalar\n"
        "3) double: mass\n"
        "4) double: epsilon link HMC\n"
        "5) double: epsilon sigma HMC\n"
        "6) double: epsilon pi HMC\n"
        "7) int: step link HMC\n"
        "8) int: step sigma HMC\n"
        "9) int: step pi HMC\n"
        "10) int: seed number\n", i-1);
        exit(1);
    }
}

/**
 * Create all the directories required for the simulation of the model.
 */
bool CreateAllRequiredDirectories(void)
{
    // initialize strings
    char string_main_dir[100], string_data_simulation[300], string_acceptance_simulation[300], string_terma_simulation[300];
    struct stat stat_main_dir = {0}, stat_data_dir = {0}, stat_acceptance_dir = {0}, stat_terma_dir = {0};

    // sprintf all the strings
    sprintf(string_main_dir, "dati_simulation%dd", DIM);
    sprintf(string_data_simulation, "dati_simulation%dd/data_fermionsun1c%df%dl", DIM, NFLAV, LTIME);
    sprintf(string_acceptance_simulation, "dati_simulation%dd/acceptance_fermionsun1c%df%dl", DIM, NFLAV, LTIME);
    sprintf(string_terma_simulation, "dati_simulation%dd/terma_fermionsun1c%df%dl", DIM, NFLAV, LTIME);

    // if the directories do not exist, create them
    if (stat(string_main_dir, &stat_main_dir) == -1) mkdir(string_main_dir, 0777);
    if (stat(string_data_simulation, &stat_data_dir) == -1) mkdir(string_data_simulation, 0777);
    if (stat(string_acceptance_simulation, &stat_acceptance_dir) == -1) mkdir(string_acceptance_simulation, 0777);
    if (stat(string_terma_simulation, &stat_terma_dir) == -1) mkdir(string_terma_simulation, 0777);
}

/**
 * Create the string for the title of the obesrvables.
 */
void StringObservables(char * string_observables, const struct Info const * info)
{
    sprintf(string_observables, "dati_simulation%dd/data_fermionsun1c%df%dl/datafermion_sun1c%df%dlbeta%.8fbetascalar%.8fmass%.8ftemp%d.txt", DIM, NFLAV, LSPACE, NFLAV, LSPACE, info->beta, info->beta_scalar, info->mass, info->sosia);
}

/**
 * Create the string for the acceptance file.
 */
void StringAcceptance(char * string_acceptance, const struct Info const * info)
{
    #ifdef METROPOLIS_LINK
    sprintf(string_acceptance, "dati_simulation%dd/acceptance_fermionsun1c%df%dl/acceptancefermion_sun1c%df%dlbeta%.8fbetascalar%.8fmass%.8ftemp%d.txt", DIM, NFLAV, LSPACE, NFLAV, LSPACE, info->beta, info->beta_scalar, info->mass, info->sosia);
    #endif
    #ifdef HYBRID_MC
    sprintf(string_acceptance, "dati_simulation%dd/acceptance_fermionsun1c%df%dl/acceptancefermion_sun1c%df%dlbeta%.8fbetascalar%.8fmass%.8ftemp%d.txt", DIM, NFLAV, LSPACE, NFLAV, LSPACE, info->beta, info->beta_scalar, info->mass, info->sosia);
    #endif
}

/**
 * Create the string for the terma file.
 */
void StringTerma(char * string_terma, char * string_terma_copy, const struct Info const * info)
{
    #ifdef METROPOLIS_LINK
    sprintf(string_terma, "dati_simulation%dd/terma_fermionsun1c%df%dl/termafermion_sun1c%df%dlbeta%.8fbetascalar%.8fmass%.8ftemp%d.txt", DIM, NFLAV, LSPACE, NFLAV, LSPACE, info->beta, info->beta_scalar, info->mass, info->sosia);
    sprintf(string_terma_copy, "dati_simulation%dd/terma_fermionsun1c%df%dl/termafermion_sun1c%df%dlbeta%.8fbetascalar%.8fmass%.8ftemp%dcopy.txt", DIM, NFLAV, LSPACE, NFLAV, LSPACE, info->beta, info->beta_scalar, info->mass, info->sosia);
    #endif
    #ifdef HYBRID_MC
    sprintf(string_terma, "dati_simulation%dd/terma_fermionsun1c%df%dl/termafermion_sun1c%df%dlbeta%.8fbetascalar%.8fmass%.8ftemp%d.txt", DIM, NFLAV, LSPACE, NFLAV, LSPACE, info->beta, info->beta_scalar, info->mass, info->sosia);
    sprintf(string_terma_copy, "dati_simulation%dd/terma_fermionsun1c%df%dl/termafermion_sun1c%df%dlbeta%.8fbetascalar%.8fmass%.8ftemp%dcopy.txt", DIM, NFLAV, LSPACE, NFLAV, LSPACE, info->beta, info->beta_scalar, info->mass, info->sosia);    
    #endif
}

/**
 * Print the title (the first three lines) of the observables file.
 * These lines start with a hastag "#" symbol
 */
void OpenAndWriteObservablesTitle(FILE * opf, char * string_title)
{
    // try to open the file observable
    opf=fopen(string_title, "a");
    if(opf == NULL)
    {
        fprintf(stderr, "Error in opening observables data, cannot open %s\n", string_title);
        exit(1);
    }

    // first line
    // generalities : gauge group, DIM, L,..
    fprintf(opf, "# Gauge=U(1), NFLAV=%d, DIM=%d, LSPACE=%d, LTIME=%d, Q_CHARGE=%d\n", NFLAV, DIM, LSPACE, LTIME, Q_CHARGE);
    
    // second line
    // simulation parameters: IMEAS, LATTICE_UPDATE, MEASURE,..
    fprintf(opf, "# H=");
    #ifdef PURE_GAUGE
    fprintf(opf, "*PURE_GAUGE");
    #endif
    #ifdef PSEUDO_FERMION
    fprintf(opf, "*PSEUDO_FERMION");
    #endif
    #ifdef SIGMA_INTERACTION
    fprintf(opf, "*SIGMA_INTERACTION");
    #endif
    fprintf(opf, ", MEASURE=%d, IMEAS=%d, TERMALIZATION=%d, LATTICE_UPDATE=%d, LINK_INIT=", MEASURE, IMEAS, TERMALIZATION_STEP, LATTICE_UPDATE);
    #ifdef START_LINK_IDENTITY
    fprintf(opf, "LINK_ID");
    #endif
    #ifdef START_LINK_RANDOM
    fprintf(opf, "LINK_RANDOM");
    #endif

    #ifdef SIGMA_INTERACTION
    fprintf(opf, ", START_SIGMA=");
    #ifdef START_SIGMA_ZERO
    fprintf(opf, "START_SIGMA_ZERO");
    #endif
    #ifdef START_SIGMA_GAUSSIAN
    fprintf(opf, "START_SIGMA_GAUSSIAN");
    #endif
    #ifdef START_WEAK_SIGMA
    fprintf(opf, "START_WEAK_SIGMA=%.4f", START_WEAK_SIGMA);
    #endif
    #endif

    #ifdef METROPOLIS_LINK
    fprintf(opf, ", SIMULATION=*METROPOLIS_LINK");
    #endif
    #ifdef HYBRID_MC
    fprintf(opf, ", SIMULATION=*HYBRID_MC");
    #endif

    // HMC update ---------------------------
    #ifdef LF2SIGMA
    fprintf(opf, ", HMC_UPDATE=*LF2SIGMA");
    #endif
    #ifdef LFMN2SIGMA
    fprintf(opf, ", HMC_UPDATE=*LFMN2SIGMA");
    #endif
    #ifdef LFMN2SIGMAVF
    fprintf(opf, ", HMC_UPDATE=*LFMN2SIGMAVF");
    #endif
    #ifdef LF4MN4FPSIGMA
    fprintf(opf, ", HMC_UPDATE=*LF4MN4FPSIGMA");
    #endif

    // start x0 ---------------------------------
    #ifdef START_X0_NULL
    fprintf(opf, ", START_CG=*X0_NULL");
    #endif
    #ifdef START_X0_PREVIOUS
    fprintf(opf, ", START_CG=*X0_PREVIOUS");
    #endif

    // error CG --------------------------------
    fprintf(opf, ", MEDIUM_PREC=%E", MEDIUM_PRECISION_INVERSION);
    fprintf(opf, ", HIGH_PREC=%E", HIGH_PRECISION_INVERSION);

    // boundary conditions
    fprintf(opf, ", BC=");
    #ifndef APBC_FERMION_TIME
    #ifndef APBC_ALL
    fprintf(opf, "*PBC");
    #endif
    #endif
    #ifdef APBC_FERMION_TIME
    fprintf(opf, "*APBC_TIME");
    #endif
    #ifdef APBC_ALL
    fprintf(opf, "*APBC_ALL");
    #endif
    fprintf(opf, "\n");

    // third line
    // observables considered in the simulation
    fprintf(opf, "# measure,\tL,\tbeta,\tbeta_scalar,\tmass,\titeration,\tdistance");
    #ifdef OBS_MONOPOLE
    fprintf(opf, ",\tmonopole_density");
    #endif

    #ifdef OBS_CONDENSATE
    fprintf(opf, ",\tcondensate_tr1,\tcondensate_tr1_abs,\tcondensate_tr2,\tcondensate_tr4,\tcondensate_mu2,\tcondensate_mu4,\tcondensate_chi,\tcondensate_chi_abs,\tcondensate_chip,\tcondensate12_nosigma,\tcondensate12_withsigma");
    #endif

    #ifdef OBS_CONDENSATE_SIGMA
    fprintf(opf, ",\tcondensate_sigma,\tcondensate_sigma_abs,\tcondensate_sigma_mu2,\tcondensate_sigma_mu4,\tcondensate_sigma_chi,\tcondensate_sigma_chip,\tcondensate_sigma_chip_pbc1,\tcondensate_sigma_chip_pbc2");
    #endif
    
    #ifdef OBS_PLAQUETTE
    fprintf(opf, ",\tplaq1,\tplaq2,\tplaq3");
    #endif
    fprintf(opf, "\n");

    // close the file
    fclose(opf);
}

/**
 * Open and write the title of the acceptance file.
 */
void OpenAndWriteAcceptanceTitle(FILE * opf, char * string_title, const struct Info const * info)
{
    opf=fopen(string_title, "a");
    if(opf == NULL)
    {
        fprintf(stderr, "Error in opening acceptance file, cannot open %s\n", string_title);
        exit(1);
    }

    // first line
    // generalities : gauge group, DIM, L,..
    fprintf(opf, "# Gauge=U(1), NFLAV=%d, DIM=%d, LSPACE=%d, LTIME=%d, Q_CHARGE=%d\n", NFLAV, DIM, LSPACE, LTIME, Q_CHARGE);

    // second line
    // simulation parameters: IMEAS, LATTICE_UPDATE, MEASURE,..
    fprintf(opf, "# H=");
    #ifdef PURE_GAUGE
    fprintf(opf, "*PURE_GAUGE");
    #endif
    #ifdef PSEUDO_FERMION
    fprintf(opf, "*PSEUDO_FERMION");
    #endif
    #ifdef SIGMA_INTERACTION
    fprintf(opf, "*SIGMA_INTERACTION");
    #endif
    fprintf(opf, ", MEASURE=%d, IMEAS=%d, LATTICE_UPDATE=%d, TERMALIZATION=%d, LINK_INIT=", MEASURE, IMEAS, TERMALIZATION_STEP, LATTICE_UPDATE);
    #ifdef START_LINK_IDENTITY
    fprintf(opf, "LINK_ID, ");
    #endif
    #ifdef START_LINK_RANDOM
    fprintf(opf, "LINK_RANDOM, ");
    #endif

    #ifdef SIGMA_INTERACTION
    fprintf(opf, ", START_SIGMA=");
    #ifdef START_SIGMA_ZERO
    fprintf(opf, "START_SIGMA_ZERO");
    #endif
    #ifdef START_SIGMA_GAUSSIAN
    fprintf(opf, "START_SIGMA_GAUSSIAN");
    #endif
    #ifdef START_WEAK_SIGMA
    fprintf(opf, "START_WEAK_SIGMA=%.4f", START_WEAK_SIGMA);
    #endif
    #endif

    #ifdef METROPOLIS_LINK
    fprintf(opf, " UPDATE=METROPOLIS_LINK,\tMETROPOLIS_LINK_THRESHOLD=%d, METROPOLIS_LINK_ACCEPTANCE=%.8f", METROPOLIS_LINK_THRESHOLD, METROPOLIS_LINK_ACCEPTANCE);
    #endif
    #ifdef HYBRID_MC
    fprintf(opf, " HYBRID_EPSILON=%.16f, HYBRID_STEPS=%d, HYBRID_PATH=%.16f", info->epsilon_link_HMC, info->steps_link_HMC, info->distance_link_HMC);
    #endif
    // boundary conditions
    fprintf(opf, ", BC=");
    #ifndef APBC_FERMION_TIME
    #ifndef APBC_ALL
    fprintf(opf, "*PBC");
    #endif
    #endif
    #ifdef APBC_FERMION_TIME
    fprintf(opf, "*APBC_TIME");
    #endif
    #ifdef APBC_ALL
    fprintf(opf, "*APBC_ALL");
    #endif
    fprintf(opf, "\n");

    // third line
    // list the observables considered
    #ifdef METROPOLIS_LINK
    fprintf(opf, "#\n");
    #endif
    #ifdef HYBRID_MC
    fprintf(opf, "# beta, beta_scalar, mass, epsilon_link, step_link, distance_link, delta_H_link, acc_link, epsilon_sigma, step_sigma, distance_sigma, delta_H_sigma, acc_sigma\n");
    #endif

    // close the file
    fclose(opf);
}

//////////////////////////////////////////////////////////////////
//                   analysis function
//////////////////////////////////////////////////////////////////

/**
 * Return the mean value of the array that contains a num_elem of elements.
 */
double Mean(const double const * array, int num_elem)
{
    double mean_value = 0.0;

    // compute the mean
    for(int i=0; i<num_elem; i++) mean_value += array[i] / num_elem;

    return mean_value;
}