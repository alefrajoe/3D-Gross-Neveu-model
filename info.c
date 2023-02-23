#include "info.h"

/**
 * Initialize the Info parameters.
 * Note that many of these variables may or may not be defined dependending on the presence or absence of the related
 * macro in "macro.h".
 * The parameter beta (the temperature) is passed to the StartInfo function.
 * The function also initializes the random seed.
 */
void StartInfo(struct Info * info, double beta, double beta_scalar, double mass, double epsilon_link_HMC, double epsilon_sigma_HMC, double epsilon_pi_HMC,  int step_link_HMC, int step_sigma_HMC, int step_pi_HMC, int sosia)
{
    info->beta = beta;
    info->beta_scalar = beta_scalar;
    info->mass = mass;
    info->sosia = sosia;
    info->measure = 0;
    #ifdef METROPOLIS_LINK
    info->max_theta = PI;
    info->metropolis_link_accepted = 0;
    info->metropolis_link_done = 0;
    info->metropolis_link_acceptance = 0.0;
    #endif

    #ifdef HYBRID_MC

    info->distance_link_HMC = epsilon_link_HMC * step_link_HMC;
    // chech possible distance errors
    if(info->distance_link_HMC < 0.0) exit(1);
    info->epsilon_link_HMC = epsilon_link_HMC;
    info->steps_link_HMC = step_link_HMC;
    info->hybrid_accepted_link = 0;
    info->hybrid_done_link = 0;
    info->hybrid_done_acceptance_link = 0.0;
    #ifdef SIGMA_INTERACTION
    info->distance_sigma_HMC = epsilon_sigma_HMC * step_sigma_HMC;
    // check possible distance errors
    if(info->distance_sigma_HMC < 0.0) exit(1);
    info->epsilon_sigma_HMC = epsilon_sigma_HMC;
    info->steps_sigma_HMC = step_sigma_HMC;
    info->hybrid_accepted_sigma = 0;
    info->hybrid_done_sigma = 0;
    info->hybrid_done_acceptance_sigma = 0.0;
    info->local_accepted_sigma = 0;
    info->local_done_sigma = 0;
    info->local_acceptance_sigma = 0.0;
    #endif
    #ifdef PI_INTERACTION
    info->distance_pi_HMC = step_pi_HMC * epsilon_pi_HMC;
    if(info->distance_pi_HMC < 0.0) exit(1);
    info->epsilon_pi_HMC = epsilon_pi_HMC;
    info->steps_pi_HMC = step_pi_HMC;
    info->hybrid_accepted_pi = 0;
    info->hybrid_done_pi = 0;
    info->hybrid_done_acceptance_pi = 0.0;
    #endif
    info->delta_H_link = 0.0;
    info->delta_H_sigma = 0.0;
    info->delta_H_pi = 0.0;

    info->story_link = 0;
    info->story_sigma = 0;
    info->story_pi = 0;
    #endif

    // check, if step is equal to one exit
    if(info->steps_sigma_HMC == 1) {printf("Cannot simulate 1 step for sigma HMC!!!\n"); exit(1);};

    // initialize random seed
    StartSeed(info->sosia);
}

/**
 * Update the info related to the metropolis update of the link variables.
 * In particular, the variables info->metropolis_link_accepted, info->metropolis_link_done and info->metropolis_link_acceptance
 * are set to zero.
 * The info are updated only if the number of Metropolis done is larger than METROPOLIS_LINK_THRESHOLD
 */ 
void UpdateMetropolisAcceptance(struct Info * info, char * file_acceptance)
{
    #ifdef METROPOLIS_LINK
    if(info->metropolis_link_done > METROPOLIS_LINK_THRESHOLD)
    {
        // compute the acceptance
        info->metropolis_link_acceptance = (double) info->metropolis_link_accepted / (double) info->metropolis_link_done;
        // compute the ration between the current acceptance and the one expected
        double ratio = info->metropolis_link_acceptance / METROPOLIS_LINK_ACCEPTANCE;


        // if the current acceptance is larger than the wanted one
        if(ratio >= 1)
        {
            // if ratio * max_theta is smaller than PI modify the max_theta, else set to PI directly
            if(ratio * info->max_theta <= PI) info->max_theta = (info->max_theta * ratio);
            else info->max_theta = PI;
        }
        else
        {
            // if ratio is zero (it should not happen) by default the range is reduced by a factor 10
            // otherwise, max_theta is multiplied by the ratio
            if(ratio == 0)  info->max_theta = (info->max_theta / 10);
            else info->max_theta = (info->max_theta * ratio);
        }

    #ifdef PRINT_LINK_ACCEPTANCE
    FILE * opf;
    opf=fopen(file_acceptance, "a");
    fprintf(opf, "Metropolis link acceptance is %d/%d=%.16f, max angle in PI units is %.16f\n", info->metropolis_link_accepted, info->metropolis_link_done, info->metropolis_link_acceptance, info->max_theta/PI);
    fclose(opf);
    #endif

    // reinitialize the parameters
    info->metropolis_link_accepted = 0;
    info->metropolis_link_done = 0;
    }
    #endif
}

/**
 * The function update the HMC parameters depending on their stories (so their trend).
 * In particular, if a story has 25 positive event, multiply the corresponding step by 1.1
 * If a story has 2 negative events, reduce the epsilon parameters by a factor 1.2
 */
void UpdateHMCParameters(struct Info * info)
{
    #ifdef AUTOREGULATION    
    // ----------- after the first 100 updates, every 100 update refine epsilon
    if(info->local_done_sigma == 100)
    {
        if(info->local_acceptance_sigma >= 0.85) {info->epsilon_sigma_HMC = info->epsilon_sigma_HMC * 1.15; info->distance_sigma_HMC = info->steps_sigma_HMC * info->epsilon_sigma_HMC; info->local_accepted_sigma = 0; info->local_done_sigma = 0; info->local_acceptance_sigma = 0.0;}
        else if(info->local_acceptance_sigma < 0.85 && info->local_acceptance_sigma >= 0.8) {info->epsilon_sigma_HMC = info->epsilon_sigma_HMC * 1.08; info->distance_sigma_HMC = info->steps_sigma_HMC * info->epsilon_sigma_HMC; info->local_accepted_sigma = 0; info->local_done_sigma = 0; info->local_acceptance_sigma = 0.0;}
        else if(info->local_acceptance_sigma < 0.8 && info->local_acceptance_sigma >= 0.75) {info->epsilon_sigma_HMC = info->epsilon_sigma_HMC / 1.08; info->distance_sigma_HMC = info->steps_sigma_HMC * info->epsilon_sigma_HMC; info->local_accepted_sigma = 0; info->local_done_sigma = 0; info->local_acceptance_sigma = 0.0;}
        else if(info->local_acceptance_sigma < 0.75) {info->epsilon_sigma_HMC = info->epsilon_sigma_HMC / 1.15; info->distance_sigma_HMC = info->steps_sigma_HMC * info->epsilon_sigma_HMC; info->local_accepted_sigma = 0; info->local_done_sigma = 0; info->local_acceptance_sigma = 0.0;}
    }
    #endif
}