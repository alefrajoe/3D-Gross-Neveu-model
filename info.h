#ifndef INFO_H
#define INFO_H

#include "macro.h"
#include "random.h"

typedef struct Info
{
    // fundamenta parameters
    double beta;
    double beta_scalar;
    double mass;
    int sosia;
    int measure;
    // if metropolis Link is active
    #ifdef METROPOLIS_LINK
    double max_theta;
    int metropolis_link_accepted;
    int metropolis_link_done;
    double metropolis_link_acceptance;
    #endif

    #ifdef HYBRID_MC
    // link + sigma HMC params
    double epsilon_link_HMC;
    double distance_link_HMC;
    int steps_link_HMC;
    double epsilon_sigma_HMC;
    double distance_sigma_HMC;
    int steps_sigma_HMC;
    double epsilon_pi_HMC;
    double distance_pi_HMC;
    int steps_pi_HMC;

    // link HMC acceptance params
    int hybrid_accepted_link;
    int hybrid_done_link;
    double hybrid_done_acceptance_link;
    // sigma HMC acceptance params
    int hybrid_accepted_sigma;
    int hybrid_done_sigma;
    double hybrid_done_acceptance_sigma;
    int local_accepted_sigma;
    int local_done_sigma;
    double local_acceptance_sigma;
    // pi HMC acceptance params
    int hybrid_accepted_pi;
    int hybrid_done_pi;
    double hybrid_done_acceptance_pi;

    // delta_H_link saved for printing the acceptance file
    double delta_H_link;
    double delta_H_sigma;
    double delta_H_pi;

    // HMC story, these quantities are useful to regulate the HMC epsilon and steps
    int story_link;
    int story_sigma;
    int story_pi;
    #endif

}Info;

// Info functions
void StartInfo(struct Info * info, double beta, double beta_scalar, double mass, double epsilon_link_HMC, double epsilon_sigma_HMC, double epsilon_pi_HMC, int step_link_HMC, int step_sigma_HMC, int step_pi_HMC, int sosia);
void UpdateMetropolisAcceptance(struct Info * info, char * file_acceptance);
void UpdateHMCParameters(struct Info * info);
#endif