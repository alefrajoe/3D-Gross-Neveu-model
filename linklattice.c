#include "linklattice.h"


//////////////////////////////////////////////////////////////////
//                  Link energy
//////////////////////////////////////////////////////////////////
/**
 * Compute the real part of the plaquette associated with the site given in the direction mu-nu.
 * The site is passed in fermion Lessicographic coordinates
 * The function is named as "real trace plaquette", but in the case of U(1) it only computes the real part of the 
 * plaquette (the trace is trivial for a 1 x 1 matrix).
 * ======================================================
 * Return : Real { \Pi_{ \mu , \nu }( site ) }
 * ======================================================
 */
double RealTracePlaquette(const struct Lattice const * lat, int site, int mu, int nu)
{
    // temp double where the value of the plaquette will be stored
    struct Link temp1, temp2;
    
    #ifdef DEBUG_LINK
    if(mu == nu) printf("In real trace plaquette two equal directions!!!\n");
    #endif

    // compute the plaquette
    LinkLink(&(lat->dlink[COOR(site, mu)]), &(lat->dlink[COOR(lat->pp[COOR(site, mu)], nu)]), &temp1);
    LinkAdjLink(&temp1, &(lat->dlink[COOR(lat->pp[COOR(site, nu)], mu)]), &temp2);
    LinkAdjLink(&temp2, &(lat->dlink[COOR(site, nu)]), &temp1);

    return creal(temp1.link);
}

/**
 * Compute the imaginary part of the plaquette associated with the site given (in fermionic coordinates) in the direction mu-nu.
 * The function is named as "imaginary trace plaquette", but in the case of U(1) it only computes the imaginary part of the 
 * plaquette (the trace is trivial for a 1 x 1 matrix).
 * ======================================================
 * Return : Imag { \Pi_{ \mu , \nu }( site ) }
 * ======================================================
 */
double ImagTracePlaquette(const struct Lattice const * lat, int site, int mu, int nu)
{
    // temp double where the value of the plaquette will be stored
    struct Link temp1, temp2;

    // compute the plaquette
    LinkLink(&(lat->dlink[COOR(site, mu)]), &(lat->dlink[COOR(lat->pp[COOR(site, mu)], nu)]), &temp1);
    LinkAdjLink(&temp1, &(lat->dlink[COOR(lat->pp[COOR(site, nu)], mu)]), &temp2);
    LinkAdjLink(&temp2, &(lat->dlink[COOR(site, nu)]), &temp1);

    return cimag(temp1.link);
}

/**
 * Sum over the volume the real part of all the plaquettes associated with the lattice.
 * Of course, as well as in the case of the computation of a single plaquette, the trace is trivial in the U(1) case.
 * ========================================================================
 * Return : \sum_{ x } \sum_{\mu > \nu}  Real { \Pi_{ \mu , \nu }( x ) }
 * ========================================================================
 */
double TotalRealTracePlaquette(const struct Lattice const * lat)
{
    double temp = 0.0;

    // for all sites in the lattice
    for(int x=0; x<VOLUME; x++)
    {
        // for all directions \mu
        for(int mu=0; mu<DIM; mu++)
        {
            // for \mu > \nu
            for(int nu=0; nu<mu; nu++)
            temp += RealTracePlaquette(lat, x, mu, nu);
        }
    }

    return temp;
}

/**
 * Sum over the volume the imaginary part of all the plaquettes associated with the lattice.
 * Of course, as well as in the case of the computation of a single plaquette, the trace is trivial in the U(1) case.
 * ========================================================================
 * Return : \sum_{ x } \sum_{\mu > \nu}  imag { \Pi_{ \mu , \nu }( x ) }
 * ========================================================================
 */
double TotalImagTracePlaquette(const struct Lattice const * lat)
{
    double temp = 0.0;

    // for all sites in the lattice
    for(int x=0; x<VOLUME; x++)
    {
        // for all directions \mu
        for(int mu=0; mu<DIM; mu++)
        {
            // for \mu > \nu
            for(int nu=0; nu<mu; nu++)
            temp += ImagTracePlaquette(lat, x, mu, nu);
        }
    }

    return temp;   
}

/**                                                                                               * _________ < __________*
 *                                                                                                |                       |
 *                                                                                                v                       ^ 
 * Compute the staple associated with the link U_{ \mu }( site ).                                 |         U_\mu(x)      |
 *                                                                    =  \sum_{ \nu \neq \mu}     *  -  -  -  >  -  -  -  *
 * The staple variable is initialized to zero into the function.                                  |                       |
 * The staple evaluated is saved into the Link structure passed to the function.                  ^                       v
 *                                                                                                |                       |
 *                                                                                                *__________<____________*
 */
void ComputeStaple(const struct Lattice const * lat, int site, int mu, struct Link * staple)
{
    struct Link temp1;

    // initialize the staple container
    staple->link = 0.0;

    // for all directions different from mu
    for(int nu=0; nu<DIM; nu++)
    if(nu != mu)
    {
        // compute forward staple
        LinkAdjLink(&(lat->dlink[COOR(lat->pp[COOR(site, mu)], nu)]), &(lat->dlink[COOR(lat->pp[COOR(site, nu)], mu)]), &temp1);
        SumLinkAdjLink(&temp1, &(lat->dlink[COOR(site, nu)]), staple);

        // backward staple
        AdjLinkAdjLink(&(lat->dlink[COOR(lat->pp[COOR(lat->mm[COOR(site, nu)], mu)], nu)]), &(lat->dlink[COOR(lat->mm[COOR(site, nu)], mu)]), &temp1);
        SumLinkLink(&temp1, &(lat->dlink[COOR(lat->mm[COOR(site, nu)], nu)]), staple);
    }
}

/**
 * Compute the sum of the imaginary part of all the plaquettes that contain the link
 * passed to the function (the link is indirectly given through the site and the direction).
 * The site is passed in Lessicographic coordinates.
 * The resulting sum is returned by the funtion.
 * =============================================================================
 * Return : \sum_{ \nu \neq \mu }     Imag{ \Pi_{ \mu \nu } ( x )}
 * =============================================================================
 */
double ComputeImagAllPlaquettesIncludingLink(const struct Lattice const * lat, int site, int mu)
{
    struct Link temp1, temp2;
    double total = 0.0;

    // for all directions different from mu
    for(int nu=0; nu<DIM; nu++)
    if(nu != mu)
    {

        // compute forward plaquettes
        LinkLink(&(lat->dlink[COOR(site, mu)]), &(lat->dlink[COOR(lat->pp[COOR(site, mu)], nu)]), &temp1);
        LinkAdjLink(&temp1, &(lat->dlink[COOR(lat->pp[COOR(site, nu)], mu)]), &temp2);
        LinkAdjLink(&temp2, &(lat->dlink[COOR(site, nu)]), &temp1);
        // pay attention to the signs with the imaginary part!!!
        total -= cimag(temp1.link);

        // backward staple
        LinkLink(&(lat->dlink[COOR(lat->mm[COOR(site, nu)], mu)]), &(lat->dlink[COOR(lat->pp[COOR(lat->mm[COOR(site, nu)], mu)], nu)]), &temp1);
        LinkAdjLink(&temp1, &(lat->dlink[COOR(site, mu)]), &temp2);
        LinkAdjLink(&temp2, &(lat->dlink[COOR(lat->mm[COOR(site, nu)], nu)]), &temp1);
        // pay attention to the signs with the imaginary parts!!!
        total += cimag(temp1.link);
    }

    return total;
}

/**
 * Compute the real part of the trace between a link and a staple (both of them are passed to the function).
 * Note that differently from the function that we are going to define in the following, the staple may not be associated
 * with the staple passed to the function.
 * ===========================================================================
 * Return :  Real { Trace { L_\mu (x) * Staple } } 
 * ===========================================================================
 */
double RealTraceLinkStaple(const struct Link const * link, const struct Link const * staple);

/**
 * Compute the imaginary part of the trace between a link and a staple (both of them are passed to the function).
 * Note that differently from the function that we are going to define in the following, the staple may not be associated
 * with the staple passed to the function.
 * ===========================================================================
 * Return :  Imag { Trace { L_\mu (x) * Staple } } 
 * ===========================================================================
 */
double ImagTraceLinkStaple(const struct Link const * link, const struct Link const * staple);


/**
 * The function returns the real trace of a link (located at site along the mu direction) 
 * multiplied by the staple associated with it.
 * This is equivalent to consider the sum of the real trace of the plaquettes that contain U_\mu(x).
 * ===========================================================================
 * Return : \sum_{ \Pi that contains U_\mu(x) } Real { \Pi_{\mu \nu} (x) } 
 * ===========================================================================
 */
double RealTraceLinkSiteStapleSameSite(const struct Lattice const * lat, int site, int mu)
{
    // initialize temp variable and define the staple
    double temp = 0.0; struct Link staple;

    // compute the staple
    ComputeStaple(lat, site, mu, &staple);
    return creal(lat->dlink[COOR(site, mu)].link * staple.link);
}

/**
 * This function saves in a double structure or array the expectation values of the gauge repart of the action.
 * In particular, the plaqutte density and its power-raised quantities are stored into the pointers passed to the function.
 */
void ObsPlaquette(const struct Lattice const * lat, double * plaq1, double * plaq2, double * plaq3)
{
    // compute the plaquette density (apart from a factor DIM(DIM-1)/2.0)
    double plaq_density = TotalRealTracePlaquette(lat) / VOLUME;

    (*plaq1) = plaq_density;
    (*plaq2) = pow(plaq_density, 2);
    (*plaq3) = pow(plaq_density, 3);
}

/**
 * Copy all links from lattice to a copy.
 */
void CopyAllLinkFromLattice(const struct Lattice * lat, complex * copy)
{
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    copy[COOR(x, d)] = lat->dlink[COOR(x, d)].link;
}

/**
 * Copy all link from a copy to the lattice.
 */
void CopyAllLinkToLattice(const complex * copy, struct Lattice * lat)
{
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    lat->dlink[COOR(x, d)].link = copy[COOR(x, d)];
}

// ///////////////////////////////////////////////////////////////////////////////////////
// //                              Metropolis Link
// ///////////////////////////////////////////////////////////////////////////////////////
// /**
//  * Construct a simple Metropolis algorithm for U(1) pure gauge theories (it is used as a test).
//  * A trial phase is picked from the domain (-max_theta, max_theta), as described into the info structure.
//  */
void MetropolisLink(struct Lattice * lat, struct Info * info, int site)
{
    #ifdef METROPOLIS_LINK
    // for all lattice directions
    for(int d=0; d<DIM; d++)
    {
        #ifdef DEBUG_LINK
        double start_energy = (- TotalRealTracePlaquette(lat));
        #endif
        // trial link for the metropolis test
        struct Link trial;

        // pick a random angle for the trial matrix
        double theta = 0.0;
        SymmetricCasuale(info->max_theta, &theta);
    
        if(Casuale() < 0.5) trial.link = lat->dlink[COOR(site, d)].link * cexp(I_UNIT * Q_CHARGE * theta); 
        else trial.link = lat->dlink[COOR(site, d)].link * cexp(- I_UNIT * Q_CHARGE * theta);

        // evaluate the staple in to compute successively the energy variation
        Link staple;
        ComputeStaple(lat, site, d, &staple);

        // compute the energy variation
        double energy = 0.0, trial_energy = 0.0;

        
        energy -= RealTraceLinkStaple(&(lat->dlink[COOR(site, d)]), &staple);
        trial_energy -= RealTraceLinkStaple(&trial, &staple); 


        // metropolis test
        // it the accept-reject test is accepted
        if(Casuale() <= exp( - info->beta * (trial_energy - energy) ) )
        {
            // update link
            lat->dlink[COOR(site, d)].link = trial.link;
            // update link info
            info->metropolis_link_accepted++;
            info->metropolis_link_done++;

            #ifdef DEBUG_LINK
            double final_energy = (- TotalRealTracePlaquette(lat));
            if(fabs((trial_energy - energy) - (final_energy - start_energy)) > 1e-13) printf("In metropolis link the energy difference is %.16f\n", (trial_energy - energy) - (final_energy - start_energy));
            #endif
        }
        else
        {
            // else it is only done
            info->metropolis_link_done++;
        }
    }
    #endif
}

/**
 * Correct all the link in the lattice.
 * The link is reshaped in order to have unit norm, so that it is a U(1) phase.
 */
void CorrectLinkLattice(struct Lattice * lat)
{
    for(int x=0; x<VOLUME; x++)
    for(int d=0; d<DIM; d++)
    CorrectLink(&(lat->dlink[COOR(x, d)]));
}