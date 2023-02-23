#include "link.h"

///////////////////////////////////////////////////////
//         general functions
///////////////////////////////////////////////////////
/**
 * The function prints a link to standard output.
 * It both prints out the element of the algebra and the group element.
 */
void PrintLink(struct Link *link)
{
    printf("Link phase is %.16f and link value is %.16f + i %.16f\n", carg(link->link), creal(link->link), cimag(link->link));
}

/**
 * Copy the link from into destination.
 */
void CopyLink(const struct Link const * from, struct Link * destination);

/**
 * Set the link to zero.
 */
void SetToZero(struct Link * link);

/**
 * The function initializes the link.
 * The way this is accomplished depends on the macro defined in "macro.h".
 * ----------------------------------------------------
 * START_LINK_IDENTITY : The phase of the link is set to zero, such that exp(I * phase) = 1.
 * START_LINK_RANDOM : The phase of the link is set to a random number in (-PI, PI).
 */
void StartLink(struct Link * link)
{
    #ifdef START_LINK_IDENTITY
    link->link = 1.0;
    #endif
    #ifdef START_LINK_RANDOM
    double phase = 0.0;
    SymmetricCasuale(PI, &phase);
    link->link = cexp(I_UNIT * Q_CHARGE * phase);
    #endif
}

/**
 * Correct numerical errors of the link, reshaping the link in order to leave it into 
 * the gauge group of definition.
 */
void CorrectLink(struct Link * link)
{
    double norm = pow(cabs(link->link), 0.5);
    
    // correct the link norm
    link->link = link->link / norm;
}

////////////////////////////////////////////////////////////////
//                Link operations
////////////////////////////////////////////////////////////////
/**
 * Save into staple the product link-link (it is equivalent to the sum of phases).
 */
void LinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * link3);
/**
 * Save into staple the product link-adjlink (it is equivalent to the difference of phases).
 */
void LinkAdjLink(const struct Link const * link1, const struct Link const * link2, struct Link * link3);
/**
 * Save into staple the product adjlink-link (it is equivalent to minus the difference of phases).
 */
void AdjLinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * link3);
/**
 * Save into staple the product adjlink-adjlink (it is equivalent to minus the sum of phases).
 */
void AdjLinkAdjLink(const struct Link const * link1, const struct Link const * link2, struct Link * link3);


/**
 * Add to staple the product link-link (it is equivalent to the sum of phases).
 */
void SumLinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple);
/**
 * Add to staple the product link-adjlink (it is equivalent to the difference of phases).
 */
void SumLinkAdjLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple);
/**
 * Add to staple the product adjlink-link (it is equivalent to minus the difference of phases).
 */
void SumAdjLinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple);
/**
 * Add to staple the product adjlink-adjlink (it is equivalent to minus the sum of phases).
 */
void SumAdjLinkAdjLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple);



/**
 * Subtract to staple the product link-link (it is equivalent to the sum of phases).
 */
void SubtractLinkLink(const struct Link const * link1, const struct Link const * link2, struct Link * staple);