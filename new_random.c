//-------------------RANDOM NUMBER GENERATOR----------------------------
//
//
// Really minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

#include "new_random.h"

uint32_t pcg32_random_r(pcg32_random_t* rng)
    {
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
    {
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
    }

// ------------ my wrapper --------------

pcg32_random_t random_state; // random number internal state

void myrand_init(unsigned long int initstate, unsigned long int initseq)
  {
  pcg32_srandom_r(&random_state, initstate, initseq);
  }

double myrand(void)
  {
  return (double) pcg32_random_r(&random_state)/(pow(2.0, 32.0)-1.0);
  }

//-------END OF RANDOM NUMBER GENERATOR----------------------------------