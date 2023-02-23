//-------------------RANDOM NUMBER GENERATOR----------------------------
//
//
// Really minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
#ifndef NEW_RANDOM_H
#define NEW_RANDOM_H

#include<math.h>
#include<stdint.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng);

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq);

pcg32_random_t random_state; // random number internal state

void myrand_init(unsigned long int initstate, unsigned long int initseq);

double myrand(void);

#endif

//-------END OF RANDOM NUMBER GENERATOR----------------------------------