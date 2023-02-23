#ifndef RANDOM_C
#define RANDOM_C

#include "macro.h"
#include "random.h"


// this is the status of the random number generator
dsfmt_t rng_status;

///////////////////////////////////////////////////////////////////
//            uniform distribution int numbers
///////////////////////////////////////////////////////////////////
/**
 * Take a random integer in [0, extrema).
 */
void UniformInt(const int extrema, int *number)
{
  (*number) = rand() % extrema;
}

///////////////////////////////////////////////////////////////////
//           uniform distribution double numbers
///////////////////////////////////////////////////////////////////
/**
 * Generate a double number in (0,1).
 */
double Casuale(void)
{
  double ris;

  ris=dsfmt_genrand_open_open(&(rng_status));

  return ris;
}

/**
 * Generate a random double in (-extrema, +extrema) from a uniform distribution.
 */
void SymmetricCasuale(const double extrema, double *number)
{   
    (*number) = (2.0 * Casuale() - 1.0) * extrema;
}


////////////////////////////////////////////////////////////////
//                      seed
////////////////////////////////////////////////////////////////
/**
 * Initialize the random number generator.
 */
void Initrand(unsigned int s)
  {
    unsigned int seed;
    
  if(s==0)
    {
    seed=((unsigned int) time(NULL)+10) % UINT_MAX;
    }
  else
    {
    seed=s;
    }

  dsfmt_init_gen_rand(&(rng_status), seed);
  }

/**
 * Start the random seed inside the code.
 * The function throws a number N = seed of Casuale, so that different runs have different random seeds.
 */
void StartSeed(int seed)
{
    // initialize the random seed and throw a N="seed" number if Casuale()
    srand(time(NULL));
    for(int i=0; i<10+seed; i++) rand();
    
    myrand_init(24912, 2413);
    Initrand(rand());
}

/////////////////////////////////////////////////////////////
//                 gaussian distribution
/////////////////////////////////////////////////////////////
/**
 * Normal gaussian random number generator (polar method, knuth vol 2, p. 117).
 */
double Gauss1()
   {
   double v1, v2, s, ris;

   do
     {
     v1=1.0-2.0*Casuale();
     v2=1.0-2.0*Casuale();
     s=v1*v1+v2*v2;
     }
   while(s >= 1);

   ris=v1*sqrt(-2*log(s)/s);
   return ris;
   }


/**
 * Normal gaussian random number generator (polar method, knuth vol 2, p. 117).
 */
void Gauss2(double *ris1, double *ris2)
{
  double v1, v2, s;

  do
    {
    v1=-1.0+2.0*Casuale();
    v2=-1.0+2.0*Casuale();
    s=v1*v1+v2*v2;
    }
  while(s >= 1);

  *ris1=v1*sqrt(-2*log(s)/s);
  *ris2=v2*sqrt(-2*log(s)/s);
}


/**
 * Generate two random numbers distrubuted according to a Gaussian distribution
 * with average value "mean" and variance equal to "variance".
 */
void Gauss2MeanVariance(double *ris1, double *ris2, double mean, double variance)
{
  double v1, v2, s;

  do
    {
    v1=-1.0+2.0*Casuale();
    v2=-1.0+2.0*Casuale();
    s=v1*v1+v2*v2;
    }
  while(s >= 1);

  *ris1=mean + sqrt(variance) * v1*sqrt(-2*log(s)/s);
  *ris2=mean + sqrt(variance) * v2*sqrt(-2*log(s)/s);
}

#endif
