#ifndef RANDOM_H
#define RANDOM_H

#include "dSFMT.h"
#include "new_random.h"


// uniform distribution int numbers
void UniformInt(const int extrema, int *number);

// uniform distribution double numbers
double Casuale(void);
void SymmetricCasuale(const double extrema, double *number);

// seed
void StartSeed(int seed);
void Initrand(unsigned int s);

// gaussian distribution
double Gauss1();
void Gauss2(double *ris1, double *ris2);
void Gauss2MeanVariance(double *ris1, double *ris2, double mean, double variance);
#endif 